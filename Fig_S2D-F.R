# Prep data for running scUTRquant ----

library(data.table)
library(dplyr)

samp = "E4_2yr"

barcodes = fread(paste0(samp,"_cells.csv"), skip = 1)
clusters = fread(paste0(samp,"_clusters.csv"), skip = 1, header = F)
embeddings = fread(paste0(samp,"_embeddings.csv"), skip = 1)
# these files were generated in Fig 1 prep.

annos_e3 = data.table(cell_id = paste0(samp,"_",barcodes$V2), cluster = clusters$V2, umap_1 = embeddings$V2, umap_2 = embeddings$V3)
annos_e4 = data.table(cell_id = paste0(samp,"_",barcodes$V2), cluster = clusters$V2, umap_1 = embeddings$V2, umap_2 = embeddings$V3)

annos = rbind(annos_e3, annos_e4)
fwrite(annos, "annos.csv")

# Analysis after run ----

library(magrittr)
library(tidyverse)
library(cowplot)
library(SingleCellExperiment)
library(scutrboot)
library(scater)

set.seed(20220627)
N_BOOTSTRAPS=10000
MIN_NCELLS=50

sce = readRDS("E3E4_2yr.txs.Rds") # scUTRquant output
sce = logNormCounts(sce) 
sce %<>% `[`(,!is.na(colData(.)$cluster))
sce %<>% `[`(!rowData(.)$is_blacklisted,)
idx_multi <- rowData(sce)$atlas.utr_type == 'multi'
idx_consistent <- rowData(sce)$is_consistent
tim_clusters = c("TIMs","Effector-hi TIMs","Serpine1+ TIMs","Early TIMs")
idx_tim = colData(sce)$cluster %in% tim_clusters
gene_id2sym <- rowData(sce) %>% as_tibble %$% setNames(gene_symbol, gene_id)
plot_pval_histogram <- function (df) {
  ggplot(df, aes(x=pval)) +
    geom_histogram(boundary=0, bins=41, 
                   color='black', fill='lightgrey', size=0.2) +
    scale_x_continuous(expand=c(0.02,0.02)) +
    scale_y_continuous(expand=c(0,0,0,100)) +
    labs(x="p-values", y="Genes") +
    theme_bw()
}


#to compare within tims, make sce[idx_multi,] into sce[idx_multi,idx_tim]. same for hom
#WD test
df_test_wd <- testTwoSample(sce[idx_multi,idx_tim], assayName="logcounts",
                            sampleKey="sample_id", geneKey="gene_id",
                            sample0="E3_2yr", sample1="E4_2yr", 
                            statistic="WD", nBootstraps=N_BOOTSTRAPS, 
                            minCellsPerGene=MIN_NCELLS, 
                            featureExclude="is_ipa",) %>%
  as_tibble() %>%
  mutate(qval=p.adjust(pval, method="BH"))
# what about diff between TIM and effector-hi?
df_test_wd <- testTwoSample(sce[idx_multi,], assayName="logcounts",
                            sampleKey="cluster", geneKey="gene_id",
                            sample0="Effector-hi TIMs", sample1="TIMs", 
                            statistic="WD", nBootstraps=N_BOOTSTRAPS, 
                            minCellsPerGene=MIN_NCELLS, 
                            featureExclude="is_ipa",) %>%
  as_tibble() %>%
  mutate(qval=p.adjust(pval, method="BH"))

plot_pval_histogram(df_test_wd)

#LUI test
#note: dLUI = sample1 - sample0
df_test_lui <- testTwoSample(sce[idx_multi & idx_consistent,idx_tim], assayName="logcounts",
                             sampleKey="sample_id", geneKey="gene_id",
                             sample0="E3_2yr", sample1="E4_2yr", 
                             statistic="UI", featureIndex="is_distal",
                             nBootstraps=N_BOOTSTRAPS, 
                             minCellsPerGene=MIN_NCELLS, 
                             featureExclude="is_ipa") %>%
  as_tibble() %>%
  mutate(qval=p.adjust(pval, method="BH"))
# same, TIM vs effector-hi
df_test_lui <- testTwoSample(sce[idx_multi & idx_consistent,], assayName="logcounts",
                             sampleKey="cluster", geneKey="gene_id",
                             sample0="Effector-hi TIMs", sample1="TIMs", 
                             statistic="UI", featureIndex="is_distal",
                             nBootstraps=N_BOOTSTRAPS, 
                             minCellsPerGene=MIN_NCELLS, 
                             featureExclude="is_ipa") %>%
  as_tibble() %>%
  mutate(qval=p.adjust(pval, method="BH"))
plot_pval_histogram(df_test_lui)

#IPA test
idx_ipa <- rowData(sce) %>%
  as_tibble %>%
  group_by(gene_id) %>%
  mutate(has_ipa=any(is_ipa), has_nipa=any(!is_ipa)) %>%
  ungroup() %>%
  select(transcript_id, is_ipa, has_ipa, has_nipa) %>%
  filter(has_ipa, has_nipa) %$%
  transcript_id

df_test_ipa <- testTwoSample(sce[idx_ipa,idx_tim], assayName="logcounts",
                             sampleKey="sample_id", geneKey="gene_id",
                             sample0="E3_2yr", sample1="E4_2yr", 
                             statistic="UI", featureIndex="is_ipa",
                             nBootstraps=N_BOOTSTRAPS, 
                             minCellsPerGene=MIN_NCELLS) %>%
  as_tibble() %>%
  mutate(qval=p.adjust(pval, method="BH"))
#again, TIM vs effector-hi
df_test_ipa <- testTwoSample(sce[idx_ipa,], assayName="logcounts",
                             sampleKey="cluster", geneKey="gene_id",
                             sample0="Effector-hi TIMs", sample1="TIMs", 
                             statistic="UI", featureIndex="is_ipa",
                             nBootstraps=N_BOOTSTRAPS, 
                             minCellsPerGene=MIN_NCELLS) %>%
  as_tibble() %>%
  mutate(qval=p.adjust(pval, method="BH"))
plot_pval_histogram(df_test_ipa)

# visualize
library(ggrastr)
options(ggrepel.max.overlaps = Inf)

MIN_DIPA=0.15
MIN_DLUI=0.15
MIN_WD=0.15
MAX_QVAL=0.05

luiVolcanoPlot <- function (df, max_q=0.05, min_dlui=0.1, max_dlui=0.5, legend=TRUE, dpi=300) {
  possible_labels = apply(expand.grid(c(sprintf("|dLUI| < %0.2f", min_dlui),
                                        sprintf("dLUI < -%0.2f", min_dlui),
                                        sprintf("dLUI > %0.2f", min_dlui)),
                                      c("NS",sprintf("q < %0.2f", max_q))),
                          1,paste,collapse=", ")
  
  df %<>% 
    mutate(sig=ifelse(qval < max_q, sprintf("q < %0.2f", max_q), "NS")) %>%
    mutate(size_effect=ifelse(abs(stat) < min_dlui, 
                              sprintf("|dLUI| < %0.2f", min_dlui),
                              ifelse(stat < 0, 
                                     sprintf("dLUI < -%0.2f", min_dlui),
                                     sprintf("dLUI > %0.2f", min_dlui)))) %>%
    mutate(sig_size=paste(size_effect, sig, sep=", ")) %>%
    mutate(outlier=abs(stat) > max_dlui,
           stat_adj=ifelse(outlier, sign(stat)*max_dlui, stat)) %>% 
    mutate(sig_size = factor(sig_size, levels = possible_labels))
  
  
  
  pmain <- df %>% 
    ggplot(aes(x=stat_adj, y=-log10(pval), label = ifelse(qval < max_q & abs(stat) > min_dlui, gene_id2sym[gene], ''))) +
    rasterize(geom_point(aes(color=sig_size, shape=outlier, size=sig_size)), dpi=dpi) +
    geom_vline(xintercept=min_dlui*c(-1,1), 
               linetype='dashed', size=0.5, color="darkgrey") +
    scale_color_manual(labels = possible_labels,
                       values=c("grey", "darkgrey", "darkgrey", 
                                "#958DBE", "#3F6CA5", "#C63732"),
                       drop = FALSE) +
    scale_x_continuous(breaks=scales::pretty_breaks(6), 
                       limits=max_dlui*c(-1,1)) +
    scale_shape_manual(values=c("FALSE"=16, "TRUE"=17)) +
    scale_size_manual(labels = possible_labels, values=c(0.3, 0.3, 0.3, 0.3, 1.2, 1.2), drop = FALSE) +
    ggrepel::geom_label_repel(size = 3.5, force = 15, box.padding = 0.5, point.padding = 0.5, segment.size = 0.1, min.segment.length = 0.1) + 
    labs(x="dLUI", 
         y="-log10(p-value)", color="Significance\n (BH Adjusted)") +
    theme(legend.title.align=0.5) + 
    guides(shape="none") +
    theme_cowplot()
  
  if (!legend) {
    pmain <- pmain + guides(color="none", size="none")
  }
  
  xdens <- axis_canvas(pmain, axis='x') + 
    geom_density(data=df, aes(x=stat, fill=sig), alpha=0.7, size=0.2) +
    scale_fill_manual(values=c("black", "#958DBE"))
  
  insert_xaxis_grob(pmain, xdens, grid::unit(0.2, 'null'), position='top') %>% 
    ggdraw()
}

ipaVolcanoPlot <- function (df, max_q=0.05, min_dipa=0.1, max_dipa=0.5, legend=TRUE,dpi=300) {
  possible_labels = apply(expand.grid(c(sprintf("|dIPA| < %0.2f", min_dipa),
                                        sprintf("dIPA < -%0.2f", min_dipa),
                                        sprintf("dIPA > %0.2f", min_dipa)),
                                      c("NS",sprintf("q < %0.2f", max_q))),
                          1,paste,collapse=", ")
  
  df %<>% 
    mutate(sig=ifelse(qval < max_q, sprintf("q < %0.2f", max_q), "NS")) %>%
    mutate(size_effect=ifelse(abs(stat) < min_dipa, 
                              sprintf("|dIPA| < %0.2f", min_dipa),
                              ifelse(stat < 0, 
                                     sprintf("dIPA < -%0.2f", min_dipa),
                                     sprintf("dIPA > %0.2f", min_dipa)))) %>%
    mutate(sig_size=paste(size_effect, sig, sep=", ")) %>%
    mutate(outlier=abs(stat) > max_dipa,
           stat_adj=ifelse(outlier, sign(stat)*max_dipa, stat)) %>%
    mutate(sig_size = factor(sig_size, levels = possible_labels))
  
  pmain <- df %>% 
    ggplot(aes(x=stat_adj, y=-log10(pval), label = ifelse(qval < max_q & abs(stat) > min_dipa, gene_id2sym[gene], ''))) +
    rasterize(geom_point(aes(color=sig_size, shape=outlier, size=sig_size)), dpi=dpi) +
    geom_vline(xintercept=min_dipa*c(-1,1), 
               linetype='dashed', size=0.5, color="darkgrey") +
    scale_color_manual(labels = possible_labels,
                       values=c("grey", "darkgrey", "darkgrey", 
                                "#958DBE", "#3F6CA5", "#C63732"),
                       drop = FALSE) +
    scale_shape_manual(values=c("FALSE"=16, "TRUE"=17)) +
    scale_size_manual(labels = possible_labels, values=c(0.3, 0.3, 0.3, 0.3, 1.2, 1.2), drop = FALSE) +
    scale_x_continuous(breaks=scales::pretty_breaks(6), 
                       limits=max_dipa*c(-1,1)) +
    ggrepel::geom_label_repel(size = 3.5, force = 15, box.padding = 0.5, point.padding = 0.5, segment.size = 0.1, min.segment.length = 0.1) + 
    labs(x="dIPA", 
         y="-log10(p-value)", color="Significance\n (BH Adjusted)") +
    guides(shape="none") +
    theme(legend.title.align=0.5) + 
    theme_cowplot()
  
  if (!legend) {
    pmain <- pmain + guides(color="none", size="none")
  }
  
  xdens <- axis_canvas(pmain, axis='x') + 
    geom_density(data=df, aes(x=stat, fill=sig), alpha=0.7, size=0.2) +
    scale_fill_manual(values=c("black", "red"))
  
  insert_xaxis_grob(pmain, xdens, grid::unit(0.2, 'null'), position='top') %>% 
    ggdraw()
}

wdVolcanoPlot <- function (df, max_q=0.05, min_wd=0.1, max_wd=0.5, legend=TRUE,dpi=300) {
  possible_labels = apply(expand.grid(c(sprintf("WD < %0.2f", min_wd),
                                        sprintf("WD > %0.2f", min_wd)),
                                      c("NS",sprintf("q < %0.2f", max_q))),
                          1,paste,collapse=", ")
  
  df %<>% 
    mutate(sig=ifelse(qval < max_q, sprintf("q < %0.2f", max_q), "NS")) %>%
    mutate(size_effect=ifelse(stat > min_wd, 
                              sprintf("WD > %0.2f", min_wd), 
                              sprintf("WD < %0.2f", min_wd))) %>%
    mutate(sig_size=paste(size_effect, sig, sep=", ")) %>%
    mutate(outlier=stat > max_wd,
           stat_adj=ifelse(outlier, max_wd, stat)) %>%
    mutate(sig_size = factor(sig_size, levels = possible_labels))
  
  
  pmain <- df %>% 
    ggplot(aes(x=stat_adj, y=-log10(pval), label = ifelse(qval < max_q & stat > min_wd, gene_id2sym[gene], ''))) +
    rasterize(geom_point(aes(color=sig_size), size=0.7), dpi=dpi) +
    geom_vline(xintercept=min_wd, linetype='dashed', size=0.5, color="darkgrey") +
    scale_color_manual(labels = possible_labels, values=c("grey", "darkgrey", "#663333", "red"), drop = FALSE) +
    scale_shape_manual(values=c("FALSE"=16, "TRUE"=17)) +
    scale_x_continuous(breaks=scales::pretty_breaks(6), 
                       limits=c(0,max_wd)) +
    ggrepel::geom_label_repel(size = 3.5, force = 15, box.padding = 0.5, point.padding = 0.5, segment.size = 0.1, min.segment.length = 0.1) + 
    labs(x="Total Isoform Change", 
         y="-log10(p-value)", color="Significance\n (BH Adjusted)") +
    guides(shape="none") +
    theme(legend.title.align=0.5) + 
    theme_cowplot()
  
  if (!legend) {
    pmain <- pmain + guides(color="none")
  }
  
  xdens <- axis_canvas(pmain, axis='x') + 
    geom_density(data=df, aes(x=stat, fill=sig), alpha=0.7, size=0.2) +
    scale_fill_manual(values=c("black", "red"))
  
  insert_xaxis_grob(pmain, xdens, grid::unit(0.2, 'null'), position='top') %>% 
    ggdraw()
}

#set dpi to 300 while testing, set to 600 before saving
## Fig. S2D ----
wdVolcanoPlot(df_test_wd, max_q=0.2, min_wd=0.1, legend = FALSE, max_wd = 0.35, dpi = 600)
## Fig. S2E ----
luiVolcanoPlot(df_test_lui, max_q=0.25, min_dlui=0.15, legend = FALSE, max_dlui = 0.3, dpi = 600)
## Fig. S2F ----
ipaVolcanoPlot(df_test_ipa, max_q=0.5, min_dipa=0.05, legend = FALSE, max_dipa = 0.2, dpi = 600)
