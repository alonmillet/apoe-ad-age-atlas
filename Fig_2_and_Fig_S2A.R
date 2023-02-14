library(data.table)
library(Seurat)
library(dplyr)
library(ComplexHeatmap)

# SCENIC ----

mglia = readRDS("mglia_only.rds") # structure generated in Fig. 1 code
dat = fread("scenic_aucell_output.csv") # generated in Python; see accompanying code
dat = as.data.frame(dat, rownames = TRUE)
colnames(dat) = gsub("\\s*\\([^\\)]+\\)","",as.character(colnames(dat)))
dat$ident = as.character(mglia$mglia_ident)

dat_dt = dat %>% as.data.table
dat_dt = dat_dt[,-1] #drop cell column
dat_dt_z = dat_dt[, lapply(.SD,scale), .SDcols = names(dat_dt)[1:(ncol(dat_dt)-1)]] #convert to z-scores
dat_dt_z$ident = dat_dt$ident
colnames(dat_dt_z) = colnames(dat_dt)
dat_dt_toplot = dat_dt_z[, lapply(.SD,mean), by = ident]
thresh = dat_dt[, lapply(.SD, sd), .SDcols = names(dat_dt)[1:(ncol(dat_dt)-1)]] %>% unlist %>% unname
rangediff = function(x){return(max(x)-min(x))}
thresh_max = dat_dt_toplot[, lapply(.SD, rangediff), .SDcols = names(dat_dt_z)[1:(ncol(dat_dt)-1)]] %>% unlist %>% unname
agg_toplot = as.matrix(dat_dt_toplot[,-1])
agg_toplot = agg_toplot[, which(thresh>0.025 & thresh_max > 1.75)]
rownames(agg_toplot) = dat_dt_toplot$ident
o1 = seriation::seriate(dist(agg_toplot), method = "GW")
o2 = seriation::seriate(dist(t(agg_toplot)), method = "GW")
ap1_fam = c("Atf1","Atf2","Atf3","Atf4","Atf5","Atf6","Atf6b","Batf","Batf3",
            "Fos","Fosb","Fosl1","Jun","Junb","Jund","Nrf1","Mafb","Maff","Mafg","Mafk")
irf_fam = c("Irf1","Irf2","Irf3","Irf4","Irf5","Irf7","Irf8","Irf9","Nfkb1","Nfkb2","Stat1","Stat2","Stat3","Stat4")
fams = ifelse(colnames(agg_toplot) %in% ap1_fam, "AP-1", ifelse(colnames(agg_toplot) %in% irf_fam, "IFN", "None"))
ha = HeatmapAnnotation(`TF Family` = fams, col = list(`TF Family` = c("AP-1" = "chartreuse2", "IFN" = "deeppink", "None" = "lightgray")),
                       annotation_legend_param = list(`TF Family` = list(title = "TF Family", title_position = "leftcenter-rot")))
ht = Heatmap(agg_toplot, name = "Aggregated\nModules", bottom_annotation = ha, cluster_rows = as.dendrogram(o1[[1]]), column_split = 6, row_split = 8,
             rect_gp = gpar(col = "black", lwd = 0.25), cluster_columns = as.dendrogram(o2[[1]]), column_gap = unit(5,"mm"), row_gap = unit(2,"mm"),
             row_title = c(), column_names_rot = 45, column_title = c("Homeostatic","Inflammatory","Interferon\nAssociated","Proliferation\nAssociated","Cd74+ Microglia\nAssociated","Il34+\nAssociated"),
             heatmap_legend_param = list(title = "Aggregated Transcription Factor Z-Scores", title_position = "leftcenter-rot", legend_height = unit(8,"cm")))
png("hmap.png",width=30,height=12,units="in",res=600)
draw(ht, heatmap_legend_side = "left", annotation_legend_side = "left", padding = unit(c(2,2,10,15),"mm"))
dev.off()

## Fig. 2A ----
agg_toplot = as.matrix(dat_dt_toplot[,-1])
agg_toplot = agg_toplot[, which(thresh>0.03 & thresh_max > 2)]
rownames(agg_toplot) = dat_dt_toplot$ident
o1 = seriation::seriate(dist(agg_toplot), method = "GW")
o2 = seriation::seriate(dist(t(agg_toplot)), method = "GW")
ap1_fam = c("Atf1","Atf2","Atf3","Atf4","Atf5","Atf6","Atf6b","Batf","Batf3",
            "Fos","Fosb","Fosl1","Jun","Junb","Jund","Nrf1","Mafb","Maff","Mafg","Mafk")
irf_fam = c("Irf1","Irf2","Irf3","Irf4","Irf5","Irf7","Irf8","Irf9","Nfkb1","Nfkb2","Stat1","Stat2","Stat3","Stat4")
fams = ifelse(colnames(agg_toplot) %in% ap1_fam, "AP-1", ifelse(colnames(agg_toplot) %in% irf_fam, "IFN", "None"))
ht = Heatmap(agg_toplot, name = "Aggregated\nModules", cluster_rows = as.dendrogram(o1[[1]]), column_split = 6, row_split = 8,
             rect_gp = gpar(col = "black", lwd = 0.25), cluster_columns = as.dendrogram(o2[[1]]), column_gap = unit(5,"mm"), row_gap = unit(2,"mm"),
             row_title = c(), column_names_rot = 45, column_title = c("Homeostatic","Il34+\nAssociated","Inflammatory","Prolif.","DAM-2","IFN"),
             heatmap_legend_param = list(title = "Aggregated Transcription Factor Z-Scores", title_position = "leftcenter-rot", legend_height = unit(8,"cm")))
png("hmap2.png",width=16,height=6,units="in",res=600)
draw(ht, heatmap_legend_side = "left", annotation_legend_side = "left", padding = unit(c(2,2,10,15),"mm"))
dev.off()

## Fig. 2B ----
agg_toplot = as.matrix(dat_dt_toplot[,-1])
agg_toplot = agg_toplot[, which(thresh>0.01 & thresh_max > 0.5)]
rownames(agg_toplot) = dat_dt_toplot$ident
agg_toplot_filt = agg_toplot[rownames(agg_toplot) %in% c("Effector-hi TIMs","TIMs","Serpine1+ TIMs","Homeostatic Microglia",
                                                         "Lars2-mid Homeostatic Microglia","Lars2-hi Homeostatic Microglia"),]
agg_toplot_filt = agg_toplot_filt[,colnames(agg_toplot) %in% c(ap1_fam,irf_fam)]
o1 = seriation::seriate(dist(agg_toplot_filt), method = "GW")
o2 = seriation::seriate(dist(t(agg_toplot_filt)), method = "GW")
fams = ifelse(colnames(agg_toplot_filt) %in% ap1_fam, "AP-1", ifelse(colnames(agg_toplot_filt) %in% irf_fam, "IFN", "None"))
ha = HeatmapAnnotation(`TF Family` = fams, col = list(`TF Family` = c("AP-1" = "chartreuse2", "IFN" = "deeppink", "None" = "lightgray")),
                       annotation_legend_param = list(`TF Family` = list(title = "TF Family", title_position = "topcenter", nrow = 1)))
ht = Heatmap(agg_toplot_filt, name = "Aggregated\nModules", bottom_annotation = ha, cluster_rows = as.dendrogram(o1[[1]]), column_split = 3, row_split = 3,
             rect_gp = gpar(col = "black", lwd = 0.25), cluster_columns = as.dendrogram(o2[[1]]), column_gap = unit(5,"mm"), row_gap = unit(2,"mm"),
             row_title = c(), column_title = c(),
             heatmap_legend_param = list(title = "Aggregated Transcription Factor Z-Scores", title_position = "leftcenter-rot", legend_height = unit(8,"cm")))
png("hmap_ap1_ifn_only.png",width=10,height=4,units="in",res=600)
draw(ht, heatmap_legend_side = "left", annotation_legend_side = "bottom", padding = unit(c(2,2,2,2),"mm"))
dev.off()

# CellPhoneDB ----

## Prep files for CellPhoneDB ----
library(Seurat)
library(SeuratObject)
library(dplyr)
library(data.table)
library(org.Hs.eg.db)

# prep counts matrix with human orthologs for cellphonedb
adapoe = readRDS("adapoe.rds")
counts = as.matrix(GetAssayData(adapoe, slot = "counts", assay = "RNA"))
counts = apply(counts, 2, function(x) (x/sum(x))*10000) #recommended cpdb normalization

# convert mouse genes to human genes
orthologs = fread("HMD_HumanPhenotype.rpt.txt") # download from MGI
genelist = data.table(mouse = rownames(counts))
genelist$index = 1:nrow(genelist)
genelist = merge.data.table(genelist, orthologs, by.x = "mouse", by.y = "V3")
genelist = genelist[order(genelist$index),]
ensembs = mapIds(org.Hs.eg.db, keys = genelist$V2 %>% as.character, keytype = "ENTREZID", column = "ENSEMBL")
genelist$ensembs = ensembs
counts = counts[genelist$index, ]
rownames(counts) = genelist$ensembs
counts_dt = as.data.table(counts, keep.rownames = TRUE)
fwrite(counts_dt, "cpdb_counts_matrix.txt", sep = "\t", 
       row.names = FALSE, col.names = TRUE, quote = FALSE)

# also prep metadata file
Idents(adapoe) = adapoe$subclust_id
md = Idents(adapoe) %>% as.data.table(keep.rownames = TRUE)
colnames(md) = c("Cell","cell_type")
fwrite(md, "cpdb_metadata.txt", sep = "\t", 
       row.names = FALSE, col.names = TRUE, quote = FALSE)

## Plot after generating data ----
library(data.table)
library(Seurat)
library(dplyr)
library(circlize)
library(ggplot2)
library(ggridges)
library(forcats)
library(ggrepel)
library(tidyverse)
library(stringr)
library(pheatmap)

adapoe = readRDS("adapoe.rds")
cluster_order = levels(adapoe$subclust_id)

##rename clusters and pairs
pvals = fread("pvalues.txt") # CellPhoneDB output
means = fread("means.txt") # CellPhoneDB output
setnames(pvals, colnames(pvals[,12:ncol(pvals)]), 
         sapply(colnames(pvals[,12:ncol(pvals)]), function(x) str_replace(x, fixed("|"), paste0(" ",sprintf("\U1F812")," "))) %>% as.character) #these are arrows
setnames(means, colnames(means[,12:ncol(means)]), 
         sapply(colnames(means[,12:ncol(means)]), function(x) str_replace(x, fixed("|"), paste0(" ",sprintf("\U1F812")," "))) %>% as.character) # these are arrows
pvals$interacting_pair = sapply(pvals$interacting_pair, function(x) str_replace(x, fixed("_"), paste0(" ",sprintf("\U00D7")," "))) # these are x signs
means$interacting_pair = sapply(means$interacting_pair, function(x) str_replace(x, fixed("_"), paste0(" ",sprintf("\U00D7")," "))) # these are x signs

## heatmap plot (code from CellPhoneDB package)
heatmaps_plot = function(meta_file, pvalues_file, count_filename, log_filename, count_network_filename, interaction_count_filename, 
                         count_network_separator, interaction_count_separator, show_rownames = T, show_colnames = T,
                         scale="none", cluster_cols = F,border_color='white', cluster_rows = F, fontsize_row=11,
                         fontsize_col = 11, main = '',treeheight_row=0, family='Arial', treeheight_col = 0,
                         col1 = "dodgerblue4", col2 = 'peachpuff', col3 = 'deeppink4', meta_sep='\t', pvalues_sep='\t', pvalue=0.05){
  #######   Network
  
  meta = read.csv(meta_file, comment.char = '', sep=meta_sep)
  
  all_intr = read.table(pvalues_file, header=T, stringsAsFactors = F, sep=pvalues_sep, comment.char = '', check.names = F)
  intr_pairs = all_intr$interacting_pair
  all_intr = all_intr[,-c(1:11)]
  
  
  split_sep = '\\|'
  join_sep = '|'
  
  pairs1_all = unique(meta[,2])
  
  pairs1 = c()
  for (i in 1:length(pairs1_all))
    for (j in 1:length(pairs1_all))
      pairs1 = c(pairs1,paste(pairs1_all[i],pairs1_all[j],sep=join_sep))
  
  all_count = matrix(ncol=3)
  colnames(all_count) = c('SOURCE','TARGET','count')
  count1 = c()
  for(i in 1:length(pairs1))
  {
    p1 = strsplit(pairs1[i], split_sep)[[1]][1]
    p2 = strsplit(pairs1[i], split_sep)[[1]][2]
    
    n1 = intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]
    pairs_rev = paste(p2, p1, sep=join_sep)
    n2 = intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]
    
    if(p1!=p2)
      count1 = length(unique(n1))+length(unique(n2))
    else
      count1 = length(unique(n1))
    
    new_count = c(p1,p2,count1)
    names(new_count) = c('SOURCE','TARGET','count')
    all_count = rbind(all_count, new_count)
  }
  
  all_count = all_count[-1,]
  write.table(all_count, count_network_filename, sep=count_network_separator, quote=F, row.names = F)
  
  #######   count interactions
  
  count1 = c()
  for(i in 1:length(pairs1))
  {
    p1 = strsplit(pairs1[i], split_sep)[[1]][1]
    p2 = strsplit(pairs1[i], split_sep)[[1]][2]
    
    n1 = intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]
    
    pairs_rev = paste(p2, p1, sep=join_sep)
    n2 = intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]
    if(p1!=p2)
      count1 = c(count1,length(unique(n1))+length(unique(n2)))
    else
      count1 = c(count1,length(unique(n1)))
    
  }
  
  if (any(count1)>0)
  {
    count_matrix = matrix(count1, nrow=length(unique(meta[,2])), ncol=length(unique(meta[,2])))
    rownames(count_matrix)= unique(meta[,2])
    colnames(count_matrix)= unique(meta[,2])
    
    all_sum = rowSums(count_matrix)
    all_sum = cbind(names(all_sum), all_sum)
    write.table(all_sum, file=interaction_count_filename, quote=F, sep=count_network_separator, row.names=F)
    
    col.heatmap <- colorRampPalette(c(col1,col2,col3 ))( 1000 )
    
    count_matrix = count_matrix[match(cluster_order,rownames(count_matrix)),][,match(cluster_order,colnames(count_matrix))]
    
    pheatmap(count_matrix, show_rownames = show_rownames, show_colnames = show_colnames, scale=scale, cluster_cols = cluster_cols,
             border_color=border_color, cluster_rows = cluster_rows, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
             main = main, treeheight_row = treeheight_row, family = family,color = col.heatmap, treeheight_col = treeheight_col, filename = count_filename)
    
    pheatmap(log(count_matrix+1), show_rownames = show_rownames, show_colnames = show_colnames, scale=scale, cluster_cols = cluster_cols,
             border_color=border_color, cluster_rows = cluster_rows, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
             main = main, treeheight_row = treeheight_row, family = family,color = col.heatmap, treeheight_col = treeheight_col, filename = log_filename)
  } else {
    stop("There are no significant results using p-value of: ", pvalue, call.=FALSE)
  }
}

heatmaps_plot(meta_file = "cpdb_metadata.txt", 
              pvalues_file = "pvalues.txt",
              count_filename = "heatmap_count.png",
              log_filename = "heatmap_log_count.png",
              count_network_filename = "network.txt",
              interaction_count_filename = "interaction_count.txt",
              count_network_separator = "\t",
              interaction_count_separator = "\t")






### Fig. 2C ----
dat = fread("deconvoluted.txt")
dat = subset(dat, !duplicated(dat[,complex_name]))
dat[complex_name == "", complex_name := protein_name]
dat = melt.data.table(dat, measure.vars = colnames(dat)[7:(length(cluster_order)+6)], variable.name = "cluster", value.name = "deconv") %>%
  setorder(-deconv)
dat$cluster = factor(dat$cluster, levels = levels(adapoe$subclust_id))
dat$rank = seq_len(nrow(dat))
dat[, stat := (deconv)/mean(deconv), by = cluster]
dat[complex_name == "", complex_name := paste0(gene_name,"_signaling")]
dat.mglia = dat[cluster %in% levels(adapoe)[1:21]]

library(ggtext)
dat.mglia[cluster %in% levels(adapoe)[1:5], type := "Homeostatic"]
dat.mglia[cluster %in% levels(adapoe)[16:18], type := "TIMs"]
dat.mglia[cluster %in% cluster_order[9:10], type := "DAMs"]

comp = dat.mglia[, mean(stat), by = .(complex_name,type)]
comp = comp[!is.na(type)]
comp$V1 = log(1.25+comp$V1) #log-norm to get rid of infs/nas
comp = dcast(comp, complex_name ~ type, value.var = "V1")
highlight.these.high = c("IL11_receptor","Dehydroepiandrosterone_bySTS","2arachidonoylglycerol_byDAGLB")
highlight.these.low = c("integrin_aMb2_complex","LeukotrieneB4_byLTA4H","TGFbeta_receptor1")

comp$score = comp$TIMs - comp$Homeostatic
comp$label = ifelse(comp$score > 0.15 | comp$score < -0.37, comp$complex_name, "")

comp %>% 
  mutate(facet = ifelse(score >= 0, "Up in TIMs", "Up in Homeostatic")) %>% 
  mutate(facet = factor(facet, levels = c("Up in TIMs", "Up in Homeostatic"))) %>%
  arrange(desc(score)) %>% 
  filter((row_number() > max(row_number()) - 15 & score<0) | (row_number() <= 15 & score>0)) %>%
  mutate(complex_name = paste("<span style = 'color: ",
                              ifelse(complex_name %in% highlight.these.low,"blue",
                                     ifelse(complex_name %in% highlight.these.high, "red","black")),
                              ";'>",
                              ifelse(complex_name %in% c(highlight.these.low,highlight.these.high), "<b>",""),
                              complex_name,
                              ifelse(complex_name %in% c(highlight.these.low,highlight.these.high), "</b>",""),
                              "</span>",sep = "")) %>%
  mutate(complex_name = factor(complex_name, levels = c(rev(complex_name[1:15]),complex_name[16:30]))) %>%
  {ggplot(.,aes(x = score, y = complex_name)) + 
      geom_segment(aes(y=complex_name,yend=complex_name,x=0,xend=score), lwd=1.5) +
      geom_point(size=5) + 
      labs(x = "TIM Score - Homeostatic Score", y = "") +
      theme_bw() +
      theme(axis.text.y = element_markdown(angle = 0, hjust = 1, size = 22.5),
            axis.text.x = element_text(size=17.5),
            axis.title = element_text(size=25),
            strip.text = element_text(size=25)) +
      facet_wrap(~facet,scales="free")} %>% 
  ggsave("tim_hom_complexes_lollipop.png",.,width=15.748,height=6.299,dpi=600)

fwrite(comp, "tim_hom_comp.csv") # for comparing to 1yr multiome in Fig. 3

### Fig. 2D ----
comp$score = comp$TIMs - comp$DAMs
comp$label = ifelse(comp$score > 0.2 | comp$score < -0.2, comp$complex_name, "")

comp %>% 
  mutate(facet = ifelse(score >= 0, "Up in TIMs", "Up in DAMs")) %>% 
  mutate(facet = factor(facet, levels = c("Up in TIMs", "Up in DAMs"))) %>%
  arrange(desc(score)) %>% 
  filter((row_number() > max(row_number()) - 15 & score<0) | (row_number() <= 15 & score>0)) %>%
  mutate(complex_name = paste("<span style = 'color: ",
                              ifelse(complex_name %in% highlight.these.low,"blue",
                                     ifelse(complex_name %in% highlight.these.high, "red","black")),
                              ";'>",
                              ifelse(complex_name %in% c(highlight.these.low,highlight.these.high), "<b>",""),
                              complex_name,
                              ifelse(complex_name %in% c(highlight.these.low,highlight.these.high), "</b>",""),
                              "</span>",sep = "")) %>%
  mutate(complex_name = factor(complex_name, levels = c(rev(complex_name[1:15]),complex_name[16:30]))) %>%
  {ggplot(.,aes(x = score, y = complex_name)) + 
      geom_segment(aes(y=complex_name,yend=complex_name,x=0,xend=score), lwd=1.5) +
      geom_point(size=5) + 
      labs(x = "TIM Score - DAM Score", y = "") +
      theme_bw() + 
      theme(axis.text.y = element_markdown(angle = 0, hjust = 1, size = 22.5),
            axis.text.x = element_text(size=17.5),
            axis.title = element_text(size=25),
            strip.text = element_text(size=25)) +
      facet_wrap(~facet,scales="free")} %>% 
  ggsave("tim_dam_complexes_lollipop.png",.,width=15.748,height=6.299,dpi=600)

fwrite(comp, "tim_dam_comp.csv") # for comparing to 1yr multiome in Fig. 3

### Fig. 2E ----
dat = fread("network.txt") # output from running heatmaps_plot function
dat = dcast(dat, SOURCE ~ TARGET, value.var = "count")
dat.m = as.matrix(dat[,-1])
rownames(dat.m) = dat$SOURCE
dat.m = dat.m[match(cluster_order,rownames(dat.m)),][,match(cluster_order,colnames(dat.m))]

dat.m_filt = dat.m
for (col in seq_len(ncol(dat.m))){
  cutoff = mean(unname(quantile(dat.m[,col]))[4:5])
  for (row in seq_len(nrow(dat.m))) {
    if (dat.m[row,col] < cutoff){
      dat.m[row,col] = 0
    }
  }
}
group = structure(c(rep("HMs",5),rep("Other Ms",2),"HMs",rep("IMs",7),rep("TIMs",3),rep("Other Ms",3),rep("Other Immune",16)),
                  names = rownames(dat.m))
grid.col = structure(c(rep(2,5),rep(3,2),2,rep(4,7),rep(5,3),rep(3,3),rep(6,16)),
                     names = rownames(dat.m))
circos.clear()
chordDiagram(dat.m_filt, group = group, grid.col = grid.col, directional = 1,
             direction.type = c("arrows","diffHeight"), link.arr.type = "big.arrow",
             diffHeight = mm_h(3), target.prop.height = mm_h(2))

#with rotated labels
png("circos.png", height = 20, width = 20, units = "in", res = 600)
circos.par(circle.margin = c(0.25,0.5,0.65,0.65))
chordDiagram(dat.m_filt, group = group, grid.col = grid.col, directional = 1,
             direction.type = c("arrows","diffHeight"), link.arr.type = "big.arrow",
             diffHeight = mm_h(3), target.prop.height = mm_h(2), 
             annotationTrack = "grid", preAllocateTracks = 1.1)
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 2)
}, bg.border = NA)
dev.off()

### Fig. S2A ----
group_barplot = structure(c(rep("Homeostatic",5),rep("Other Microglia",2),"Homeostatic",rep("Inflammatory",7),
                            rep("TIMs",3),rep("Other Microglia",3),rep("Non-Microglia",16)),
                          names = rownames(dat.m))
dat_sums = rowSums(dat.m)
dat_sums = data.table(cluster = factor(names(dat_sums),levels = cluster_order), ints = dat_sums, group = group_barplot)
f <- function(x) {list(mean(x), sd(x)/sqrt(length(x)))} 
dat_sums[,c("mean","sem") := f(ints),by=group]
dat_sums = unique(dat_sums, by = c("group","mean","sem"))
dat_sums$group = factor(dat_sums$group, levels = c("Homeostatic","Inflammatory","Other Microglia","TIMs","Non-Microglia"))
(ggplot(dat_sums, aes(x=group,y=mean,fill=group)) + 
    geom_bar(stat = "identity", color = "black", linewidth = 1.5) + 
    geom_errorbar(aes(ymin=mean-sem,ymax=mean+sem),width=0.2, linewidth=1.5,position=position_dodge(0.9)) +
    scale_fill_manual(values = c("red","cornflowerblue","mediumpurple1","orange","chartreuse")) +
    theme_bw() + 
    labs(y = "Total Interactions", x= "", fill = "Group") +
    theme(axis.title = element_text(size=30), axis.text.y = element_text(size=25), 
          axis.text.x = element_text(angle = 45, hjust = 1, size = 25),
          legend.title = element_text(size=30), legend.text = element_text(size=25), 
          legend.spacing = unit(0.1,"in")) +
    guides(fill = guide_legend(byrow = TRUE))
) %>% ggsave("interactions_barplot.png",.,height=16,width=12,dpi=600)

# Compass ----
## Prep data for running Compass ----
library(Seurat)
library(Signac)
library(dplyr)
library(data.table)

mglia = readRDS("mglia_only.rds")

md = mglia@meta.data %>% as.data.table
md$samp = 1:nrow(md)
md$agg_id = md$mglia_ident
md = split(md, md$agg_id)
for(i in 1:length(md)){
  num_to_pseudobulk = ceiling(nrow(md[[i]])/50) # pseudobulk into 50 cells per "cell"
  md[[i]]$pseudobulk_id = paste0(md[[i]]$agg_id,"_",rep(1:num_to_pseudobulk))
}
md = do.call("rbind",md)
setorder(md, "samp")
mglia$pseudobulk_id = md$pseudobulk_id
out = AverageExpression(mglia, group.by = "pseudobulk_id", assays = "RNA", slot = "counts")$RNA
write.table(out, file = "dat_matrix.tsv", sep = "\t", quote = FALSE) # this is the Compass input

type = unname(colnames(out))
type = sapply(strsplit(type, "_"), function(x) paste(x[1], collapse = "-"))
metadata = data.frame("type" = type)
rownames(metadata) = unname(colnames(out))
write.table(metadata, file = "cell_metadata.csv", sep = "\t", quote = FALSE)

## Plot after Compass run ----
library(data.table)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(matrixStats)
library(stringr)
library(Seurat)

mglia = readRDS("mglia_only.rds")
dat = fread("reactions.tsv") # output from Compass
penalties = as.matrix(dat[,-1])
rownames(penalties) = dat$V1
cell_metadata = fread("cell_metadata.csv") # generated above
cell_metadata$sample = paste0(str_split(cell_metadata$V1,"_") %>% sapply(.,"[[",1),"_",str_split(cell_metadata$V1,"_") %>% sapply(.,"[[",2))
cell_metadata$type = str_split(cell_metadata$V1,"_") %>% sapply(.,"[[",3)
cell_metadata$joined = paste0(str_split(cell_metadata$V1,"_") %>% sapply(.,"[[",1),"_",cell_metadata$type)
cell_metadata$genotype = str_split(cell_metadata$V1,"_") %>% sapply(.,"[[",1)
grps = cell_metadata$joined
rxn_metadata = fread("compass_reaction_metadata.csv") # accompanies Compass documentation

reaction.consistencies = -log1p(penalties)
reaction.consistencies = reaction.consistencies[which(apply(reaction.consistencies, 1, max) - apply(reaction.consistencies, 1, min) > 1e-3),]
reaction.consistencies = reaction.consistencies - min(reaction.consistencies)
reaction.consistencies = scale(reaction.consistencies)

rc.means = t(reaction.consistencies) %>% as.data.table(keep.rownames = TRUE)
rc.means$rn = cell_metadata$joined
rc.means = rc.means[, lapply(.SD,mean), by = rn]
rc.means.mat = as.matrix(rc.means[,-1])
rownames(rc.means.mat) = rc.means$rn
metrics = abs(colSds(rc.means.mat)/colMeans(rc.means.mat))
filt.mean.mat = rc.means.mat[,which(metrics>1)]
# rename reactions
rxns = strsplit(colnames(filt.mean.mat),"_(?!.*_)",perl = TRUE) %>% sapply("[[",1) %>% as.data.table
rxns$dir = strsplit(colnames(filt.mean.mat),"_(?!.*_)",perl = TRUE) %>% sapply("[[",2)
rxns$dir = ifelse(rxns$dir == "pos", "(+)","(-)")
rxns = merge(rxns, rxn_metadata, by.x = ".", by.y = "reaction_no_direction")
rxns[, new_name := paste0(stringr::str_trunc(reaction_name,43)," ",dir)]
colnames(filt.mean.mat) = rxns$new_name
# prep for annotations
col_fun = randomcoloR::distinctColorPalette(length(unique(rxns$subsystem)))
names(col_fun) = unique(rxns$subsystem)
confs = ifelse(rxns$confidence %in% c(0,4), "High", "Low")
templist = c()
for(i in 1:length(levels(mglia))){
  templist[3*i-2] = levels(mglia)[i]
  templist[3*i-1] = levels(mglia)[i]
  templist[3*i] = levels(mglia)[i]
}
row_order = paste0(rep(c("E2_","E3_","E4_"),21),templist)
row_order = row_order[row_order %in% rownames(filt.mean.mat)]

o2 = seriation::seriate(dist(t(filt.mean.mat)), method = "GW")
ht2 = Heatmap(filt.mean.mat, name = "COMPASS\nReaction\nScores",
              row_order = row_order,
              rect_gp = gpar(col = "black", lwd = 0.25), cluster_columns = as.dendrogram(o2[[1]]), row_gap = unit(2,"mm"),
              row_title = c(), column_title = c(), column_names_max_height = unit(35,"cm"),
              bottom_annotation = HeatmapAnnotation(Subsystem = rxns$subsystem, Confidence = confs, 
                                                    col = list(Subsystem = col_fun, Confidence = c("High" = "lightblue","Low" = "grey")), 
                                                    annotation_name_side = "left"),
              top_annotation = HeatmapAnnotation(Density = anno_boxplot(filt.mean.mat, height = unit(2, "cm")), 
                                                 annotation_name_side = "left", annotation_name_rot = 0))
png("hmap_means_by_genoandclust.png",width=20,height=15,units="in",res=600)
draw(ht2)
dev.off()

## Fig. 2F ----
library(ggblend)
library(ggtext)
rc = t(reaction.consistencies) %>% as.data.table(keep.rownames = TRUE)
rc$rn = strsplit(cell_metadata$sample,"_") %>% sapply(.,"[[",1)
means = rc[rn %in% c("E3","E4"), lapply(.SD,mean), by = rn][,-1] %>% as.matrix
mean_diffs = colDiffs(means) %>% as.vector # this is E4 - E3, as we want
sds = (rc[rn %in% c("E3","E4"), lapply(.SD, sd), .SDcols = colnames(rc)[-1]][1] %>% as.matrix)[1,]
cohends = mean_diffs/sds
cohends = data.table(rxn = names(cohends), cohend = cohends)
cohends$rxnrename = strsplit(cohends$rxn,"_(?!.*_)",perl = TRUE) %>% sapply("[[",1) %>% as.data.table
cohends$dir = strsplit(cohends$rxn,"_(?!.*_)",perl = TRUE) %>% sapply("[[",2)
cohends$dir = ifelse(cohends$dir == "pos", "(+)","(-)")
cohends = merge(cohends, rxn_metadata, by.x = "rxnrename", by.y = "reaction_no_direction")
cohends[, new_name := paste0(stringr::str_trunc(rxnrename,43)," ",dir)]
# pull in padj values with presto
presto_res = presto::wilcoxauc(X = t(rc[rn %in% c("E3","E4"),-1]), y = rc[rn %in% c("E3","E4"),rn]) %>% as.data.table
presto_res = presto_res[group == "E4"]
cohends$padj = presto_res$padj
# order and color
highlight.these.low = c("ROS detoxification","Fatty acid oxidation","Pyruvate metabolism","CoA synthesis","Oxidative phosphorylation")
highlight.these.high = c("Phenylalanine metabolism","Histidine metabolism","Lysine metabolism","Glutathione metabolism")
cohends$subsystem = paste("<span style = 'color: ",
                          ifelse(cohends$subsystem %in% highlight.these.low,"blue",
                                 ifelse(cohends$subsystem %in% highlight.these.high, "red","black")),
                          ";'>",
                          ifelse(cohends$subsystem %in% c(highlight.these.low,highlight.these.high), "<b>",""),
                          cohends$subsystem,
                          ifelse(cohends$subsystem %in% c(highlight.these.low,highlight.these.high), "</b>",""),
                          "</span>",sep = "")
cohends$subsystem = factor(cohends$subsystem, levels = cohends[, median(cohend), by = subsystem][order(V1)]$subsystem)
cohends[, color := ifelse(cohend >= 0, "E4", "E3")]
cohends[, color := ifelse(cohend >= 0, "E4", "E3")]
cohends$color = factor(cohends$color, levels = c("E4","E3"))
cohends[, shading := padj < 0.1]
cohends$shading = factor(cohends$shading, levels = c(TRUE,FALSE))

these.subsystems = c(head(levels(cohends$subsystem),n=30),tail(levels(cohends$subsystem),n=15))
cohends[subsystem %in% these.subsystems[1:30], facet := "Down in *APOE4*"]
cohends[subsystem %in% these.subsystems[31:45], facet := "Up in *APOE4*"]
cohends$facet = factor(cohends$facet, levels = c("Up in *APOE4*","Down in *APOE4*"))

(ggplot(cohends[subsystem %in% these.subsystems], aes(x=cohend,y=subsystem)) + 
    geom_point(aes(fill=color, alpha = shading),size=5,shape=21,color="black",lwd=0.25) %>% blend("source") + 
    scale_fill_manual(values = c("E4" = "red", "E3" = "blue")) + 
    scale_alpha_discrete(range=c(1,0.4)) +
    theme_bw() + 
    labs(x = "Cohen's d", y= "", fill = "Enriched in:", alpha = "Adjusted P-Value\n< 0.1?") +
    facet_wrap(~facet, scales = "free") +
    theme(strip.text = ggtext::element_markdown(size=25), axis.text.y = ggtext::element_markdown(size=22.5),
          axis.text.x = element_text(size=17.5),
          axis.title = element_text(size=25),
          legend.title = element_text(size=25),
          legend.text = element_text(size=22.5),
          legend.spacing.y = unit(0.2, 'in')) + 
    guides(fill=guide_legend(override.aes = list(size = 6),byrow=TRUE),
           alpha=guide_legend(override.aes = list(size = 6, fill = "black"),byrow=TRUE))
) %>% 
  ggsave("cohen_plot_e4e3.png",.,height=11.772,width=35.315,dpi=600)