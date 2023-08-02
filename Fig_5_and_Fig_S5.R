library(ggplot2)
library(dplyr)
library(data.table)
library(ggbreak)

setwd("xenium")

e3 = fread('cooccurrence_raw_e3.csv')
e4 = fread('cooccurrence_raw_e4.csv')

x = list(e3,e4)
mean = (Reduce(`+`, x) / length(x)) %>% as.data.table
raw = melt.data.table(mean,id.vars="interval",measure.vars=colnames(mean)[3:24])

maxes = raw[interval<1000,max(abs(1-value)),by=variable]

# Fig. 5C ----
use.these = maxes[V1>0.15,variable]
ggplot(raw[variable %in% use.these], aes(x=interval,y=value,color=variable)) + 
  geom_point() +
  geom_smooth(method="loess") +
  theme_classic(base_size=24) +
  labs(x="Distance",y="Increased Odds over Expectation",color="Cluster") +
  scale_x_continuous(expand = c(0,0)) +
  guides(color = guide_legend(ncol=2)) +
  theme(legend.position = "top", plot.margin = margin(r=1,unit="cm"))
ggsave("densities.png",dpi=600,width=8,height=8)

# Fig. S5E ----
use.these = maxes[V1>0.05,variable]
ggplot(raw[interval < 1000 & variable %in% use.these],aes(x=reorder(variable,-value),y=value)) +
  geom_boxplot(aes(fill=variable)) +
  theme_bw(base_size = 16) +
  labs(x="",y="Increased Odds\nover Expectation\nin First 1000 Pixels") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("densities_boxplot.png",dpi=600,width=8,height=6)

# Fig. 5D ----
library(scales)
e3 = melt.data.table(e3,id.vars="interval",measure.vars=colnames(e3)[3:24])
e4 = melt.data.table(e4,id.vars="interval",measure.vars=colnames(e4)[3:24])

pcalcer = function(clust){
  pval = t.test(e3[interval<3000 & variable == clust,value],e4[interval<3000 & variable == clust,value])$p.value
  if(pval>0.05){
    return("ns")
  } else if(pval>0.01) {
    return("*")
  } else if(pval>0.001) {
    return("**")
  } else {
    return("***")
  }
}

e3_means = e3[interval<1000,mean(value,na.rm=T),by=variable]
e4_means = e4[interval<1000,mean(value,na.rm=T),by=variable]
merged_means = merge.data.table(e3_means,e4_means,by="variable")
colnames(merged_means) = c("cluster","ApoE3","ApoE4")
merged_means$diffs = merged_means$ApoE4 - merged_means$ApoE3
merged_means = merged_means[abs(diffs)>0.05]
merged_means[, means := mean(c(ApoE3,ApoE4)), by = cluster]
merged_means$cluster = forcats::fct_reorder(merged_means$cluster,merged_means$means)
merged_means[, pval := pcalcer(cluster), by = cluster]
merged_means[, label_pos := mean(c(ApoE3,ApoE4)), by = cluster]
melted_merged_means = melt.data.table(merged_means,id.vars="cluster",measure.vars=c('ApoE3','ApoE4'))
merged_means = merged_means[order(means)]
colors = hue_pal()(2)
cols_for_plot = ifelse(merged_means$diffs > 0, colors[2], colors[1])

ggplot(data = merged_means) +
  geom_vline(xintercept = 1, lwd = 1, linetype = "dashed") +
  geom_segment(aes(y=cluster,yend=cluster,x=ApoE3,xend=ApoE4), lwd = 1.5) +
  geom_point(data = melted_merged_means, aes(y=cluster,x=value,fill=variable), shape=21, size=7.5, stroke = 1.5) +
  geom_text(aes(x=label_pos,y=cluster,label=pval), vjust = -0.5, size = 6) +
  theme_bw(base_size = 24) +
  labs(x="Increased Odds over Expectation\nin First 1000 Pixels",y="",fill="Genotype") +
  scale_x_break(c(1.2, 1.39)) +
  expand_limits(x=1.5) +
  theme(axis.text.y = element_text(color = cols_for_plot))
ggsave("e4vse3_barbells.png",dpi=600,width=15,height=11.25)

# now switching to metadata
md = fread("xenium_md.csv")
library(stringr)
# Fig. 5B ----
freqs = md[,.N,by=.(subclust,sample,genotype,superclust)]
freqs[, freq := 100*N/sum(N), by = .(sample,genotype,superclust)]
good_samps = freqs[subclust == "TIMs" & N > 250] # only keep samples with a good representation of TIMs
good_samps %>% 
  mutate(genotype = str_replace(genotype, "E", "ApoE")) %>%
  ggplot(aes(x=genotype,y=freq,fill=genotype)) + 
  geom_bar(stat = "summary", fun = "mean") +
  #stat_summary(geom = "errorbar", fun.data = mean_se, width = 0.4) +
  theme_bw(base_size = 24) +
  labs(x="",y="Frequency of TIMs (% of Microglia)") +
  theme(legend.position = "none")
ggsave("tims_per_genotype.png",dpi=600,width=8,height=8)

# make spatial plots
library(patchwork)
md = fread("xenium_md.csv")
order = c("L1 Neurons","L2 Neurons","L3 Neurons","L4 Neurons","L5 Neurons","L6 Neurons","Mixed Border Neurons","Mixed Neurons",
          "CNTN2-hi Oligodendrocytes","EFHD1-hi Oligodendrocytes","ERMN-hi Oligodendrocytes","AQP4-hi Astrocytes","GJA1-hi Astrocytes",
          "APP-hi Astrocytes","APOE-hi Astrocytes",'Homeostatic Microglia','APP-hi Homeostatic Microglia',
          "FLT1-hi Inflammatory Microglia","APOE-hi Inflammatory Microglia","TIMs","OPCs","VLMCs")
md$subclust = factor(md$subclust, levels = order)
md$size = ifelse(md$subclust == "TIMs", "A", "B")

# Fig. 5A ----
# zoomed views
p1 = ggplot(md[x_centroid > 3000 & x_centroid < 8000 & y_centroid > 5000 & y_centroid < 10000 & sample == "E4_2"]) + 
  geom_point(aes(x=x_centroid,y=y_centroid, color=subclust, size=size)) +
  scale_size_manual(values = c(8,0.5)) +
  theme_classic(base_size=24) + 
  guides(size="none", color=guide_legend(ncol=1, override.aes = list(size=5))) +
  labs(x="Horizontal Pixels",y="Vertical Pixels",color="Cluster",title="") +
  scale_color_manual(values = scales::hue_pal()(22), limits = levels(md$subclust), labels = levels(md$subclust)) +
  theme(legend.position = "none") +
  scale_y_reverse()
ggsave("e4_spatial_zoom.png",p1,dpi=600,height=8,width=8)
p2 = ggplot(md[x_centroid > 1500 & x_centroid < 6500 & y_centroid > 0 & y_centroid < 5000 & sample == "E3_2"]) + 
  geom_point(aes(x=x_centroid,y=y_centroid, color=subclust, size=size)) +
  scale_size_manual(values = c(8,0.5)) +
  theme_classic(base_size=24) + 
  guides(size="none", color=guide_legend(ncol=1, override.aes = list(size=5))) +
  labs(x="Horizontal Pixels",y="Vertical Pixels",color="Cluster",title="") +
  scale_color_manual(values = scales::hue_pal()(22), limits = levels(md$subclust), labels = levels(md$subclust)) +
  theme(legend.position = "none") +
  scale_y_reverse()
ggsave("e3_spatial_zoom.png",p2,dpi=600,height=8,width=8)
# zoomed-out views with highlight
inset_e3 = ggplot(md[sample == "E3_2"]) +
  geom_point(aes(x=x_centroid,y=y_centroid, color=subclust), size = 0.05) +
  geom_rect(aes(xmin=1500,xmax=6500,ymin=0,ymax=5000),color="black", fill = NA, linetype = "dashed", lwd = 2) +
  theme_classic(base_size=24) + 
  guides(color=guide_legend(ncol=1, override.aes = list(size=5))) +
  labs(x="Horizontal Pixels",y="Vertical Pixels",color="Cluster",title="ApoE3 (Donor 211)") +
  theme(legend.position = "none", plot.margin = margin(r=1.5,unit='cm')) +
  scale_y_reverse()
ggsave("e3_inset_highlight.png",inset_e3,dpi=600,width=8,height=8)
inset_e4 = ggplot(md[sample == "E4_2"]) +
  geom_point(aes(x=x_centroid,y=y_centroid, color=subclust), size = 0.05) +
  geom_rect(aes(xmin=3000,xmax=8000,ymin=5000,ymax=10000),color="black", fill = NA, linetype = "dashed", lwd = 2) +
  theme_classic(base_size=24) + 
  guides(color=guide_legend(ncol=1, override.aes = list(size=5))) +
  labs(x="Horizontal Pixels",y="Vertical Pixels",color="Cluster",title="ApoE4 (Donor 173)") +
  theme(legend.position = "none", plot.margin = margin(r=1.5,unit='cm')) +
  scale_y_reverse()
ggsave("e4_inset_highlight.png",inset_e4,dpi=600,width=8,height=8)
# just make legend
library(cowplot)
library(grid)
library(gridExtra) 
legend_source = ggplot(md) +
  geom_point(aes(x=x_centroid,y=y_centroid, color=subclust)) +
  labs(color = "Cluster") +
  guides(color=guide_legend(ncol=1, override.aes = list(size=8))) +
  theme_classic(base_size=24)
legend = cowplot::get_legend(legend_source)  
png("spatial_legend.png",res=600,width=5,height=8,units='in')
grid.newpage()
grid.draw(legend)
dev.off()

# Fig. S5B ----
spatplotter = function(samp,donor){
  p = ggplot(md[sample == samp]) + 
    geom_point(aes(x=x_centroid,y=y_centroid, color=subclust), size = 0.05) +
    theme_classic(base_size=24) + 
    guides(color=guide_legend(ncol=1, override.aes = list(size=5))) +
    labs(x="Horizontal Pixels",y="Vertical Pixels",color="Cluster",title=donor) +
    scale_y_reverse()
  return(p)
}
samps = c("E3_1","E3_2","E3_3","E4_1","E4_2","E4_3")
donors = c("ApoE3 (Donor 148)","ApoE3 (Donor 211)","ApoE3 (Donor 251)","ApoE4 (Donor 83)","ApoE4 (Donor 173)","ApoE4 (Donor 191)")

wrap_plots(mapply(spatplotter, samp = samps, donor = donors, SIMPLIFY = F)) + plot_layout(guides="collect") &
  scale_color_manual(values = scales::hue_pal()(22), limits = levels(md$subclust), labels = levels(md$subclust)) &
  theme(plot.margin = margin(0.25,1,0.25,0.25, unit = "cm"), legend.margin = margin(l = 3, unit='cm'))
ggsave("all_sections.png",dpi=600,width=25,height=10)

# Fig. S5C-D ----
# made entirely in squidpy, see accompanying python code 

# Fig. S5A ----
library(Seurat)
library(data.table)
library(dplyr)

setwd("DLPFC_2") # see code for Figure 4

dlpfc2_md = fread("dlpfc_md.csv")
og = readRDS("dlpfc_mglia.rds")

panel = fread("Xenium_hBrain_v1_metadata.csv")
subset.matrix = og@assays$RNA@counts[panel$Genes,]
subset.matrix = subset.matrix[,-which(colSums(subset.matrix) == 0)]
ogsub = CreateSeuratObject(subset.matrix)
rownames(dlpfc2_md) = paste0(dlpfc2_md$libraryBatch,"_",dlpfc2_md$cellBarcode)
ogsub = AddMetaData(ogsub, dlpfc2_md)
Idents(ogsub) = ogsub$predicted.id
ogsub = SCTransform(ogsub, vst.flavor = "v2", method = "glmGamPoi", vars.to.regress = c("libraryBatch","individualID"))
saveRDS(ogsub,"dlpfc2_for_xenium.rds")

ogsub = subset(ogsub, predicted.id %in% c("Homeostatic Microglia","DAM-2","TIMs"))
ogsub$predicted.id = factor(ogsub$predicted.id, levels = c("Homeostatic Microglia","DAM-2","TIMs"))
(VlnPlot(ogsub, c("RNASET2","PTPRC","GPR183"), ncol = 1) + NoLegend()) & labs(x="")
ggsave("markers_from_dlpfc2.png",dpi=600,width=6,height=14)
# adjoined to violin plots made in squidpy, see accompanying python code
