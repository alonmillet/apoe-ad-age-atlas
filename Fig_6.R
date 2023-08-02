setwd("abeta")

library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(data.table)
library(patchwork)
library(ggplot2)

# import and merge data ----
path = "cellranger_outputs/abeta/outs/per_sample_outs/"
dir_list = list.files("cellranger_outputs/abeta/outs/per_sample_outs")

for(i in seq_along(dir_list)){
  tmp = Read10X(paste0(path,dir_list[i],"/count/sample_feature_bc_matrix"))
  tmpseu = CreateSeuratObject(counts = tmp$`Gene Expression`)
  tmpseu[["HTO"]] = CreateAssayObject(counts = tmp$`Multiplexing Capture`)
  tmpseu$orig.ident = dir_list[i]
  tmpseu = RenameCells(object = tmpseu, new.names = paste0(dir_list[i], "_", rownames(tmpseu[[]])))
  # MiQC filtering
  tmpseu[["percent.mt"]] = PercentageFeatureSet(tmpseu, pattern = "^mt-")
  tmpseu = RunMiQC(tmpseu, percent.mt = "percent.mt", nFeature_RNA = "nFeature_RNA", 
                   posterior.cutoff = 0.75, model.slot = "flexmix_model", model.type = "spline")
  tmpseu = subset(tmpseu, miQC.keep == "keep")
  # Merge as output
  if (i == 1){
    abeta = tmpseu
  } else {
    abeta = merge(x = abeta, y = tmpseu)
  }
}

# run SCT ----
abeta.list = SplitObject(abeta, split.by = "orig.ident")
for (i in names(abeta.list)) {
  abeta.list[[i]] = SCTransform(abeta.list[[i]], verbose = TRUE, vst.flavor = "v2")
}
all_genes = lapply(abeta.list, rownames) %>% Reduce(intersect, .)
abeta.features = SelectIntegrationFeatures(object.list = abeta.list, nfeatures = 3000)
abeta.list = PrepSCTIntegration(object.list = abeta.list, anchor.features = abeta.features)
abeta.anchors = FindIntegrationAnchors(object.list = abeta.list, normalization.method = "SCT", 
                                       anchor.features = abeta.features)

abeta.integrated = IntegrateData(anchorset = abeta.anchors, normalization.method = "SCT", features.to.integrate = abeta.features)

# clustering ----
abeta.integrated = RunPCA(abeta.integrated, verbose = FALSE)
abeta.integrated = RunUMAP(abeta.integrated, reduction = "pca", dims = 1:45)
abeta.integrated = FindNeighbors(abeta.integrated, dims = 1:45, verbose = TRUE)
abeta.integrated = FindClusters(abeta.integrated, verbose = TRUE, resolution = 1.4)
abeta.integrated$genotype = stringr::str_sub(abeta.integrated$orig.ident, 1, 2)
abeta.integrated$uptake = strsplit(abeta.integrated$orig.ident,"_",fixed=T) %>% sapply(.,"[[",2)
abeta.integrated$genotype = factor(abeta.integrated$genotype, levels = c("E2","E3","E4"))
abeta.integrated$uptake = factor(abeta.integrated$uptake, levels = c("hi","lo"))
abeta.integrated$orig.ident = factor(abeta.integrated$orig.ident, levels = dir_list)
saveRDS(abeta.integrated, "abeta.rds")

# project onto atlas ----
atlas = readRDS("~/adapoe.rds")
atlas = SCTransform(atlas, verbose=T, vst.flavor="v2") # unify atlas SCT models
atlas = RunUMAP(atlas, reduction = "pca", dims = 1:45, return.model = T) # rerun UMAP, saving model

anchors = FindTransferAnchors(reference = atlas, query = abeta.integrated, dims = 1:30, reference.reduction = "pca", normalization.method = "SCT")
predictions = TransferData(anchorset = anchors, refdata = atlas$subclust_id, dims = 1:30)

abeta.transfer = MapQuery(anchorset = anchors, reference = atlas, query = abeta.integrated, refdata = "subclust_id", 
                          reference.reduction = "pca", reduction.model = "umap")
abeta.transfer$predicted.id = factor(abeta.transfer$predicted.id, levels = levels(atlas$subclust_id))
saveRDS(abeta.transfer, "abeta_transfer.rds")

# annotate ----
markers = FindAllMarkers(abeta.transfer, assay = "RNA", min.pct = 0.5, logfc.threshold = 0.5, only.pos = TRUE)
idents = c(rep("Microglia",9),"Cd8 T Cells","Microglia","Neutrophils","Microglia","Immature T Cells","Macrophages","Microglia",
           "Monocytes","Mature B Cells","ILC2s","F13a1+ Monocytes","Dendritic Cells","Proliferating Immune Cells","Microglia",
           "GD T Cells","Activated NK Cells","Activated T Cells","Microglia","Macrophages","NK Cells","Cd8 T Cells","Interferon Induced T Cells",
           "Immature B Cells","Proliferating Immune Cells","Basophils")
names(idents) = levels(abeta.transfer)
abeta.transfer = RenameIdents(abeta.transfer, idents)
abeta.transfer$clust_id = Idents(abeta.transfer) %>% as.character
abeta.transfer$microglia = abeta.transfer$clust_id == "Microglia"
DimPlot(abeta.transfer, label = T)

# subset to microglia ----
abeta.mglia = subset(abeta.transfer, microglia == T)
abeta.mglia = SCTransform(abeta.mglia, method = "glmGamPoi", vst.flavor = "v2", vars.to.regress = c("percent.mt"))
abeta.mglia = RunPCA(abeta.mglia, npcs = 50, verbose = FALSE)
abeta.mglia = RunUMAP(abeta.mglia, dims = 1:46, verbose = TRUE)
abeta.mglia = FindNeighbors(abeta.mglia, dims = 1:46, verbose = TRUE)
abeta.mglia = FindClusters(abeta.mglia, verbose = TRUE, resolution = 0.8)
markers = FindAllMarkers(abeta.mglia, assay = "RNA", min.pct = 0.3, logfc.threshold = 0.3, only.pos = TRUE) %>% as.data.table
mglia.ident = c("Poised Homeostatic Microglia","Homeostatic Microglia","Ribosome-hi Homeostatic Microglia","Lars2-mid Homeostatic Microglia","DAM-1","Effector-hi TIMs",
                "Bri3-Negative Homeostatic Microglia","Interferon Induced Microglia","Effector-lo TIMs","DAM-2","Lars2-hi Homeostatic Microglia",
                "Cd74+ Microglia")
names(mglia.ident) = levels(abeta.mglia)
abeta.mglia = RenameIdents(abeta.mglia, mglia.ident)
my_mglia_levels = c("Homeostatic Microglia","Lars2-mid Homeostatic Microglia","Lars2-hi Homeostatic Microglia","Bri3-Negative Homeostatic Microglia",
                    "Ribosome-hi Homeostatic Microglia","Poised Homeostatic Microglia","DAM-1","DAM-2","Interferon Induced Microglia",
                    "Effector-lo TIMs","Effector-hi TIMs","Cd74+ Microglia")
levels(abeta.mglia) = my_mglia_levels
abeta.mglia$mglia_ident = Idents(abeta.mglia)
DimPlot(abeta.mglia, label = TRUE) + NoLegend()
# paste info back into larger seurat structure
abeta.transfer$subclust_id = as.character(abeta.transfer$clust_id)
abeta.transfer$subclust_id[Cells(abeta.mglia)] = as.character(Idents(abeta.mglia)) 
new.levels = c(levels(abeta.mglia$mglia_ident),"Proliferating Immune Cells","Neutrophils","Immature B Cells","Mature B Cells","Monocytes",
               "F13a1+ Monocytes","Dendritic Cells","Macrophages","Cd8 T Cells","GD T Cells","Interferon Induced T Cells","Activated T Cells",
               "Immature T Cells","NK Cells","Activated NK Cells","ILC2s","Basophils")
abeta.transfer$subclust_id = factor(abeta.transfer$subclust_id, levels = new.levels)
VlnPlot(abeta.transfer, "P2ry12", group.by = 'subclust_id') + NoLegend()
Idents(abeta.transfer) = abeta.transfer$subclust_id;DimPlot(abeta.transfer, label = TRUE, repel = TRUE) + NoLegend()
saveRDS(abeta.transfer,"abeta_transfer.rds")
saveRDS(abeta.mglia,"abeta_mglia.rds")

# visualize ----
## fracplot ----
abeta.mglia = readRDS("abeta_mglia.rds")
md = abeta.mglia@meta.data %>% as.data.table
setkey(md, mglia_ident, orig.ident)
dat = md[CJ(mglia_ident, orig.ident, unique = TRUE),.N,by=.EACHI]
dat[, frac := N/sum(N), by = orig.ident]
dat$geno = strsplit(dat$orig.ident %>% as.character,"_",fixed=T) %>% sapply(.,"[[",1)
dat$uptake = strsplit(dat$orig.ident %>% as.character,"_",fixed=T) %>% sapply(.,"[[",2)
dat$geno = factor(dat$geno, levels = c("E2","E3","E4"))
dat$uptake = factor(dat$uptake, levels = c("hi","lo"))
dat.cast = dcast(dat, geno + mglia_ident ~ uptake, value.var = "N")
dat.cast$rat = dat.cast$hi / (dat.cast$hi + dat.cast$lo) # what frac of the cells are in hi vs lo?
dat.cast[, tot := sum(hi) + sum(lo), by = .(geno)] # how many total cells do we have per geno?
dat.cast$frac_of_all = (dat.cast$hi + dat.cast$lo) / dat.cast$tot # what fraction of cells from this geno are in this cluster?
# calculate chisq pvals
chisq.calc = function(clust){
  tab = dat.cast[mglia_ident == clust, .(hi,lo)]
  tab$expect = rowMeans(tab)
  tab$iter = 1:nrow(tab)
  chisq = tab[,((hi-expect)^2)/expect,by=iter] %>% sum
  return(pchisq(2*chisq, df = 2,lower.tail=F))
}
these.clusters = c("Homeostatic Microglia","Poised Homeostatic Microglia",
                   "Ribosome-hi Homeostatic Microglia","Interferon Induced Microglia",
                   "Effector-lo TIMs","Effector-hi TIMs")
chisqs = data.table(mglia_ident = these.clusters,
                    chisqs = lapply(these.clusters, chisq.calc) %>% unlist %>% formatC(.,format='e',digits=2))
formatted = c()
for(i in 1:nrow(chisqs)){
  formatted[i] = bquote(paste(chi[.05]^2  == .(chisqs$chisqs[i]))) %>% deparse
}
chisqs$formatted = formatted
chisqs$mglia_ident = factor(chisqs$mglia_ident, levels = these.clusters)
coeff = 0.5 # for second axis scaling
fracplot = ggplot(dat.cast[mglia_ident %in% these.clusters], aes(x = geno)) +
  geom_bar(aes(y = rat, fill = geno), stat = "identity") + 
  geom_point(aes(y = frac_of_all/coeff), size = 6) + 
  geom_hline(yintercept = 0.5, linetype = "dashed", lwd = 1) +
  scale_y_continuous(name = "Fraction of Cells in High-Uptake Pool (             )",
                     sec.axis = sec_axis(~.*coeff, name = "Fraction of Cells in Cluster from this Genotype (\U25CF)")) +
  # leave empty spaces in first axis label to add unicode-12 symbols that don't render in R
  facet_grid(cols = vars(mglia_ident), labeller = label_wrap_gen(width = 10, multi_line = TRUE)) + 
  geom_text(data = chisqs, aes(x = 2, y = 0.8, label = formatted), parse = T, size = 4.5) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") + 
  labs(x = "")
ggsave("uptake_fracs.png",fracplot,dpi=600,width=12,height=7)

## Fig. 6D ----
fracplot_sep = ggplot(dat.cast[mglia_ident %in% these.clusters], aes(x = geno)) +
  geom_bar(aes(y = rat, fill = geno), stat = "identity", width = 0.6) + 
  # 
  geom_hline(yintercept = 0.5, linetype = "dashed", lwd = 1) +
  #scale_y_continuous(name = "Fraction of Cells in High-Uptake Pool (             )",
  #                   sec.axis = sec_axis(~.*coeff, name = "Fraction of Cells in Cluster from this Genotype (\U25CF)")) +
  # leave empty spaces in first axis label to add unicode-12 symbols that don't render in R
  facet_grid(cols = vars(mglia_ident), labeller = label_wrap_gen(width = 10, multi_line = TRUE)) + 
  geom_text(data = chisqs, aes(x = 2, y = 0.79, label = formatted), parse = T, size = 4.5) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") + 
  labs(x = "", y="Fraction of Cells in\nHigh-Uptake Pool")
dat.cast = dat.cast[, frac_of_all_normed:= 100*frac_of_all / mean(frac_of_all), by = mglia_ident]
dotplot_sep = ggplot(dat.cast[mglia_ident %in% these.clusters]) +
  geom_segment(aes(xend = geno, x = geno, yend = frac_of_all_normed, y = 0), lwd = 1.5) +
  geom_point(aes(x = geno, y = frac_of_all_normed, color=geno), size = 8) +
  facet_grid(cols = vars(mglia_ident)) +
  theme_bw(base_size=14) +
  theme(legend.position = "none", strip.background = element_blank(), strip.text.x = element_blank()) +
  labs(x = "", y = "Frequency of Cluster\n(% of Cross-Genotype Mean)")
ggsave("uptake_fracs.png",(fracplot_sep/dotplot_sep) + plot_layout(widths = c(1.25,1)),dpi=600,width=12,height=9)

## Fig. 6C ----
dat = abeta.mglia@reductions$umap@cell.embeddings %>% as.data.frame
dat$uptake = abeta.mglia$uptake
dat$geno = abeta.mglia$genotype
dat$samp = abeta.mglia$orig.ident
dat$clust = abeta.mglia$mglia_ident
abeta.mglia$label_densityplot = ifelse(abeta.mglia$mglia_ident %in% c("Lars2-hi Homeostatic Microglia","Lars2-mid Homeostatic Microglia"), 
                                       "Homeostatic Microglia",as.character(abeta.mglia$mglia_ident)) %>% 
  stringr::str_wrap(., width=18)

p1 = DimPlot(abeta.mglia, label = T, repel = T, label.box = T, label.size = 5.5, group.by = "label_densityplot") + 
  ggtitle("Joint UMAP") +
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  NoLegend()
p2 = (ggplot(dat, aes(x=UMAP_1,y=UMAP_2)) + 
        geom_point(color = "gray", alpha = 0.4) +
        stat_density_2d(geom = "density_2d_filled", contour = T, aes(fill = after_stat(level)), data = subset(dat, uptake == "hi"),
                        linewidth = 0.5, color = "ivory", alpha = 0.8) + 
        labs(title = "High Uptake")) +
  (ggplot(dat, aes(x=UMAP_1,y=UMAP_2)) + 
     geom_point(color = "gray", alpha = 0.4) +
     stat_density_2d(geom = "density_2d_filled", contour = T, aes(fill = after_stat(level)), data = subset(dat, uptake == "lo"),
                     linewidth = 0.5, color = "ivory", alpha = 0.8) + 
     labs(title = "Low Uptake")) & 
  scale_x_continuous(limits = c(-6.3,7.3)) &
  scale_y_continuous(limits = c(-5.3,5.3)) &
  theme_classic(base_size = 20) &
  theme(plot.title = element_text(hjust = 0.5)) &
  NoLegend() & 
  scale_fill_distiller(palette="Spectral", direction=-1)
p3 = (ggplot(dat, aes(x=UMAP_1,y=UMAP_2)) + 
        geom_point(color = "gray", alpha = 0.4) +
        stat_density_2d(geom = "density_2d_filled", contour = T, aes(fill = after_stat(level)), data = subset(dat, geno == "E2"),
                        linewidth = 0.5, color = "ivory", alpha = 0.8) + 
        labs(title = "ApoE2")) +
  (ggplot(dat, aes(x=UMAP_1,y=UMAP_2)) + 
     geom_point(color = "gray", alpha = 0.4) +
     stat_density_2d(geom = "density_2d_filled", contour = T, aes(fill = after_stat(level)), data = subset(dat, geno == "E3"),
                     linewidth = 0.5, color = "ivory", alpha = 0.8) + 
     labs(title = "ApoE3")) +
  (ggplot(dat, aes(x=UMAP_1,y=UMAP_2)) + 
     geom_point(color = "gray", alpha = 0.4) +
     stat_density_2d(geom = "density_2d_filled", contour = T, aes(fill = after_stat(level)), data = subset(dat, geno == "E4"),
                     linewidth = 0.5, color = "ivory", alpha = 0.8) + 
     labs(title = "ApoE4")) & 
  scale_x_continuous(limits = c(-6.3,7.3)) &
  scale_y_continuous(limits = c(-5.3,5.3)) &
  theme_classic(base_size = 20) &
  theme(plot.title = element_text(hjust = 0.5)) &
  NoLegend() & 
  scale_fill_distiller(palette="Spectral", direction=-1)

(p1 / p3) %>% ggsave("umap_density_bygeno.png",.,dpi=600,width=10,height=10)
((p1 | p2) + plot_layout(widths = c(1.5,2))) %>% ggsave("umap_density_byuptake.png",.,dpi=600,width=15.4,height=6.1)
(((p1 | p2) + plot_layout(widths = c(1,2))) / p3) %>% ggsave("umap_density_combo.png",.,dpi=600,width=14,height=14)

## Fig. 6B ----
library(stringr)
abeta.transfer$label_densityplot = str_wrap(abeta.transfer$clust_id, width=10)
(DimPlot(abeta.transfer, label = T, repel = T, label.size = 5, group.by = "label_densityplot") + NoLegend() + ggtitle("")) %>% 
  ggsave("umap_allcells.png",.,dpi=600,width=6,height=11.25)
