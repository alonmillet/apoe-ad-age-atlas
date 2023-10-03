setwd("aducanumab")

library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(data.table)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(stringr)

# import and merge data ----
path = "cellranger_outputs/aducanumab/outs/per_sample_outs/"
dir_list = list.files("cellranger_outputs/aducanumab/outs/per_sample_outs")

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
    aduc = tmpseu
  } else {
    aduc = merge(x = aduc, y = tmpseu)
  }
}

# run SCT ----
aduc.list = SplitObject(aduc, split.by = "orig.ident")
for (i in names(aduc.list)) {
  aduc.list[[i]] = SCTransform(aduc.list[[i]], verbose = TRUE, vst.flavor = "v2")
}
all_genes = lapply(aduc.list, rownames) %>% Reduce(intersect, .)
aduc.features = SelectIntegrationFeatures(object.list = aduc.list, nfeatures = 3000)
aduc.list = PrepSCTIntegration(object.list = aduc.list, anchor.features = aduc.features)
aduc.anchors = FindIntegrationAnchors(object.list = aduc.list, normalization.method = "SCT", 
                                      anchor.features = aduc.features)

aduc.integrated = IntegrateData(anchorset = aduc.anchors, normalization.method = "SCT", features.to.integrate = aduc.features)

# clustering ----
aduc.integrated = RunPCA(aduc.integrated, verbose = FALSE)
ElbowPlot(aduc.integrated, ndims = 50)
aduc.integrated = RunUMAP(aduc.integrated, reduction = "pca", dims = 1:45)
aduc.integrated = FindNeighbors(aduc.integrated, dims = 1:45, verbose = TRUE)
aduc.integrated = FindClusters(aduc.integrated, verbose = TRUE, resolution = 1.4)
aduc.integrated$genotype = stringr::str_sub(aduc.integrated$orig.ident, 1, 2)
aduc.integrated$tx = strsplit(aduc.integrated$orig.ident,"_",fixed=T) %>% sapply(.,"[[",2)
aduc.integrated$genotype = factor(aduc.integrated$genotype, levels = c("E2","E3","E4"))
aduc.integrated$tx = factor(aduc.integrated$tx, levels = c("aducanumab","isocon"))
aduc.integrated$orig.ident = factor(aduc.integrated$orig.ident, levels = dir_list)
saveRDS(aduc.integrated, "aduc.rds")

# project onto atlas ----
atlas = readRDS("~/adapoe.rds")
atlas = SCTransform(atlas, verbose=T, vst.flavor="v2") # unify atlas SCT models
atlas = RunUMAP(atlas, reduction = "pca", dims = 1:45, return.model = T) # rerun UMAP, saving model

anchors = FindTransferAnchors(reference = atlas, query = aduc.integrated, dims = 1:30, reference.reduction = "pca", normalization.method = "SCT")
predictions = TransferData(anchorset = anchors, refdata = atlas$subclust_id, dims = 1:30)

aduc.transfer = MapQuery(anchorset = anchors, reference = atlas, query = aduc.integrated, refdata = "subclust_id", 
                         reference.reduction = "pca", reduction.model = "umap")
aduc.transfer$predicted.id = factor(aduc.transfer$predicted.id, levels = levels(atlas$subclust_id))
saveRDS(aduc.transfer, "aduc_transfer.rds")

# annotate ----
markers = FindAllMarkers(aduc.transfer, assay = "RNA", min.pct = 0.5, logfc.threshold = 0.5, only.pos = TRUE)
idents = c(rep("Microglia",12),"NK Cells","Late Neutrophils","Microglia","Cd8 T Cells",rep("Microglia",4),"Dendritic Cells",
           "Early Neutrophils","Monocytes","Mature B Cells","F13a1+ Monocytes","Macrophages","GD T Cells","DP T Cells","Activated NK Cells",
           "Proliferating Immune Cells","Immature B Cells","Erythroid Cells","Immature T Cells")
names(idents) = levels(aduc.transfer)
aduc.transfer = RenameIdents(aduc.transfer, idents)
aduc.transfer$clust_id = Idents(aduc.transfer) %>% as.character
aduc.transfer$microglia = aduc.transfer$clust_id == "Microglia"
DimPlot(aduc.transfer, label = T)

# subset to microglia ----
aduc.mglia = subset(aduc.transfer, microglia == T)
aduc.mglia = SCTransform(aduc.mglia, method = "glmGamPoi", vst.flavor = "v2", vars.to.regress = c("percent.mt"))
aduc.mglia = RunPCA(aduc.mglia, npcs = 50, verbose = FALSE)
aduc.mglia = RunUMAP(aduc.mglia, dims = 1:46, verbose = TRUE)
aduc.mglia = FindNeighbors(aduc.mglia, dims = 1:46, verbose = TRUE)
aduc.mglia = FindClusters(aduc.mglia, verbose = TRUE, resolution = 0.6)
markers = FindAllMarkers(aduc.mglia, assay = "RNA", min.pct = 0.3, logfc.threshold = 0.3, only.pos = TRUE) %>% as.data.table
mglia.ident = c("Homeostatic Microglia","DAMs","Poised-like Homeostatic Microglia",
                "Lars2-mid Homeostatic Microglia","Effector-hi TIMs","Effector-lo TIMs","Ccl3/4+ Inflammatory Microglia","Lars2-hi Homeostatic Microglia",
                "Bri3-Negative Homeostatic Microglia","Lars2-hi Homeostatic Microglia","Interferon Induced Microglia",
                "Poised-like Homeostatic Microglia","DAMs","Rgs1-hi Homeostatic Microglia","Jchain+ Microglia") 
names(mglia.ident) = levels(aduc.mglia)
aduc.mglia = RenameIdents(aduc.mglia, mglia.ident)
my_mglia_levels = c("Homeostatic Microglia","Lars2-mid Homeostatic Microglia","Lars2-hi Homeostatic Microglia","Bri3-Negative Homeostatic Microglia",
                    "Rgs1-hi Homeostatic Microglia","Poised-like Homeostatic Microglia","DAMs","Interferon Induced Microglia","Ccl3/4+ Inflammatory Microglia",
                    "Effector-lo TIMs","Effector-hi TIMs","Jchain+ Microglia")
levels(aduc.mglia) = my_mglia_levels
aduc.mglia$mglia_ident = Idents(aduc.mglia)
DimPlot(aduc.mglia, label = TRUE) + NoLegend()
## paste info back into larger seurat structure ----
aduc.transfer$subclust_id = as.character(aduc.transfer$clust_id)
aduc.transfer$subclust_id[Cells(aduc.mglia)] = as.character(Idents(aduc.mglia)) 
new.levels = c(levels(aduc.mglia$mglia_ident),"Proliferating Immune Cells","Early Neutrophils","Late Neutrophils","Immature B Cells","Mature B Cells","Monocytes","F13a1+ Monocytes",
               "Dendritic Cells","Macrophages","Cd8 T Cells","DP T Cells","GD T Cells","Immature T Cells","NK Cells","Activated NK Cells","Erythroid Cells")
aduc.transfer$subclust_id = factor(aduc.transfer$subclust_id, levels = new.levels)
VlnPlot(aduc.transfer, "P2ry12", group.by = 'subclust_id') + NoLegend()
Idents(aduc.transfer) = aduc.transfer$subclust_id;DimPlot(aduc.transfer, label = TRUE, repel = TRUE) + NoLegend()
# run ALRA imputation (used for some packages e.g. scriabin)
aduc.mglia = RunALRA(aduc.mglia, assay = "RNA")
aduc.transfer = RunALRA(aduc.transfer, assay = "RNA")
DefaultAssay(aduc.mglia) = "SCT"
DefaultAssay(aduc.transfer) = "SCT"
# save
saveRDS(aduc.transfer,"aduc_transfer.rds")
saveRDS(aduc.mglia,"aduc_mglia.rds")

# Fig. 7A ----
aduc.transfer$label_densityplot = str_wrap(aduc.transfer$clust_id, width=10)
(DimPlot(aduc.transfer, label = T, repel = T, label.size = 5, group.by = "label_densityplot") + NoLegend() + 
    theme(plot.title = element_blank(),
          axis.text = element_text(size=16),
          axis.title = element_text(size=20))) %>% 
  ggsave("umap_allcells.png",.,dpi=600,width=8,height=7.27)

# Fig. 7B ----
aduc.mglia = readRDS("aduc_mglia.rds")
dat = aduc.mglia@reductions$umap@cell.embeddings %>% as.data.frame
dat$tx = aduc.mglia$tx
dat$geno = aduc.mglia$genotype
dat$samp = aduc.mglia$orig.ident
dat$clust = aduc.mglia$mglia_ident
aduc.mglia$label_densityplot = ifelse(aduc.mglia$mglia_ident %in% c("Lars2-hi Homeostatic Microglia","Lars2-mid Homeostatic Microglia","Rgs1-hi Homeostatic Microglia"), 
                                      "Homeostatic Microglia",as.character(aduc.mglia$mglia_ident)) %>% 
  stringr::str_wrap(., width=10)

p1 = DimPlot(aduc.mglia, label = T, repel = T, label.box = T, group.by = "label_densityplot", label.size = 5) + 
  ggtitle("Joint UMAP") +
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  NoLegend()
p2 = (ggplot(dat, aes(x=UMAP_1,y=UMAP_2)) + 
        geom_point(color = "gray", alpha = 0.4) +
        stat_density_2d(geom = "density_2d_filled", contour = T, aes(fill = after_stat(level)), data = subset(dat, tx == "isocon"),
                        linewidth = 0.5, color = "ivory", alpha = 0.8) + 
        labs(title = "Isotype Control")) +
  (ggplot(dat, aes(x=UMAP_1,y=UMAP_2)) + 
     geom_point(color = "gray", alpha = 0.4) +
     stat_density_2d(geom = "density_2d_filled", contour = T, aes(fill = after_stat(level)), data = subset(dat, tx == "aducanumab"),
                     linewidth = 0.5, color = "ivory", alpha = 0.8) + 
     labs(title = "Aducanumab")) & 
  scale_x_continuous(limits = c(-7.3,5.3)) &
  scale_y_continuous(limits = c(-8.3,6.3)) &
  theme_classic(base_size = 20) &
  theme(plot.title = element_text(hjust = 0.5)) &
  NoLegend() & 
  scale_fill_distiller(palette="Spectral", direction=-1)

p2_alt = (ggplot(dat, aes(x=UMAP_1,y=UMAP_2)) + 
            geom_point(color = "gray", alpha = 0.4) +
            stat_density_2d(geom = "density_2d_filled", contour = T, aes(fill = after_stat(level)), data = subset(dat, tx == "isocon"),
                            linewidth = 0.5, color = "ivory", alpha = 0.8) + 
            labs(title = "Isotype Control")) /
  (ggplot(dat, aes(x=UMAP_1,y=UMAP_2)) + 
     geom_point(color = "gray", alpha = 0.4) +
     stat_density_2d(geom = "density_2d_filled", contour = T, aes(fill = after_stat(level)), data = subset(dat, tx == "aducanumab"),
                     linewidth = 0.5, color = "ivory", alpha = 0.8) + 
     labs(title = "Aducanumab")) & 
  scale_x_continuous(limits = c(-7.3,5.3)) &
  scale_y_continuous(limits = c(-8.3,6.3)) &
  theme_classic(base_size = 20) &
  theme(plot.title = element_text(hjust = 0.5)) &
  NoLegend() & 
  scale_fill_distiller(palette="Spectral", direction=-1)

(p1 | p2_alt) %>% ggsave("umap_density_bytx.png",.,dpi=600,width=10,height=10)

## split by genotype ----
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
  scale_x_continuous(limits = c(-7.3,5.3)) &
  scale_y_continuous(limits = c(-8.3,6.3)) &
  theme_classic(base_size = 20) &
  theme(plot.title = element_text(hjust = 0.5)) &
  NoLegend() & 
  scale_fill_distiller(palette="Spectral", direction=-1)

(p1 / p3) %>% ggsave("umap_density_bygeno.png",.,dpi=600,width=10,height=10)
(((p1 | p2) + plot_layout(widths = c(1,2))) / p3) %>% ggsave("umap_density_combo.png",.,dpi=600,width=18,height=12.9)

# Fig. S6A ----
poised = FindMarkers(aduc.mglia, ident.1 = "Homeostatic Microglia", ident.2 = "Poised-like Homeostatic Microglia", 
                     test.use = "MAST", assay = "RNA", min.pct = 0.01, logfc.threshold = 0.01) %>% as.data.table(keep.rownames = T)
poised$logp = -log10(poised$p_val_adj)
pval_cutoff = quantile(poised$logp, 0.975)
logfc_cutoff = quantile(abs(poised$avg_log2FC), 0.975)
poised[, rank := frank(logp) + frank(abs(avg_log2FC))]
poised[, in_zone := ifelse(abs(avg_log2FC) >= logfc_cutoff & logp >= pval_cutoff, TRUE, FALSE)]
poised$label = FALSE
poised[avg_log2FC > 0, label := ifelse(rank >= .SD$rank[kit::topn(.SD$rank, 12, decreasing = T)[12]], TRUE, FALSE)]
poised[avg_log2FC < 0, label := ifelse(rank >= .SD$rank[kit::topn(.SD$rank, 12, decreasing = T)[12]], TRUE, FALSE)]
poised[, group := ifelse(abs(avg_log2FC) >= logfc_cutoff & logp >= pval_cutoff, "Hit", "Filter")]
poised$group = factor(poised$group, levels = c("Hit","Filter"))
vol = ggplot(poised, aes(x = avg_log2FC, y = logp)) +
  geom_point(aes(color = group)) +
  scale_color_manual(values = c("royalblue","gray")) +
  geom_text_repel(label = ifelse(poised$label==TRUE,poised$rn,''),
                  size = 5, force = 15, box.padding = 0.5, point.padding = 0.5, segment.size = 0.1, 
                  min.segment.length = 0.1, max.iter = 1e7, max.overlaps = Inf) +
  geom_vline(xintercept = -logfc_cutoff) + 
  geom_vline(xintercept = logfc_cutoff) + 
  geom_hline(yintercept = pval_cutoff) +
  geom_text(data = data.frame(xpos = c(-(logfc_cutoff+0.05), logfc_cutoff+0.05), 
                              ypos = c(0, 0), hjustvar = c(1,0), vjustvar = c(0,0), 
                              annotateText = c(paste0("Poised-like Homeostatic ",sprintf('\U2190')),
                                               paste0(sprintf('\U2192')," Homeostatic"))),
            aes(x=xpos, y=ypos, hjust=hjustvar, vjust=vjustvar, label=annotateText), size = 5) +
  theme_bw(base_size = 14) + 
  labs(x = "Log-2 Fold Change", y = "Log-10 Adjusted P-Value") + 
  theme(legend.position = "none")
ggsave("poised_vol.png",vol,dpi=600,width=7,height=6.5)
# poised look a lot like IR2 from https://www.cell.com/immunity/pdf/S1074-7613(18)30485-0.pdf
# respond to inflammation/injury after demyelination

# make heatmap showing distribution of clusters ----
md = aduc.transfer@meta.data %>% as.data.table
setkey(md, subclust_id, tx)
dat = md[microglia == TRUE][CJ(subclust_id, tx, unique = TRUE),.N,by=.EACHI]
dat[, frac := N/sum(N), by = tx]
dat = dcast.data.table(dat, tx ~ subclust_id, value.var = "frac")
datm = as.matrix(dat[,-1])
rownames(datm) = dat$tx
datm.scaled = scale(t(datm))
colnames(datm.scaled) = c("Aduc.","Iso. Con.")
rownames(datm.scaled) = stringr::str_wrap(rownames(datm.scaled),width=20)
png("mglia_freqs_zscore.png", width = 5, height = 6.67, res = 600, units = "in")
hm = Heatmap(datm.scaled, name = "Z-Scored\nFraction\nof Cells", heatmap_width = unit(4,"in"), heatmap_height = unit(6.67,"in"),
             row_order = rownames(datm.scaled), column_order = c("Iso. Con.","Aduc."),
             rect_gp = gpar(col = "black", lwd = 0.5), column_names_rot = 0, column_names_centered = T,
             row_split = factor(c(rep("Homeostatic",6),rep("DAMs",1),rep("Inflam.",2),
                                  rep("TIMs",2),rep("Other",1)),
                                levels = c("Homeostatic","DAMs","Inflam.","TIMs","Other")),
             row_title = "", column_title = "", cluster_row_slices = FALSE)
draw(hm)#, padding = unit(c(5, 20, 2, 2), "mm"))
dev.off()

## Fig. 7C ----
md = aduc.transfer@meta.data %>% as.data.table
setkey(md, subclust_id, tx)
dat = md[microglia == TRUE][CJ(subclust_id, tx, unique = TRUE),.N,by=.EACHI]
dat[, frac := N/sum(N), by = tx]
dat = dat[subclust_id %in% c("Homeostatic Microglia","Effector-lo TIMs",
                             "Poised-kike Homeostatic Microglia","DAMs","Effector-hi TIMs")]
dat$subclust_id = factor(dat$subclust_id, levels = rev(c("Homeostatic Microglia","Effector-lo TIMs",
                                                         "Poised-like Homeostatic Microglia","DAMs","Effector-hi TIMs")))
(ggplot(dat, aes(x=subclust_id,y=frac,group=tx,fill=tx)) + 
    geom_bar(position="dodge",stat="identity",color="black", width = 0.7) + 
    scale_fill_manual(values = rev(scales::hue_pal()(2)), 
                      guide = guide_legend(reverse = TRUE), labels = c("Aducanumab","Isotype Control")) +
    coord_flip() + 
    theme_bw(base_size = 14) +
    labs(y="Frequency of Cells",x="",fill="Treatment") +
    theme(legend.position = "top", legend.direction = "vertical")) %>% 
  ggsave("mglia_freqs_barplot.png",.,dpi=600,width=5,height=6.67)

# Run CellChat ----
setwd("aducanumab/CellChat")

library(Seurat)
library(data.table)
library(dplyr)
library(CellChat)
library(patchwork)

aduc.transfer = readRDS("aducanumab/aduc_transfer.rds")

# prep individual cellchat objects ----
aduc.transfer$subclust_id_overwrite = ifelse(aduc.transfer$subclust_id == "Immature B Cells", 
                                             "Mature B Cells", as.character(aduc.transfer$subclust_id))
aduc.transfer$subclust_id_overwrite = ifelse(aduc.transfer$subclust_id_overwrite == "Immature T Cells",
                                             "DP T Cells", aduc.transfer$subclust_id_overwrite)
aduc.transfer = subset(aduc.transfer, subset = subclust_id_overwrite != "Effector-hi TIMs")
aduc.transfer$subclust_id_overwrite = factor(aduc.transfer$subclust_id_overwrite, 
                                             levels = levels(aduc.transfer$subclust_id)[-c(11,16,25)])

groups = unique(aduc.transfer$orig.ident)
for(i in seq_along(groups)){
  temp = subset(aduc.transfer, orig.ident == groups[i])
  temp = createCellChat(object = temp, group.by = "subclust_id_overwrite", assay = "RNA")
  temp@DB = CellChatDB.mouse
  temp = subsetData(temp)
  future::plan("multisession", workers = 4)
  temp = identifyOverExpressedGenes(temp)
  temp = identifyOverExpressedInteractions(temp)
  temp@idents = droplevels(temp@idents, exclude = setdiff(levels(temp@idents),unique(temp@idents)))
  temp = computeCommunProb(temp)
  temp = computeCommunProbPathway(temp)
  temp = aggregateNet(temp)
  saveRDS(temp, paste0(groups[i],"_cellchat.rds")) # run on base aduc.transfer, without merging idents
}

# use liftcellchat to bring in isocon and aducanumab into same merged object ----
c1 = readRDS("E2_aducanumab_cellchat.rds")
c2 = readRDS("E2_isocon_cellchat.rds")
c3 = readRDS("E3_aducanumab_cellchat.rds")
c4 = readRDS("E3_isocon_cellchat.rds")
c5 = readRDS("E4_aducanumab_cellchat.rds")
c6 = readRDS("E4_isocon_cellchat.rds")

c2 = liftCellChat(c2, levels(c1@idents))
c4 = liftCellChat(c4, levels(c1@idents))
c6 = liftCellChat(c6, levels(c1@idents))

cc.list = list(e2_aduc = c1, e3_aduc = c3, e4_aduc = c5, e2_isocon = c2, e3_isocon = c4, e4_isocon = c6)
cc.list = lapply(cc.list, netAnalysis_computeCentrality) #takes ~1min total
names(cc.list) = c("E2+Aducanumab","E3+Aducanumab","E4+Aducanumab","E2+Isocon","E3+Isocon","E4+Isocon")
cc = mergeCellChat(cc.list, add.names = c("E2+Aducanumab","E3+Aducanumab","E4+Aducanumab","E2+Isocon","E3+Isocon","E4+Isocon"))

## Fig. 7D ----
count_and_weight_get = function(samplename){
  count = cc@net[[samplename]]$count %>% sum
  weight = cc@net[[samplename]]$weight %>% sum
  return(c(count,weight))
}
outs = lapply(names(cc.list),count_and_weight_get) %>% do.call(cbind,.) %>% t
colnames(outs) = c("count","weights")
outs = as.data.table(outs)
outs$avgweight = outs$weights/outs$count
outs$samp = c("E2 +\nAduc.","E3 +\nAduc.","E4 +\nAduc.","E2 +\nIso. Con.","E3 +\nIso. Con.","E4 +\nIso. Con.")
outs = melt.data.table(outs, id.vars = "samp", measure.vars = c("avgweight","count"))
outs[, toplot := 100*value / max(value), by = variable]
outs[variable == "avgweight", value := round(value, digits = 4)]
outs$samp = factor(outs$samp, levels = rev(c("E2 +\nAduc.","E3 +\nAduc.","E4 +\nAduc.","E2 +\nIso. Con.","E3 +\nIso. Con.","E4 +\nIso. Con.")))
ints = ggplot(outs, aes(x=samp,y=toplot,fill=variable)) + 
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(y = toplot-10,label=value), position = position_dodge(width = 0.9), size = 7) +
  scale_fill_manual(labels = c("Mean Strength","Interaction Count"), 
                    values = rev(scales::hue_pal()(2)),
                    guide = guide_legend(reverse = TRUE)) +
  labs(x = "", y="Percent of Maximum", fill = "Feature") +
  theme_bw(base_size = 25) +
  coord_flip() +
  theme(legend.position = "top", legend.direction = "vertical")
ggsave("interactions_and_weights.png",ints,dpi=600,width=8,height=11)

## Fig. 7H ----
library(ComplexHeatmap)
net = list(cc@net$`E2+Aducanumab`$weight - cc@net$`E2+Isocon`$weight,
           cc@net$`E3+Aducanumab`$weight - cc@net$`E3+Isocon`$weight,
           cc@net$`E4+Aducanumab`$weight - cc@net$`E4+Isocon`$weight)
mat = Reduce("+",net)/length(net)
color.use = scPalette(ncol(mat))
color.heatmap = c("#2166ac","#b2182b")
names(color.use) = colnames(mat)
color.heatmap.use = colorRamp3(c(min(mat), 0, max(mat)), c(color.heatmap[1], "#f7f7f7", color.heatmap[2]))
colorbar.break = c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",min(mat, na.rm = T)))+1), 0, round(max(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",max(mat, na.rm = T)))+1))
df = data.frame(group = colnames(mat)); rownames(df) = colnames(mat)
col_annotation = HeatmapAnnotation(df = df, col = list(group = color.use),which = "column",
                                   show_legend = FALSE, show_annotation_name = FALSE,
                                   simple_anno_size = grid::unit(0.2, "cm"))
row_annotation = HeatmapAnnotation(df = df, col = list(group = color.use), which = "row",
                                   show_legend = FALSE, show_annotation_name = FALSE,
                                   simple_anno_size = grid::unit(0.2, "cm"))

ha1 = rowAnnotation(Strength = anno_barplot(rowSums(abs(mat)), border = FALSE,gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)
ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(abs(mat)), border = FALSE,gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)
mat[mat == 0] = NA
colnames(mat)[c(16,22,25)] = paste0("***",colnames(mat)[c(16,22,25)],"***")
rownames(mat)[11] = "***Effector-hi TIMs***"
ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", name = "Relative Values",
              bottom_annotation = col_annotation, left_annotation=row_annotation, top_annotation = ha2, right_annotation = ha1,
              cluster_rows = T,cluster_columns = T, row_split = 5, column_split = 4, row_gap = unit(2, "mm"), column_gap = unit(2, "mm"), border=T,
              row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = 8),column_names_gp = gpar(fontsize = 8),
              column_title = 'Targets (Receiver)',column_title_gp = gpar(fontsize = 10),column_names_rot = 45, column_title_side = "bottom",
              row_title = "Sources (Sender)",row_title_gp = gpar(fontsize = 10),row_title_rot = 90, row_dend_side = "right",
              row_labels = gt_render(rownames(mat)), column_labels = gt_render(colnames(mat)),
              heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"),title_position = "leftcenter-rot",
                                          border = NA,
                                          legend_height = unit(20, "mm"),labels_gp = gpar(fontsize = 8),grid_width = unit(2, "mm"))
)
png("interactions_hmap.png",height=6,width=7.54,units="in",res=600)
draw(ht1)
dev.off()

## circos plots ----
weight.max = getMaxWeight(cc.list, attribute = c("idents","count"))
par(mfrow = c(2,3), xpd=TRUE)
for (i in 1:length(cc.list)) {
  netVisual_circle(cc.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(cc.list)[i]))
}

## Fig. S6B ----
num.link = sapply(cc.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax = c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg = list()
for (i in 1:length(cc.list)) {
  gg[[i]] = netAnalysis_signalingRole_scatter(cc.list[[i]], title = names(cc.list)[i], weight.MinMax = weight.MinMax)
  gg[[i]]$layers = c(geom_hline(yintercept = 15, linetype = "dashed"),
                     geom_vline(xintercept = 7.5, linetype = "dashed"),
                     gg[[i]]$layers)
}
plot = patchwork::wrap_plots(plots = gg, guides = "collect") & xlim(0,12) & 
  ylim(0,23) & theme(plot.title = element_text(size=15))
ggsave("cluster_interaction_shifts.png",plot,dpi=600,width=12,height=8)

## look at specific clusters between conditions ----
library(ggmagnify)
# Fig. 7G ----
(netAnalysis_signalingChanges_scatter(cc, idents.use = "Cd8 T Cells", 
                                      comparison = c(1,3), signaling.exclude = "MHC-I") +
    theme(legend.direction = "vertical", legend.background = element_rect(fill=NA)) + 
    ggtitle("Signaling Pathways in Cd8 T Cells,\nE2+Aducanumab vs. E4+Aducanumab") +
    theme(legend.position = "bottom", legend.direction = "vertical") +
    geom_magnify(aes(from = abs(outgoing) <= 0.01 & abs(incoming) <= 0.015),
                 to = c(xmin = 0.3, xmax = 0.65, ymin = 0.035, ymax = 0.245), shadow = T))%>%
  ggsave("cd8_e2e4_pathways.png",.,dpi=600,width=5,height=7)

# Fig. S6C ----
sca1 = (netAnalysis_signalingChanges_scatter(cc, idents.use = "Effector-hi TIMs", 
                                             comparison = c(1,3), signaling.exclude = "MHC-I") +
          theme(legend.direction = "vertical", legend.background = element_rect(fill=NA)) + 
          ggtitle("Signaling Pathways in Effector-hi TIMs,\nE2+Aducanumab vs. E4+Aducanumab"))
sca2 = (netAnalysis_signalingChanges_scatter(cc, idents.use = "Effector-lo TIMs", 
                                             comparison = c(1,3), signaling.exclude = "MHC-I") +
          theme(legend.direction = "vertical", legend.background = element_rect(fill=NA)) + 
          ggtitle("Signaling Pathways in Effector-lo TIMs,\nE2+Aducanumab vs. E4+Aducanumab"))
sca3 = (netAnalysis_signalingChanges_scatter(cc, idents.use = "Immature B Cells", 
                                             comparison = c(1,3), signaling.exclude = "MHC-I") +
          theme(legend.direction = "vertical", legend.background = element_rect(fill=NA)) + 
          ggtitle("Signaling Pathways in Immature B Cells,\nE2+Aducanumab vs. E4+Aducanumab"))
(wrap_plots(list(sca1,sca2,sca3)) + plot_layout(guides = "collect")) %>% 
  ggsave("supp_e2e4_pathways.png",.,dpi=600,width=14,height=4)

cc = computeNetSimilarityPairwise(cc, type = "functional")
cc = netEmbedding(cc, type = "functional")
cc = netClustering(cc, type = "functional", do.parallel = F)
netVisual_embeddingPairwise(cc, type = "functional", label.size = 3.5)

rankSimilarity(cc, type = "functional", comparison2 = c(4,1))
rankSimilarity(cc, type = "functional", comparison2 = c(5,2))
rankSimilarity(cc, type = "functional", comparison2 = c(6,3))

## Fig. 7F ----
rankSimilarity_custom = function(grp1, grp2, tit){
  Y = cc@netP$similarity[["functional"]]$dr[["1-2-3-4-5-6"]]
  group = sub(".*--", "", rownames(Y))
  data1 = Y[group %in% grp1, ]
  data2 = Y[group %in% grp2, ]
  rownames(data1) = sub("--.*", "", rownames(data1))
  rownames(data2) = sub("--.*", "", rownames(data2))
  pathway.show = as.character(intersect(rownames(data1), rownames(data2)))
  data1 = data1[pathway.show, ]
  data2 = data2[pathway.show, ]
  euc.dist = function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
  dist = NULL
  for(i in 1:nrow(data1)) dist[i] = euc.dist(data1[i,],data2[i,])
  df = data.frame(name = pathway.show, dist = dist, row.names = pathway.show)
  df = df[order(df$dist), , drop = F]
  df$name = factor(df$name, levels = as.character(df$name))
  df = tail(df, n=10)
  gg = ggplot(df, aes(x=name,y=dist)) + 
    geom_segment(aes(x=name,xend=name,y=0,yend=dist)) +
    geom_point(size = 5, shape = 21, aes(fill=dist)) +
    theme_classic(base_size = 18) + 
    coord_flip() +
    scale_fill_gradient2(limits = c(0,8), low="blue",mid="white",high="red") +
    labs(x="",y="", fill = "",
         title = tit)
  return(gg)
}

wrap = wrap_plots(list(rankSimilarity_custom("E2+Aducanumab","E2+Isocon","ApoE2") + labs(x="Pathway"),
                       rankSimilarity_custom("E3+Aducanumab","E3+Isocon","ApoE3"),
                       rankSimilarity_custom("E4+Aducanumab","E4+Isocon","ApoE4"))) + plot_layout(guides = "collect")
lab = ggplot(data.frame(l = "Euclidean Distance on Functional Similarity Manifold", x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), size = 6) + theme_void() + coord_cartesian(clip = "off")
(wrap/lab + plot_layout(heights = c(30,1))) %>% 
  ggsave("euclidean_distance.png",.,dpi=600,width=11,height=6.5)

## Fig. 7E ----
temp1 = rankNet(cc, mode = "comparison", stacked = T, do.stat = TRUE, comparison = c(1,3), return.data = T)
temp1 = temp1$signaling.contribution %>% as.data.table
temp2 = rankNet(cc, mode = "comparison", stacked = T, do.stat = TRUE, comparison = c(1,2), return.data = T)
temp2 = temp2$signaling.contribution %>% as.data.table
temp3 = rankNet(cc, mode = "comparison", stacked = T, do.stat = TRUE, comparison = c(2,3), return.data = T)
temp3 = temp3$signaling.contribution %>% as.data.table

temp = rbind(temp1,temp2,temp3) %>% unique
temp[name == "GAS", name := "GAS6"] #replace name for clarity
temp = temp[, .(contmean = mean(contribution.scaled)), by = .(name,group)]
temp$group = factor(temp$group, levels = c("E4+Aducanumab","E3+Aducanumab","E2+Aducanumab"))
temp[, pctdiff := (max(contmean)-min(contmean))/min(contmean), by = name]
temp[, scoresum := sum(contmean), by = name]

flowscore = ggplot(temp[pctdiff>0.75 & scoresum >= 4], aes(x = contmean, y = forcats::fct_reorder(name,scoresum), fill = group)) + 
  geom_bar(stat="identity",position="stack",color="black",lwd=0.4) + 
  theme_classic(base_size = 16) + 
  scale_fill_manual(values = rev(scales::hue_pal()(3)), guide = guide_legend(reverse = TRUE)) +
  scale_x_continuous(expand = c(0,0)) + 
  labs(x="Information Flow Score", y = "Pathway", fill = "Group") +
  theme(legend.position = "bottom", legend.direction = "vertical")

ggsave("information_flow.png",flowscore,dpi=600,width=5,height=8)

