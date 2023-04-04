library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)
library(data.table)
library(dplyr)
library(patchwork)
library(sctransform)
library(ggplot2)
library(Nebulosa)
library(ComplexHeatmap)

mtx_directory = "path/to/mtx/files"

dir_list = list.files(mtx_directory)
adapoe = NULL

for (i in seq_along(dir_list)){
  # Load in data
  tmp = Read10X(paste0(mtx_directory,"/", dir_list[i]))
  if(stringr::str_sub(dir_list[i], start = -3) == "2yr"){
    tmpseu = CreateSeuratObject(counts = tmp)
    tmpseu$age = "2yr"
  } else {
    tmpseu = CreateSeuratObject(counts = tmp$`Gene Expression`)
    tmpseu[["HTO"]] = CreateAssayObject(counts = tmp$`Multiplexing Capture`)
    tmpseu$age = stringr::str_sub(dir_list[i], start = -4)
  }
  tmpseu$orig.ident = dir_list[i]
  tmpseu = RenameCells(object = tmpseu, new.names = paste0(dir_list[i], "_", rownames(tmpseu[[]])))
  # MiQC filtering
  tmpseu[["percent.mt"]] = PercentageFeatureSet(tmpseu, pattern = "^mt-")
  tmpseu = RunMiQC(tmpseu, percent.mt = "percent.mt", nFeature_RNA = "nFeature_RNA", 
                   posterior.cutoff = 0.75, model.slot = "flexmix_model", model.type = "spline")
  tmpseu = subset(tmpseu, miQC.keep == "keep")
  # Merge as output
  if (i == 1){
    adapoe = tmpseu
  } else {
    adapoe = merge(x = adapoe, y = tmpseu)
  }
}

# run SCTransform pipeline and integration ----
adapoe.list = SplitObject(adapoe, split.by = "orig.ident")
for (i in names(adapoe.list)) {
  adapoe.list[[i]] = SCTransform(adapoe.list[[i]], verbose = TRUE)
}
all_genes = lapply(adapoe.list, rownames) %>% Reduce(intersect, .)
adapoe.features = SelectIntegrationFeatures(object.list = adapoe.list, nfeatures = 3000)
adapoe.list = PrepSCTIntegration(object.list = adapoe.list, anchor.features = adapoe.features)
adapoe.anchors = FindIntegrationAnchors(object.list = adapoe.list, normalization.method = "SCT", 
                                        anchor.features = adapoe.features)

adapoe.integrated = IntegrateData(anchorset = adapoe.anchors, normalization.method = "SCT", features.to.integrate = adapoe.features)

# clustering ----
adapoe.integrated = RunPCA(adapoe.integrated, verbose = FALSE)
ElbowPlot(adapoe.integrated, ndims = 50)
adapoe.integrated = RunUMAP(adapoe.integrated, reduction = "pca", dims = 1:45)
adapoe.integrated = FindNeighbors(adapoe.integrated, dims = 1:45, verbose = TRUE)
adapoe.integrated = FindClusters(adapoe.integrated, verbose = TRUE, resolution = 1.4)
adapoe.integrated$genotype = stringr::str_sub(adapoe.integrated$orig.ident, 1, 2)
saveRDS(adapoe.integrated, "mtx_directory/adapoe.rds")
adapoe = adapoe.integrated

# annotation time! ----
markers = FindAllMarkers(adapoe, assay = "RNA", min.pct = 0.5, logfc.threshold = 0.5, only.pos = TRUE)
fwrite(markers, "mtx_directory/markers.txt", sep = "\t")
idents = c(rep("Microglia",8),"Cd8 T Cells",rep("Microglia",3),"Mature B Cells",rep("Microglia",2),"Low Quality","Monocytes","Macrophages",rep("Microglia",2),
           "Late Neutrophils","DP T Cells","Microglia","NK Cells","Activated NK Cells","Dendritic Cells","F13a1+ Monocytes","Immature B Cells","Microglia","Early Neutrophils","Microglia",
           "Microglia","ILC2s","GD T Cells","Microglia","Granulocyte-Monocyte Progenitors","Migratory DCs") #note: mig-DCs: https://elifesciences.org/articles/56890
names(idents) = levels(adapoe)
adapoe = RenameIdents(adapoe, idents)
adapoe$clust_id = Idents(adapoe) %>% as.character
DimPlot(adapoe, label = TRUE, repel = TRUE) + NoLegend()
adapoe$microglia = adapoe$clust_id == "Microglia"
adapoe = subset(adapoe, clust_id != "Low Quality")
saveRDS(adapoe, "mtx_directory/adapoe.rds")

# subcluster microglia ----
seurat = subset(adapoe, microglia == TRUE)
seurat = SCTransform(seurat, method = "glmGamPoi", vars.to.regress = c("percent.mt","orig.ident"))
seurat = RunPCA(seurat, npcs = 50)
ElbowPlot(seurat, ndims = 50)
seurat = RunUMAP(seurat, dims = 1:46, verbose = TRUE)
seurat = FindNeighbors(seurat, dims = 1:46, verbose = TRUE)
seurat = FindClusters(seurat, verbose = TRUE, resolution = 0.6)
DimPlot(seurat, label = TRUE) + NoLegend()
markers = FindAllMarkers(seurat, assay = "RNA", min.pct = 0.1, logfc.threshold = 0.1, only.pos = TRUE)
fwrite(markers, "mtx_directory/mglia_markers.txt", sep = "\t")
mglia.ident = c("Homeostatic Microglia","TIMs","Poised Homeostatic Microglia","Hspa1+ Stressed Microglia","Lars2-mid Homeostatic Microglia",
                "Ier2/5+ Inflammatory Microglia","DAM-2","Ccl3/4+ Inflammatory Microglia","Lars2-hi Homeostatic Microglia","Lars2-mid Homeostatic Microglia",
                "Rgs1-hi Homeostatic Microglia","DAM-1","Adamts1+ Inflammatory Microglia","Interferon Induced Microglia","Effector-hi TIMs",
                "Bri3-Negative Homeostatic Microglia","MHCII+ Microglia","Serpine1+ TIMs","Hspb1+ Stressed Microglia","Cd74+ Microglia","Cycling Microglia")
names(mglia.ident) = levels(seurat)
seurat = RenameIdents(seurat, mglia.ident)
my_mglia_levels = c("Homeostatic Microglia","Lars2-mid Homeostatic Microglia","Lars2-hi Homeostatic Microglia","Bri3-Negative Homeostatic Microglia", "Rgs1-hi Homeostatic Microglia",
                    "Hspa1+ Stressed Microglia","Hspb1+ Stressed Microglia","Poised Homeostatic Microglia","DAM-1","DAM-2","Interferon Induced Microglia","Ier2/5+ Inflammatory Microglia",
                    "Adamts1+ Inflammatory Microglia","Ccl3/4+ Inflammatory Microglia","TIMs","Effector-hi TIMs","Serpine1+ TIMs",
                    "MHCII+ Microglia","Cd74+ Microglia","Cycling Microglia")
levels(seurat) = my_mglia_levels
seurat$mglia_ident = Idents(seurat)
saveRDS(seurat, "mtx_directory/mglia_only.rds")

# paste subclust info back into main seurat structure ----
adapoe$subclust_id = as.character(adapoe$clust_id)
adapoe$subclust_id[Cells(seurat)] = as.character(Idents(seurat)) 
new.levels = c(levels(seurat$mglia_ident), "Early Neutrophils","Late Neutrophils","Immature B Cells","Mature B Cells","Monocytes","F13a1+ Monocytes",
               "Dendritic Cells","Migratory DCs","Macrophages","Cd8 T Cells","DP T Cells","GD T Cells","NK Cells","Activated NK Cells","ILC2s","Granulocyte-Monocyte Progenitors")
adapoe$subclust_id = factor(adapoe$subclust_id, levels = new.levels)
VlnPlot(adapoe, "Serpine1", group.by = 'subclust_id') + NoLegend()
Idents(adapoe) = adapoe$clust_id;DimPlot(adapoe, label = TRUE, repel = TRUE) + NoLegend()
Idents(adapoe) = adapoe$subclust_id;DimPlot(adapoe, label = TRUE, repel = TRUE) + NoLegend()
DimPlot(seurat, label = TRUE, repel = TRUE) + NoLegend()

DimPlot(mglia, group.by = "age") %>% 
  ggsave("mtx_directory/mglia_umap_byage.png",.,
         dpi=600,width=2400,height=2400,units="px")
DimPlot(mglia, group.by = "genotype") %>% 
  ggsave("mtx_directory/mglia_umap_bygeno.png",.,
         dpi=600,width=2400,height=2400,units="px")

saveRDS(adapoe, "mtx_directory/adapoe.rds")

# run CIARA to identify rare cell types ----
library(CIARA)
mat = adapoe@assays$RNA@counts
coordinate_umap = adapoe@reductions$umap@cell.embeddings %>% as.data.frame
ciara = cluster_analysis_integrate_rare(mat, "adapoe",0.1,5,30)
norm_ciara = as.matrix(Seurat::GetAssayData(ciara, slot = "data", assay = "RNA"))
orig_clusters = Idents(adapoe)
knn_ciara = as.matrix(ciara@graphs$RNA_nn)

background = get_background_full(norm_ciara, threshold = 1, n_cells_low = 3, 
                                 n_cells_high = 20)
result = CIARA(norm_ciara, knn_ciara, background, cores_number = 8,  
               p_value = 0.001, odds_ratio = 2, local_region = 1, 
               approximation = FALSE) #note: takes ~2d to run!
saveRDS(result, "ciara_res.Rda")
ciara_genes = row.names(result)[result[, 1]<1]
ciara_genes_top = row.names(result)[order(as.numeric(result[,1]))]
norm_ciara_filt = norm_ciara[ciara_genes,]
norm_ciara_mglia_filt = norm_ciara_mglia[ciara_genes,]
## use AUCell to define new cluster from relevant genes
library(AUCell)
geneset = list(ciara = ciara_genes_top[c(1,2,4,5,7,8,9,11,12,13,14,15,16)])
cell_ranking = AUCell_buildRankings(norm_ciara_mglia, nCores = 4)
ciara_sig = AUCell_calcAUC(geneSets = geneset, rankings = cell_ranking, nCores = 4, aucMaxRank = ceiling(0.05*nrow(cell_ranking)))
cells_assignment = AUCell_exploreThresholds(ciara_sig, plotHist=TRUE, assign=TRUE)
custom_thresh = names(which(getAUC(ciara_sig)["ciara",]>0.15))
coordinate_umap_mglia[custom_thresh,]
length(custom_thresh)
## export aucell scores
ciara_aucell = ciara_sig@assays@data$AUC[1,]
## paste new cluster idents into mglia...
new_labels = data.table(cell = Cells(mglia), ident = as.character(mglia$mglia_ident))
new_labels[cell %in% custom_thresh, ident := "Il34+ Microglia"]
levels = levels(mglia$mglia_ident)
levels = append(levels, "Il34+ Microglia", after = 11)
new_labels$ident = factor(new_labels$ident, levels = levels)
mglia$mglia_ident = new_labels$ident
Idents(mglia) = mglia$mglia_ident
## ...and into adapoe
adapoe$subclust_id = as.character(adapoe$clust_id)
adapoe$subclust_id[Cells(mglia)] = as.character(Idents(mglia)) 
new.levels = c(levels(mglia$mglia_ident), "Early Neutrophils","Late Neutrophils","Immature B Cells","Mature B Cells","Monocytes","F13a1+ Monocytes",
               "Dendritic Cells","Migratory DCs","Macrophages","Cd8 T Cells","DP T Cells","GD T Cells","NK Cells","Activated NK Cells","ILC2s","Granulocyte-Monocyte Progenitors")
adapoe$subclust_id = factor(adapoe$subclust_id, levels = new.levels)
Idents(adapoe) = adapoe$subclust_id
## save both back, all done!
saveRDS(adapoe, "mtx_directory/adapoe.rds")
saveRDS(mglia, "mtx_directory/mglia_only.rds")  

# factorize stuff ----
adapoe$orig.ident = factor(adapoe$orig.ident, levels = c("E2_10wk","E3_10wk","E4_10wk","E2_20wk","E3_20wk","E4_20wk","E3_2yr","E4_2yr"))
seurat$orig.ident = factor(seurat$orig.ident, levels = c("E2_10wk","E3_10wk","E4_10wk","E2_20wk","E3_20wk","E4_20wk","E3_2yr","E4_2yr"))
adapoe$age = factor(adapoe$age, levels = c("10wk","20wk","2yr"))
seurat$age = factor(seurat$age, levels = c("10wk","20wk","2yr"))
adapoe$genotype = factor(adapoe$genotype, levels = c("E2","E3","E4"))
seurat$genotype = factor(seurat$genotype, levels = c("E2","E3","E4"))

md = adapoe@meta.data %>% as.data.table
setkey(md, subclust_id, age)
dat = md[microglia == TRUE][CJ(subclust_id, age, unique = TRUE),.N,by=.EACHI]
dat[, frac := N/sum(N), by = age]
(ggplot() + 
  geom_bar(data = dat, aes(x = subclust_id, y = frac, fill = age), color = "black", stat = "identity", position = "dodge") + 
  labs(x = "", y= "Fraction of Cells", fill = "Age") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),plot.margin = margin(5,5,5,50),
        text = element_text(size=rel(5)), legend.text = element_text(size=rel(3)))) %>% 
  ggsave("mtx_directory/fracs_mglia_byage.png",.,
         dpi=600,width=16000,height=6000,units="px")

# Figure 1B ----
(DimPlot(adapoe, label = TRUE, repel = TRUE, label.size = 8, group.by = "clust_id") + NoLegend() + 
    theme(plot.title = element_blank(),
          axis.text = element_text(size=20),
          axis.title = element_text(size=30),
          legend.text = element_text(size=22.5))) %>% 
  ggsave("adapoe_umap.png",.,width=14.812,height=17.774,dpi=600)

# Figure 1C ----
(DimPlot(mglia, group.by = "mglia_ident") +
    theme(plot.title = element_blank(),
          axis.text = element_text(size=20),
          axis.title = element_text(size=30),
          legend.text = element_text(size=22.5),
          legend.spacing.y = unit(0.1, 'in')) +
    guides(color=guide_legend(ncol=1,override.aes = list(size = 5),byrow = TRUE))) %>% 
  ggsave("mglia_umap.png",.,width=14.06,height=10.935,dpi=600)

# Figure 1D ----
geno = DimPlot(mglia, group.by = "genotype") + theme(plot.title = element_blank(),
                                                     axis.text = element_text(size=20),
                                                     axis.title = element_text(size=30),
                                                     legend.text = element_text(size=22.5),
                                                     legend.spacing.y = unit(0.1, 'in')) + guides(color=guide_legend(override.aes = list(size = 5),byrow=TRUE))
mglia$age_forplots = ifelse(mglia$age == "2yr", "96wk", as.character(mglia$age))
mglia$age_forplots = factor(mglia$age_forplots, levels = c("10wk","20wk","96wk"))
age = DimPlot(mglia, group.by = "age_forplots") + theme(plot.title = element_blank(),
                                                        axis.text = element_text(size=20),
                                                        axis.title = element_text(size=30),
                                                        legend.text = element_text(size=22.5),
                                                        legend.spacing.y = unit(0.1, 'in')) + guides(color=guide_legend(override.aes = list(size = 5),byrow=TRUE))
(geno|age) %>% ggsave("mglia_umap_faceted.png",.,width=12.167,height=6.083,dpi=600)

# Figure 1E ----
library(limma)
library(edgeR)
library(ggrepel)

options(ggrepel.max.overlaps = Inf)
mglia$grouping = ifelse(mglia$mglia_ident %in% c("TIMs","Serpine1+ TIMs", "Effector-hi TIMs"), "TIMs", "Other")
## calculate DEG statistics
expr = mglia@assays$RNA@counts
d0 = DGEList(expr)
d0 = calcNormFactors(d0)
d = d0[-which(apply(cpm(d0),1,max)<1),]
mm = model.matrix(~0 + grouping, data = mglia@meta.data)
rm(d0) # for memory clearance
gc()
y = voom(d, mm, plot = F)
fit = lmFit(y,mm)
contr = makeContrasts(groupingTIMs - groupingOther, levels=colnames(coef(fit)))
tmp = contrasts.fit(fit, contr)
tmp = eBayes(tmp, trend = TRUE)
saveRDS(tmp, "limma_tim_v_other_mglia_ebayes.rds")
## plot using B statistic
top.table = topTable(tmp, sort.by = "B", n = Inf)
top.table$logB = ifelse(top.table$B < 1, 0, log2(top.table$B))
top.table$gene = rownames(top.table)
top.table = as.data.table(top.table)
b_cutoff = log(100) # this is equivalent to 100-to-1 log-odds of differential expression
logfc_cutoff = 0.5
top.table[, in_zone := ifelse(abs(logFC) >= logfc_cutoff & logB >= b_cutoff, TRUE, FALSE)]
top.table$label = FALSE
top.table[logFC > 0, label := ifelse(logFC >= .SD$logFC[kit::topn(.SD$logFC, 15, decreasing = T)[15]], TRUE, FALSE)]
top.table[logFC < 0, label := ifelse(-logFC >= -.SD$logFC[kit::topn(-.SD$logFC, 15, decreasing = T)[15]], TRUE, FALSE)]
top.table[, group := ifelse(abs(logFC) >= logfc_cutoff & logB >= b_cutoff, "Hit", "Filter")]
top.table$group = factor(top.table$group, levels = c("Hit","Filter"))
volcano = ggplot(top.table, aes(x = logFC, y = logB)) +
  geom_point(aes(color = top.table$group)) +
  scale_color_manual(values = c("blue","gray")) +
  geom_text_repel(label = ifelse(top.table$label==TRUE,top.table$gene,''),
                  size = 8, force = 10, box.padding = 0.2, point.padding = 0.2, segment.size = 0.1, min.segment.length = 0.1, max.iter = 1e7) +
  geom_vline(xintercept = -logfc_cutoff) + 
  geom_vline(xintercept = logfc_cutoff) + 
  geom_hline(yintercept = b_cutoff) +
  geom_text(data = data.frame(xpos = c(-(logfc_cutoff+0.2), logfc_cutoff+0.2), 
                              ypos = c(b_cutoff, b_cutoff), hjustvar = c(1,0), vjustvar = c(-1,-1), 
                              annotateText = c(paste0("Non-TIM Microglia ",sprintf('\U2190')),paste0(sprintf('\U2192')," TIMs"))),
            aes(x=xpos, y=ypos, hjust=hjustvar, vjust=vjustvar, label=annotateText), size = 8) +
  theme_bw() + 
  labs(x = "Log-2 Fold Change", y = "Log-2 B-Value") + 
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=30)) +
  NoLegend()
ggsave("limma_volcano_tim_v_other_mglia.png",volcano,dpi=600,width=11.726,height=11.726,units="in")

# Figure 1F ----
md = adapoe@meta.data %>% as.data.table
setkey(md, subclust_id, orig.ident)
dat = md[microglia == TRUE][CJ(subclust_id, orig.ident, unique = TRUE),.N,by=.EACHI]
dat[, frac := N/sum(N), by = orig.ident]
dat = dcast.data.table(dat, orig.ident ~ subclust_id, value.var = "frac")
datm = as.matrix(dat[,-1])
datm.filtered = datm[,c(1,8,9,10,16,17,18)]
(datm.filtered %>% as.data.table(keep.rownames = T) %>% 
    melt.data.table(id.vars = "rn") %>% 
    mutate(rn = replace(rn, rn == "E3_2yr", "E3_96wk")) %>% 
    mutate(rn = replace(rn, rn == "E4_2yr", "E4_96wk")) %>% 
    mutate(rn = factor(rn, levels = c(levels(adapoe$orig.ident)[1:6],"E3_96wk","E4_96wk"))) %>% 
    mutate(variable = factor(stringr::str_wrap(variable,10),
                             levels=stringr::str_wrap(levels(adapoe)[c(1,8,9,10,16,17,18)],10))) %>% 
    ggplot(aes(x=variable,y=100*value,fill=rn)) + 
    geom_bar(stat = "identity", position = "dodge",color="black", linewidth=1) + 
    scale_fill_manual(values = c("lightskyblue","chartreuse","coral","royalblue","springgreen3","tomato3","seagreen4","red4")) +
    labs(x = "",y="Frequency of Microglia (%)",fill="Sample") + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = 25),
          axis.text.y = element_text(size=20),
          axis.title = element_text(size=30),
          legend.text = element_text(size=17.5),
          legend.title = element_text(size=20),
          legend.spacing.y = unit(0.05, 'in')) + 
    guides(fill=guide_legend(ncol=1,override.aes = list(size = 8),byrow = TRUE))) %>% 
  ggsave("mglia_freqs_filtered_barplot.png",.,width=17.256,height=7.888,dpi=600)

# Prep files for Figure 1G (generated in python) ----
seurat = subset(adapoe, orig.ident == "E4_2yr") # repeat this for E3_2yr and the 20wk samples.
cells = gsub(".*_(.+)-.*","\\1",Cells(seurat)) %>% as.data.frame
embeddings = seurat@reductions$umap@cell.embeddings
rownames(embeddings) = gsub(".*_(.+)-.*","\\1",rownames(embeddings))
clusters = seurat$mglia_ident %>% as.data.frame
rownames(clusters) = gsub(".*_(.+)-.*","\\1",rownames(clusters))
colnames(clusters) = "Cluster"
fwrite(cells, "cells.csv", row.names = TRUE)
fwrite(embeddings, "embeddings.csv", row.names = TRUE)
fwrite(clusters, "clusters.csv", row.names = TRUE)

# Figure S1A ----
genes_toplot = c("P2ry12","Tmem119","Lars2","Bri3","Rgs1","Hspa1a","Hspb1","Numb","Trem2","Cst7","Apoe",
                 "Ifit3","Il34","Ier2","Adamts1","Ccl3","Fos","Arhgap45","Serpine1","H2-Aa","Cd74",
                 "Stmn1","Camp","Il1b","Rag1","Cd79a","Ace","F13a1","Cd209a","Ccr7","Adgre1",
                 "Cd8a","Cd4","Trdv4","Klra9","Gzma","Gata3","Fcnb")
(DotPlot(adapoe,features=genes_toplot,assay="SCT", dot.scale = 15) + 
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, size = 40),
          axis.text.y = element_text(size=40), axis.title = element_text(size=40),
          legend.title = element_text(size=40), legend.text = element_text(size=35),legend.spacing.y = unit(0.1, 'in')) +
    guides(size = guide_legend(byrow = TRUE, order = 1), color = guide_colorbar(barheight = 12, order = 2)) + 
    labs(color = "Average Expression", size = "Percent Expressed")
) %>% ggsave("markers_dotplot.pdf",.,height=20,width=40,dpi=600)

# Figure S1C ----
library(ggforce)

mglia$ciara_val = ciara_aucell
center.coords = mglia@reductions$umap@cell.embeddings[rownames(mglia@reductions$umap@cell.embeddings) %in% 
                                                        Cells(subset(mglia, subset = mglia_ident == "Il34+ Microglia")),] %>% colMeans

(FeaturePlot(mglia, "ciara_val", pt.size = 8) + 
    ggtitle("AUCell Scores for CIARA-Derived Il34+ Microglia Geneset") +
    geom_circle(aes(x0 = center.coords[1], y0 = center.coords[2], r = 0.5), color = "red",
                linewidth = 3, inherit.aes = FALSE) +
    theme(plot.title = element_text(size=55),
          axis.text = element_text(size=30),
          axis.title = element_text(size=40),
          legend.text = element_text(size=35))+
    guides(color = guide_colorbar(barheight=15))
) %>% ggsave("ciara_aucell_vals.png",.,width=26.139,height=18.432,dpi=600)

# Figure S1D ----
library(ranger)

## generate predictions
counts = Matrix::t(mglia@assays$RNA@counts)
idents = mglia$mglia_ident %>% as.character
## function adapted from https://github.com/vilasmenon/Microglia_Olah_et_al_2020/blob/master/code_for_submission_functions.R
cell_by_cell_prediction=function(dataset,clusterids,crossval=4,iterations=10,filename="rf_pred.rda") {
  allf=unique(clusterids)
  predlist=list()
  for (ii in 1:(length(allf)-1)) {
    testcols1=rownames(dataset)[which(clusterids==allf[ii])]
    for (jj in (ii+1):length(allf)) {
      testcols2=rownames(dataset)[which(clusterids==allf[jj])]
      outmat1=matrix(0,nrow=iterations,ncol=length(testcols1))
      colnames(outmat1)=testcols1
      outmat2=matrix(0,nrow=iterations,ncol=length(testcols2))
      colnames(outmat2)=testcols2
      numtrain=min(round(length(testcols1)*3/4),round(length(testcols2)*3/4))
      alldone1=c()
      alldone2=c()
      set.seed(ii*iterations+jj)
      for (kk in 1:iterations) {
        sampids1=sample((1:length(testcols1))%%crossval)
        sampids2=sample((1:length(testcols2))%%crossval)
        for (mm in unique(sampids1)) {
          trainset=c(sample(testcols1[sampids1!=mm],min(numtrain,length(which(sampids1!=mm)))),sample(testcols2[sampids2!=mm],min(numtrain,length(which(sampids2!=mm)))))
          testset=c(testcols1[sampids1==mm],testcols2[sampids2==mm])
          ttt=as.factor(rep(c(allf[ii],allf[jj]),times=c(min(numtrain,length(which(sampids1!=mm))),min(numtrain,length(which(sampids2!=mm))))))
          #predval=randomForest(x = as.matrix(dataset[trainset,]), y = ttt)
          predval2 = ranger(x = dataset[trainset,], y = ttt, num.trees = 100)
          #outpred=predict(predval,as.matrix(dataset[testset,]))
          outpred2 = predict(predval2, dataset[testset,])
          outpred2 = outpred2$predictions
          #names(outpred) = testset
          names(outpred2)=testset
          outmat1[kk,testcols1[sampids1==mm]]=as.character(outpred2[testcols1[sampids1==mm]])
          outmat2[kk,testcols2[sampids2==mm]]=as.character(outpred2[testcols2[sampids2==mm]])
        }
        print(paste0(Sys.time(), ": Clust ", ii, " (",allf[ii],")", " by Clust ", jj," (",allf[jj],")", ", Iteration ", kk,"/",iterations))
      }
      nam=paste0(allf[ii],"-",allf[jj])
      predlist[[nam]]=list()
      predlist[[nam]][[1]]=outmat1
      predlist[[nam]][[2]]=outmat2
      save(predlist,file=filename)
    }
  }
}

cell_predictions=cell_by_cell_prediction(dataset=counts,clusterids=idents,crossval=4,iterations=25,
                                         filename="rf_predictions_all_clusters.rda")

## plot
load("rf_predictions_all_clusters.rda")
idents = unique(idents)

conf_mat = matrix(nrow = length(idents), ncol = length(idents))
rownames(conf_mat) = idents
colnames(conf_mat) = idents
minsum = function(tab){
  if(length(tab) == 1){
    return(0)
  } else {
    return(min(tab)/sum(tab))
  }
}

mat.calc = function(mat){
  most.common.preds = c()
  for(col in 1:ncol(mat)){
    most.common.preds = append(most.common.preds, names(sort(table(mat[,col]), decreasing = TRUE)[1]))
  }
  return(minsum(table(most.common.preds)))
}

#take care of row 1 separately, element is simply column - 1
for(col in 2:length(idents)){
  #conf_mat[1,col] = minsum(table(predlist[[col-1]][[1]]))
  conf_mat[1,col] = mat.calc(predlist[[col-1]][[1]])
}
#all rows beyond row 1, element is column - 1 + upper triangular from 15 by row count
for(row in 2:length(idents)){
  for(col in (row+1):length(idents)){
    pl.el = (col - 1) + sum((length(idents)-1)-(1:(row-1)))
    #conf_mat[row,col] = minsum(table(predlist[[pl.el]][[1]]))
    conf_mat[row,col] = mat.calc(predlist[[pl.el]][[1]])
  }
}
#now lower diagonal, same idea but switching row/col indices and pulling from [[2]] of predlist
for(row in 2:length(idents)){
  #conf_mat[row,1] = minsum(table(predlist[[row-1]][[2]]))
  conf_mat[row,1] = mat.calc(predlist[[row-1]][[2]])
}
for(col in 2:length(idents)){
  for(row in (col+1):length(idents)){
    pl.el = (row - 1) + sum((length(idents)-1)-(1:(col-1)))
    #conf_mat[row,col] = minsum(table(predlist[[pl.el]][[2]]))
    conf_mat[row,col] = mat.calc(predlist[[pl.el]][[2]])
  }
}
diag(conf_mat) = NA
conf_mat = conf_mat[match(levels(mglia$mglia_ident),rownames(conf_mat)),]
conf_mat = conf_mat[,match(levels(mglia$mglia_ident),colnames(conf_mat))]

dev.off()
png("clustering_robustness/25_iterations_100_tree.png",res=600,width=10,height=10,units="in")
ht = Heatmap(conf_mat, cluster_rows = FALSE, cluster_columns = FALSE, rect_gp = gpar(col = "black", lwd = 0.5),
             column_split = factor(c(rep("Homeostatic",8),rep("DAMs",2),rep("Inflammatory",5),
                                     rep("TIMs",3),rep("Other",3)),
                                   levels = c("Homeostatic","DAMs","Inflammatory","TIMs","Other")),
             row_split = factor(c(rep("Homeostatic",8),rep("DAMs",2),rep("Inflammatory",5),
                                  rep("TIMs",3),rep("Other",3)),
                                levels = c("Homeostatic","DAMs","Inflammatory","TIMs","Other")),
             column_gap = unit(2, "mm"), row_gap = unit(2, "mm"),
             cell_fun = function(j,i,x,y,width,height,fill){
               if(i==j){
                 grid.text("X",x=x,y=y)
               }
             },
             row_title = "True Assignment", column_title = "Predicted Assignment",
             heatmap_legend_param = list(title = "Confusion", 
                                         title_position = "topcenter", direction = "horizontal"))
draw(ht, padding = unit(c(15,2,2,10),"mm"), heatmap_legend_side = "top")
dev.off()

# Figure S1E ----
library(rstatix)
library(ggpubr)
md$istim = ifelse(md$subclust_id %in% c("TIMs","Effector-hi TIMs","Serpine1+ TIMs"),"TIM","Non-TIM")
md$istim = factor(md$istim, levels = c("TIM","Non-TIM"))
stat.test = md %>% wilcox_test(miQC.probability ~ istim) %>% add_significance() %>% add_xy_position(x = "istim")
(ggplot(md, aes(x=istim,y=miQC.probability)) +
    geom_boxplot(aes(fill=istim), outlier.shape = NA, size = 2) + 
    stat_pvalue_manual(stat.test %>% mutate(y.position=0.45), label = "p = {p}", vjust = -0.25, size = 15) +
    theme_bw() + 
    ylim(0,0.5) +
    theme(legend.position = "none", axis.text = element_text(size=40), axis.title = element_text(size=45)) +
    labs(x="",y="miQC Score")) %>% 
  ggsave("tim_miqc_nooutlier.png",.,dpi=600,height=20,width=20)
