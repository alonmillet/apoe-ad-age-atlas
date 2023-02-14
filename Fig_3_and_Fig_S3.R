library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(dplyr)
library(data.table)

# Generate annotated Signac structure ----
counts = Read10X_h5("filtered_feature_bc_matrix.h5") # multiome h5 output from CellRanger
fragpath = "atac_fragments.tsv.gz"

# get gene annotations for mouse genome
annotation = GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) = "UCSC"

# create a Seurat object containing the RNA data
seu = CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA")

# create ATAC assay and add it to the object
seu[["ATAC"]] = CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation)

# QC
DefaultAssay(seu) = "ATAC"
seu = NucleosomeSignal(seu)
seu = TSSEnrichment(seu)
VlnPlot(
  object = seu,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0)

# pre-filtering: 5335 barcodes
seu_filt = subset(
  x = seu,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 10000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 100 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)
# post-filtering: 5081 barcodes

# call peaks
peaks = CallPeaks(seu_filt, macs2.path = "path/to/your/macs2")
# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks = keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks = subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)
# quantify counts in each peak
macs2_counts = FeatureMatrix(
  fragments = Fragments(seu_filt),
  features = peaks,
  cells = colnames(seu_filt))
# create a new assay using the MACS2 peak set and add it to the Seurat object
seu_filt[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation)

# processing data
DefaultAssay(seu_filt) = "RNA"
seu_filt = SCTransform(seu_filt)
seu_filt = RunPCA(seu_filt)
DefaultAssay(seu_filt) = "peaks"
seu_filt = FindTopFeatures(seu_filt, min.cutoff = 5)
seu_filt = RunTFIDF(seu_filt)
seu_filt = RunSVD(seu_filt)

# unified projection
seu_filt = FindMultiModalNeighbors(
  object = seu_filt,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40), #remove first component of lsi, as it captures technical variation (can see in depthcor(seu_filt))
  modality.weight.name = "RNA.weight",
  verbose = TRUE)
seu_filt = RunUMAP(
  object = seu_filt,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE)

# add cellranger metadata
md = fread("per_barcode_metrics.csv") # from CellRanger
og_md = seu_filt@meta.data %>% as.data.table(keep.rownames = TRUE)
md = merge.data.table(og_md, md, by.x = "rn", by.y = "barcode")
md = setDF(md)
rownames(md) = md$rn
md$rn = NULL
seu_filt@meta.data = md

# calculate ATAC QC stats for peaks
seu_filt = NucleosomeSignal(seu_filt)
seu_filt = TSSEnrichment(seu_filt, fast = FALSE) #run false for TSSPlot functionality
seu_filt$high.tss = ifelse(seu_filt$TSS.enrichment > 4, 'High', 'Low')
TSSPlot(seu_filt, group.by = 'high.tss') + NoLegend()
seu_filt$nucleosome_group = ifelse(seu_filt$nucleosome_signal > 1.5, 'NS > 1.5', 'NS < 1.5')
FragmentHistogram(seu_filt, group.by = "nucleosome_group", region = "chr1-1-100000000")

# clustering
DefaultAssay(seu_filt) = "RNA"
seu_filt = FindNeighbors(seu_filt, dims = 1:30)
seu_filt = FindClusters(seu_filt, resolution = 1.2, algorithm = 3, graph.name = "wsnn")
DimPlot(seu_filt, label = TRUE, repel = TRUE, reduction = "umap")

saveRDS(seu_filt, "multiome.rds")

test = FindMarkers(seu_filt, ident.1 = 0, ident.2 = 2, test.use = "LR", latent.vars = "atac_peak_region_fragments")
CoveragePlot(seu_filt, region = rownames(test)[2], 
             features = ClosestFeature(seu_filt, rownames(test)[2])$gene_name, 
             expression.assay = "SCT", extend.upstream = 30000, extend.downstream = 30000)

# link peaks
DefaultAssay(seu_filt) = "peaks"
seu_filt = RegionStats(seu_filt, genome = BSgenome.Mmusculus.UCSC.mm10)
seu_filt = LinkPeaks(seu_filt, peak.assay = "peaks", expression.assay = "SCT") # this takes ~4 hours!

idents.mglia = c(0,1,2,3,5,10,11,22)
CoveragePlot(
  object = seu_filt,
  region = "Olfm4",
  features = "Olfm4",
  expression.assay = "SCT",
  idents = idents.mglia,
  extend.upstream = 10000,
  extend.downstream = 10000)

new_idents = c("Homeostatic Microglia","Arhgap15-hi Homeostatic Microglia","mt-Enriched Microglia","DAM-1","DAM-2","TIMs","Siglech-hi Microglia",
               "Interferon Induced Microglia","Monocytes","F13a1+ Monocytes","Macrophages","Early Neutrophils","Inflammatory Neutrophils",
               "B Cells 1","B Cells 2","B Cells 3","B Cells 4","IgM+ B Cells","Naive CD4s","Treg CD4s","Tem CD4s","Trm CD4s","Astrocytes")
names(new_idents) = c(0,1,5,3,2,10,11,22,17,12,15,7,9,8,13,6,18,21,4,14,16,20,19)
seu_filt = RenameIdents(seu_filt, new_idents)
DimPlot(seu_filt, label = TRUE, repel = TRUE) + NoLegend()
seu_filt$clust_ident = Idents(seu_filt)
seu_filt$microglia = ifelse(seu_filt$clust_ident %in% new_idents[1:8],TRUE,FALSE)

# recalculate peaks by cluster.
peaks = CallPeaks(seu_filt, macs2.path = "path/to/your/macs2", group.by = "clust_ident")
# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks = keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks = subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)
# quantify counts in each peak
macs2_counts = FeatureMatrix(
  fragments = Fragments(seu_filt),
  features = peaks,
  cells = colnames(seu_filt))
# create a new assay using the MACS2 peak set and add it to the Seurat object
seu_filt[["peaks_perclust"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation)
# relink peaks
DefaultAssay(seu_filt) = "peaks_perclust"
seu_filt = RegionStats(seu_filt, genome = BSgenome.Mmusculus.UCSC.mm10)
seu_filt = LinkPeaks(seu_filt, peak.assay = "peaks_perclust", expression.assay = "SCT") # this takes ~4 hours!
saveRDS(seu_filt, "multiome.rds")

# Export to h5ad for MIRA analysis ----
library(SeuratDisk)
seu_filt = readRDS("multiome.rds")
DefaultAssay(seu_filt) = "RNA"
seu_filt_RNA = DietSeurat(seu_filt,counts = TRUE,data = TRUE,scale.data = TRUE,
                          features = NULL,assays = "RNA",dimreducs = NULL,graphs = NULL,misc = TRUE)
DefaultAssay(seu_filt) = "peaks_perclust"
seu_filt_peak = DietSeurat(seu_filt,counts = TRUE,data = TRUE,scale.data = TRUE,
                           features = NULL,assays = "peaks_perclust",dimreducs = NULL,graphs = NULL,misc = TRUE)
SaveH5Seurat(seu_filt_RNA, "multiome_RNA.h5Seurat")
SaveH5Seurat(seu_filt_peak, "multiome_ATAC.h5Seurat")
Convert("multiome_RNA.h5Seurat", dest = "h5ad")
Convert("multiome_ATAC.h5Seurat", dest = "h5ad")
# also convert wholesale
SaveH5Seurat(seu_filt, "multiome.h5Seurat")
Convert("multiome.h5Seurat", dest = "h5ad")

# Ran MIRA here; see accompanying code

# post-MIRA, importing back in our joint UMAP generated from topic modeling
library(SeuratDisk)
library(Seurat)
library(Signac)
library(dplyr)
library(data.table)
library(presto)

Convert("multiome_RNA_processed.h5ad", dest = "h5Seurat", overwrite = TRUE) # output from MIRA
mira = LoadH5Seurat("multiome_RNA_processed.h5seurat", assay = "RNA")
seu_filt = readRDS("multiome.rds")

embds = mira@reductions$umap@cell.embeddings
seu_filt@reductions[["mira_umap"]] = CreateDimReducObject(embeddings = embds, key = "MIRAUMAP_", assay = DefaultAssay(seu_filt))
# also adding our topics metadata to our original seurat structure
mira_md = mira@meta.data
topics_md = mira_md[,which(grepl("topic",colnames(mira_md),fixed=TRUE))]
seu_filt = AddMetaData(seu_filt, topics_md)
# now let's annotate our clusters in this space!
seu_filt$mira_clusts = mira@meta.data$leiden
DimPlot(seu_filt, reduction = "mira_umap", group.by = "mira_clusts", repel = TRUE, label = TRUE) + NoLegend()
Idents(seu_filt) = seu_filt$mira_clusts
markers = FindAllMarkers(seu_filt, assay = "SCT", min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE, test.use = "LR") %>% 
  as.data.table(keep.rownames = TRUE)

DefaultAssay(seu_filt) = "peaks_perclust"
seu_filt$temp_subset = ifelse(seu_filt$mira_clusts %in% c(1,5,9), "south", ifelse(seu_filt$mira_clusts %in% c(3,4,14),"north","other"))
test = wilcoxauc(subset(seu_filt, subset = mira_clusts %in% c(1,5,9,3,4,14)), 
                 'temp_subset', assay = "data", seurat_assay = "peaks_perclust") %>% 
  as.data.table
these = test[auc > quantile(test$auc, probs = 0.999)] # get top 0.1%
genes = ClosestFeature(seu_filt, these$feature)
these$closest_feature = genes$gene_name
VlnPlot(subset(seu_filt,microglia==TRUE), these$feature, assay = "peaks_perclust") + NoLegend()

library(ComplexHeatmap)
library(seriation)
mat = subset(seu_filt, subset = microglia == TRUE)@meta.data[,which(grepl("topic",colnames(seu_filt@meta.data),fixed=TRUE))] %>% as.data.table
mat$clust_ident = subset(seu_filt, subset = microglia == TRUE)$clust_ident
mat_avg = mat[, lapply(.SD,mean), by = clust_ident]
rows = mat_avg$clust_ident
mat_avg = mat_avg[, -c(1)] %>% as.matrix
rownames(mat_avg) = rows
o1 = seriate(dist(mat_avg), method = "TSP")
o2 = seriate(dist(t(mat_avg)), method = "TSP")
Heatmap(mat_avg, row_order = get_order(o1), column_order = get_order(o2))
VlnPlot(seu_filt, "ATAC_topic_2") + NoLegend() #enriched in Runx factors, specific to TIMs and Siglech-hi
VlnPlot(seu_filt, "ATAC_topic_3") + NoLegend() #enriched in C/EBPs and AP-1 proteins, specific to DAM-2

saveRDS(seu_filt, "multiome.rds")

# Further analysis post-MIRA ----
library(Seurat)
library(Signac)
library(data.table)
library(dplyr)
library(ggplot2)
library(presto)

seu_filt = readRDS("multiome.rds")
# Prepare microglia-only structure
DefaultAssay(seu_filt) = "RNA"
seu_filt = SCTransform(seu_filt)
seu_filt = RunPCA(seu_filt)
DefaultAssay(seu_filt) = "peaks_perclust"
seu_filt = FindTopFeatures(seu_filt, min.cutoff = 5)
seu_filt = RunTFIDF(seu_filt)
seu_filt = RunSVD(seu_filt)
# make integrated UMAP
seu_filt = FindMultiModalNeighbors(object = seu_filt,
                                   reduction.list = list("pca", "lsi"), 
                                   dims.list = list(1:50, 2:40),
                                   modality.weight.name = "RNA.weight",
                                   verbose = TRUE)
seu_filt = RunUMAP(object = seu_filt,
                   nn.name = "weighted.nn",
                   assay = "RNA",
                   verbose = TRUE, return.model = TRUE)
DefaultAssay(seu_filt) = "RNA"
seu_filt = FindNeighbors(seu_filt, dims = 1:30)
seu_filt = FindClusters(seu_filt, resolution = 1.2, algorithm = 3, graph.name = "wsnn")
clust.names = c("Homeostatic Microglia", "Homeostatic Microglia", "Selplg-lo Microglia","DAM-2","Homeostatic Microglia",
                "mt-Enriched Microglia", "TIMs","Siglech-hi Microglia","mt-Depleted Microglia","DAM-1")
names(clust.names) = levels(seu_filt$seurat_clusters)
seu_filt = RenameIdents(seu_filt, clust.names)
seu_filt$mglia_ident = Idents(seu_filt)
seu_filt$mglia_ident = factor(seu_filt$mglia_ident, levels = c("Homeostatic Microglia","Selplg-lo Microglia","mt-Enriched Microglia","mt-Depleted Microglia",
                                                               "DAM-1","DAM-2","TIMs","Siglech-hi Microglia"))
Idents(seu_filt) = seu_filt$mglia_ident

# I re-map fragment path due to moving files between computers for this analysis
Fragments(seu_filt@assays$peaks_perclust) = NULL
Fragments(seu_filt@assays$peaks) = NULL
frags = CreateFragmentObject(path = "atac_fragments.tsv.gz", cells = colnames(seu_filt), validate.fragments = TRUE)
Fragments(seu_filt@assays$peaks_perclust) = frags
Fragments(seu_filt@assays$peaks) = frags

# transfer motif scoring from MIRA to our seurat object
mira_motifs = fread("MIRA_motif_scores.csv") # generated by exporting motif_scores in MIRA
mira_motifs_dat = mira_motifs[,c(-1)] %>% as.matrix %>% t
colnames(mira_motifs_dat) = mira_motifs$V1
mira_assay = CreateAssayObject(mira_motifs_dat)
seu_filt[["MIRAmotifs"]] = mira_assay
# plot our scores in UMAP space
FeaturePlot(seu_filt, reduction = "umap",
            features = colnames(seu_filt@meta.data)[57:67], ncol = 4) %>% ggsave("Rplots.pdf",.,width = 20, height = 10, units = "in")
FeaturePlot(seu_filt, reduction = "umap",
            features = colnames(seu_filt@meta.data)[68:87], ncol = 5) %>% ggsave("Rplots.pdf",.,width = 20, height = 15, units = "in")
# plot some of the top hits from topics 2, 3, 16
(VlnPlot(seu_filt, c("RUNX1","RUNX2","RUNX3"), assay = "MIRAmotifs", 
         group.by = "subclust_id", idents = levels(seu_filt$subclust_id)[1:8]) + NoLegend()) %>% ggsave("Rplots.pdf",.,height=10,width=15,units="in")
(VlnPlot(seu_filt, c("USF2","USF1","CEBPE","CEBPB"), assay = "MIRAmotifs", ncol = 4,
         group.by = "subclust_id", idents = levels(seu_filt$subclust_id)[1:8]) + NoLegend()) %>% ggsave("Rplots.pdf",.,height=10,width=20,units="in")
(VlnPlot(seu_filt, c("FOSL1","JUND","JUN.VAR.2"), assay = "MIRAmotifs", ncol = 3,
         group.by = "subclust_id", idents = levels(seu_filt$subclust_id)[1:8]) + NoLegend()) %>% ggsave("Rplots.pdf",.,height=10,width=15,units="in")
# calculate motifs enriched in TIMs manually
tims_mira = RunPresto(seu_filt, ident.1 = "TIMs", ident.2 = levels(seu_filt$subclust_id)[c(1:5,7:8)],
                      min.pct = 0.1, logfc.threshold = 0.5, assay = "MIRAmotifs")
# plot some interesting ones (skip #1 bc it's a worm motif)
(VlnPlot(seu_filt, rownames(tims_mira)[2:11], assay = "MIRAmotifs", ncol = 5,
         group.by = "subclust_id", idents = levels(seu_filt$subclust_id)[1:8]) + NoLegend()) %>% ggsave("Rplots.pdf",.,height=20,width=25,units="in")
# prep some QC plots
VlnPlot(
  object = seu_filt,
  features = c('nCount_RNA','nFeature_RNA','TSS.enrichment','nucleosome_signal'),
  pt.size = 0.1,ncol = 4, group.by = "orig.ident") %>% ggsave("Rplots.pdf",.,height=10,width=30,units="in")
VlnPlot(
  object = seu_filt,
  features = c('nCount_peaks_perclust','nFeature_peaks_perclust'),
  pt.size = 0.1,ncol = 2, group.by = "subclust_id") %>% ggsave("Rplots.pdf",.,height=10,width=20,units="in")
TSSPlot(seu_filt, assay = "ATAC") %>% ggsave("Rplots.pdf",.,scale=1.5)
# plot TIM markers so we know TIMs are real
(VlnPlot(seu_filt, c("Fos","Jun","Junb"), assay = "SCT", 
         group.by = "subclust_id", idents = levels(seu_filt$subclust_id)[1:8]) + NoLegend()) %>% 
  ggsave("Rplots.pdf",.,height = 8, width = 12, units = "in")
# use presto to identify diff expressed peaks
seu_filt$temp_subset = ifelse(seu_filt$subclust_id == "TIMs", "TIMs", ifelse(seu_filt$microglia == TRUE,"Non-TIMs","Non-Microglia"))
test = wilcoxauc(subset(seu_filt, subset = microglia == TRUE), 'temp_subset', assay = "data", seurat_assay = "peaks_perclust") %>% as.data.table
these = test[padj < 1e-10 & abs(logFC) > 0.10 & group == "TIMs"]
genes = ClosestFeature(seu_filt, these$feature)
these$closest_feature = genes$gene_name
# plot those diff expressed peaks
(VlnPlot(subset(seu_filt,microglia==TRUE), these$feature, assay = "peaks_perclust") + NoLegend()) %>% ggsave("Rplots.pdf",.,scale=4)
(CoveragePlot(object = subset(seu_filt, subset = microglia == TRUE),
              region = these[5,feature],
              features = these[5,closest_feature],
              expression.assay = "SCT",
              annotation = TRUE,
              peaks = TRUE,
              tile = FALSE,
              links = TRUE,
              extend.upstream = 10000,
              extend.downstream = 10000)) %>% ggsave("Rplots.pdf",.,width=10,height=10,units="in")
# let's choose which genes might be juicy to peek at via presto again
tims_SCT = RunPresto(seu_filt, ident.1 = "TIMs", ident.2 = levels(seu_filt$subclust_id)[c(1:5,7:8)],
                     min.pct = 0.25, logfc.threshold = 0.5, assay = "SCT")
highlight = StringToGRanges("chr18-34859500-34871000")
highlight$color = "burlywood3"
# General summary plots ----
## Fig. S3C ----
p1 = CoveragePlot(
  object = seu_filt, idents = levels(seu_filt)[1:8], expression.assay = "SCT",
  features = "Egr1", region = "Egr1", annotation = TRUE,
  peaks = TRUE, links = TRUE, tile = FALSE,
  extend.upstream = 10000, extend.downstream = 10000,
  region.highlight = highlight)
highlight = StringToGRanges("chr1-134070500-134079900")
highlight$color = "burlywood3"
p2 = CoveragePlot(
  object = seu_filt, idents = levels(seu_filt)[1:8], expression.assay = "SCT",
  features = "Btg2", region = "Btg2", annotation = TRUE,
  peaks = TRUE, links = TRUE, tile = FALSE,
  extend.upstream = 5000, extend.downstream = 5000,
  region.highlight = highlight)
(p1|p2) %>% ggsave("tim_markers_covplots.png",.,width=15,height=6,dpi=600)

# make some more plots showing off concordances
## Fig. S3D ----
highlight = StringToGRanges("chr7-19698500-19699500")
highlight$color = "burlywood3"
p1 = CoveragePlot(
  object = seu_filt,idents = levels(seu_filt)[1:8],expression.assay = "SCT",
  features = "Apoe", region = "Apoe",annotation = TRUE,
  peaks = TRUE,links = TRUE,tile = FALSE,
  extend.upstream = 5000,extend.downstream = 5000,
  region.highlight = highlight)
highlight = StringToGRanges("chr11-96452000-96470000")
highlight$color = "burlywood3"
p2 = CoveragePlot(
  object = seu_filt,idents = levels(seu_filt)[9:23],expression.assay = "SCT",
  features = "Skap1", region = "Skap1",annotation = TRUE,
  peaks = TRUE,links = TRUE,tile = FALSE,
  extend.upstream = 15000,extend.downstream = 15000,
  region.highlight = highlight)
highlight = StringToGRanges("chr9-110419200-110429000")
highlight$color = "burlywood3"
p3 = CoveragePlot(
  object = seu_filt,idents = levels(seu_filt)[9:23],expression.assay = "SCT",
  features = "Ngp", region = "Ngp",annotation = TRUE,
  peaks = TRUE,links = TRUE,tile = FALSE,
  extend.upstream = 5000,extend.downstream = 5000,
  region.highlight = highlight)
highlight = StringToGRanges("chr11-44610000-44625200")
highlight$color = "burlywood3"
p4 = CoveragePlot(
  object = seu_filt,idents = levels(seu_filt)[9:23],expression.assay = "SCT",
  features = "Ebf1", region = "Ebf1",annotation = TRUE,
  peaks = TRUE,links = TRUE,tile = FALSE,
  extend.upstream = 15000,extend.downstream = 15000,
  region.highlight = highlight)
wrap_plots(list(p1,p2,p3,p4)) %>% ggsave("other_markers_covplots.png",.,width=15,height=12,dpi=600)

# make zoomed insets
highlight = StringToGRanges("chr11-96452000-96470000")
highlight$color = "burlywood3"
CoveragePlot(
  object = seu_filt,idents = levels(seu_filt)[9:23],expression.assay = "SCT",
  region = highlight, links = FALSE, peaks = FALSE, annotation = FALSE,
  region.highlight = highlight) %>% ggsave("other_markers_covplots_zoomed_skap1.png",.,width=5,height=3,dpi=600)
highlight = StringToGRanges("chr11-44610000-44625200")
highlight$color = "burlywood3"
CoveragePlot(
  object = seu_filt,idents = levels(seu_filt)[9:23],expression.assay = "SCT",
  region = highlight, links = FALSE, peaks = FALSE, annotation = FALSE,
  region.highlight = highlight) %>% 
  ggsave("other_markers_covplots_zoomed_ebf1.png",.,width=5,height=3,dpi=600)


# compare frequency of TIMs across our atlases!
## Fig. S3B ----
integrated_mglia = readRDS("mglia_only.rds")
multiome_mglia = readRDS("multiome_mglia.rds")

int.md = integrated_mglia@meta.data %>% as.data.table
mult.md = multiome_mglia@meta.data %>% as.data.table

int.md$is_tim = ifelse(int.md$mglia_ident %in% c("TIMs","Serpine1+ TIMs","Effector-hi TIMs"), TRUE, FALSE)
setkey(int.md, is_tim, orig.ident)
freqs = int.md[CJ(is_tim, orig.ident, unique = TRUE),.N,by=.EACHI] %>% 
  dcast(orig.ident ~ is_tim, value.var = "N") %>% 
  mutate(frac = 100*`TRUE`/(`TRUE` + `FALSE`))
mult.md$is_tim = ifelse(mult.md$mglia_ident == "TIMs", TRUE, FALSE)
mult.md$orig.ident = "E4_1yr_Multiome"
setkey(mult.md, is_tim, orig.ident)
freqs_2 = mult.md[CJ(is_tim, orig.ident, unique = TRUE),.N,by=.EACHI] %>% 
  dcast(orig.ident ~ is_tim, value.var = "N") %>% 
  mutate(frac = 100*`TRUE`/(`TRUE` + `FALSE`))
freqs = rbind(freqs,freqs_2) %>% 
  mutate(orig.ident = factor(orig.ident, levels = c("E2_10wk","E3_10wk","E4_10wk",
                                                    "E2_20wk","E3_20wk","E4_20wk",
                                                    "E4_1yr_Multiome","E3_2yr","E4_2yr")))
library(ggforce)
(ggplot(freqs,aes(x=orig.ident,y=frac,fill=orig.ident)) + 
    geom_bar(stat = "identity", color = "black") + 
    labs(x = "Sample", y = "Frequency of TIMs (%)") + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    theme(plot.margin = margin(5,5,5,20)) + 
    facet_zoom(y = orig.ident %in% c("E2_10wk","E3_10wk","E4_10wk","E2_20wk","E3_20wk","E4_20wk"), zoom.size = 0.5) +
    Seurat::NoLegend()) %>% 
  ggsave("tim_freq_per_sample.png",.,height=7,width=7,dpi=600)  

# DE analysis ----
# now let's run signac's built-in motif scanning functionality
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)

# get motif position frequency matrix from JASPAR and add to object
pfm = getMatrixSet(x = JASPAR2020, opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
seu_filt = AddMotifs(object = seu_filt, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pfm) #takes 5min or so

# compare TIMs and all other microglia
da_peaks = FindMarkers(seu_filt, ident.1 = "TIMs", ident.2 = levels(seu_filt$subclust_id)[c(1:5,7:8)],
                       only.pos = TRUE, test.use = "LR", min.pct = 0.05, 
                       latent.vars = "nCount_peaks_perclust", assay = "peaks_perclust") # slow... ~12min
top.da.peak = rownames(da_peaks[da_peaks$p_val < 0.005, ])
# get background
open.peaks = AccessiblePeaks(seu_filt, idents = levels(seu_filt$subclust_id)[1:8])
meta.feature = GetAssayData(seu_filt, assay = "peaks_perclust", slot = "meta.features")
peaks.matched = MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[top.da.peak, ],
  n = 50000)
# find enriched motifs
enriched.motifs = FindMotifs(object = seu_filt, features = top.da.peak, background = peaks.matched)
MotifPlot(object = seu_filt,motifs = head(rownames(enriched.motifs),n=10)) %>% ggsave("Rplots.pdf",.)

# run chromVAR
seu_filt = RunChromVAR(object = seu_filt, genome = BSgenome.Mmusculus.UCSC.mm10)
DefaultAssay(seu_filt) = "chromvar"
FeaturePlot(seu_filt, features = enriched.motifs$motif[1], min.cutoff = 'q10', max.cutoff = 'q90', pt.size = 0.1) %>% 
  ggsave("Rplots.pdf",.)

saveRDS(seu_filt, "multiome.rds")

motif_names = GetMotifData(seu_filt, assay = "peaks_perclust", slot = "motif.names") %>% unlist
chromvar.diff = FindMarkers(seu_filt, ident.1 = "TIMs", ident.2 = c("DAM-2","Homeostatic Microglia"),
                            only.pos = TRUE, mean.fxn = rowMeans, fc.name = "avg_diff") # almost instant! :tada:
chromvar.diff$motif_names = motif_names[rownames(chromvar.diff)] %>% unname
chromvar.diff %>% head(n=16)

library(patchwork)
chromvar_plotter = function(i) {(FeaturePlot(seu_filt, features = rownames(chromvar.diff)[i], min.cutoff = 'q10', 
                                             max.cutoff = 'q90', pt.size = 0.1) + ggtitle(chromvar.diff$motif_names[i]))}
wrap_plots(lapply(1:16, chromvar_plotter)) %>% ggsave("Rplots.pdf",.,width=20,height=20,units="in")

# Further summary plots ----
## Fig. S3A ----
markers = FindAllMarkers(seu_filt, assay = "SCT", min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC) %>% print(n=Inf)

genes_toplot = c("Tanc2","P2ry12","Selplg","mt-Nd3","mt-Nd1","Ctsd","Apoe","Egr1","Plek","Siglech","Alcam","F13a1","Mrc1",
                 "Ngp","S100a9","Ebf1","Bank1","Ighm","Skap1","St6galnac3","Atp1a2")
(DotPlot(seu_filt,features=genes_toplot,assay="SCT")+ theme(axis.text.x = element_text(angle = 45, hjust=1))) %>% ggsave("markers_dotplot.pdf",.,height=8,width=14,dpi=600)

## Fig. 3A ----
(DimPlot(seu_filt,group.by = "subclust_id", label = TRUE, repel = TRUE, label.size = 8, pt.size = 2) + NoLegend() + 
    theme(plot.title = element_blank(),
          axis.text = element_text(size=20),
          axis.title = element_text(size=30),
          legend.text = element_text(size=22.5))) %>% 
  ggsave("all_umap.png",.,width=15.975,height=15.975,dpi=600)

# footprinting
## Fig. S3E ----
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)
motif_names = GetMotifData(seu_filt, assay = "peaks_perclust", slot = "motif.names") %>% unlist
chromvar.diff = FindMarkers(seu_filt, ident.1 = "DAM-2", ident.2 = c("TIMs","Homeostatic Microglia"), assay = "chromvar",
                            only.pos = TRUE, mean.fxn = rowMeans, fc.name = "avg_diff") # almost instant! :tada:
chromvar.diff$motif_names = motif_names[rownames(chromvar.diff)] %>% unname
chromvar.diff %>% head(n=16)

# for TIMs: NFKB2
# for DAM-2: CEBPD
seu_filt = Footprint(seu_filt, motif.name = c("NFKB2","CEBPD"), genome = BSgenome.Mmusculus.UCSC.mm10)
temp = Idents(seu_filt) %>% as.character
temp[temp == "Homeostatic Microglia"] = "Hom"
Idents(seu_filt) = factor(temp, levels = c("Hom",levels(seu_filt)[-1]))
fprints = PlotFootprint(seu_filt, features = c("CEBPD","NFKB2"), idents = c("Hom","DAM-2","TIMs"), label = FALSE) &
  theme(plot.title = element_text(size=30), axis.title = element_text(size=25),axis.text = element_text(size=20), 
        legend.title = element_blank(), legend.text = element_text(size=20)) &
  guides(color = guide_legend(override.aes = list(linewidth=5)))
(fprints + patchwork::plot_layout(ncol=2)) %>% ggsave("footprints.png",.,width=17.5,height=6,dpi=800)
Idents(seu_filt) = seu_filt$subclust_id

## Fig. 3B ----
library(SeuratWrappers)
library(ggrepel)
res = RunPresto(seu_filt, ident.1 = "TIMs", ident.2 = "DAM-2", group.by = "subclust_id", logfc.threshold = 0, min.pct = 0, assay = "chromvar") %>% as.data.table(keep.rownames = TRUE)
motif_names = GetMotifData(seu_filt, assay = "peaks_perclust", slot = "motif.names") %>% unlist
res$motif = motif_names[res$rn]
res$logp = -log10(res$p_val_adj)
fc_cutoff = 0.75
p_cutoff = -log(0.05)
res[, group := ifelse(abs(avg_log2FC) >= fc_cutoff & logp >= p_cutoff, "Hit", "Filter")]
res$group = factor(res$group, levels = c("Hit","Filter"))
res$label = FALSE
res[, in_zone := ifelse(abs(avg_log2FC) >= fc_cutoff & logp >= p_cutoff, TRUE, FALSE)]
res[avg_log2FC > 0 & in_zone, label := TRUE]
res[avg_log2FC < 0 & in_zone, label := TRUE]

options(ggrepel.max.overlaps = Inf)
volcano = ggplot(res, aes(x=avg_log2FC,y=logp)) + 
  geom_point(aes(color=res$group), size=5) +
  scale_color_manual(values = c("blue","gray")) +
  geom_text_repel(label = ifelse(res$label==TRUE,res$motif,''),
                  size = 10, force = 10, box.padding = 0.2, point.padding = 0.2, segment.size = 0.1, min.segment.length = 0.1, max.iter = 1e7) +
  geom_vline(xintercept = -fc_cutoff) + 
  geom_vline(xintercept = fc_cutoff) + 
  geom_hline(yintercept = p_cutoff) +
  geom_text(data = data.frame(xpos = c(-(fc_cutoff+0.2), fc_cutoff+0.2), 
                              ypos = c(-5, -5), hjustvar = c(1,0), vjustvar = c(-1,-1), 
                              annotateText = c(paste0("DAM-2 ",sprintf('\U2190')),paste0(sprintf('\U2192')," TIMs"))),
            aes(x=xpos, y=ypos, hjust=hjustvar, vjust=vjustvar, label=annotateText), size = 15) +
  theme_bw() + 
  theme(plot.title = element_blank(),
        axis.text = element_text(size=20),
        axis.title = element_text(size=30)) +
  labs(x = "Log-2 Fold Change", y = "Log-10 P-Value") + 
  NoLegend()

ggsave("chromvar_volcano_tim_dam.png",volcano,dpi=600,width=15.975,height=15.975,units="in")

# make nice mira topics plot
## Fig. 3C ----
library(patchwork)
DefaultAssay(seu_filt) = "MIRAmotifs"
Idents(seu_filt) = seu_filt$subclust_id
dat = seu_filt@meta.data %>% as.data.table
dat = dat[subclust_id %in% levels(seu_filt)[1:8]]
topic2 = fread("~/MIRA/atac_topic_2_tfs.csv") # output from MIRA
topic3 = fread("~/MIRA/atac_topic_3_tfs.csv") # output from MIRA
topic15 = fread("~/MIRA/atac_topic_15_tfs.csv") # output from MIRA

miraplotter = function(x,tb){
  vln = ggplot(dat, aes(x = subclust_id, y = !!sym(paste0("ATAC_topic_",x)), fill = subclust_id)) +
    geom_boxplot() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust=1, size = 22.5),
          axis.text.y = element_text(size = 20),
          axis.title = element_text(size=25),
          plot.title = element_text(size=25)) + 
    labs(x = "", y = "", title = paste0("MIRA Topic ",x)) + 
    NoLegend()
  tb$logp = -log10(tb$pval)
  tb = head(tb,n=15) %>% 
    mutate(parsed_name = forcats::fct_reorder(parsed_name,logp))
  topic = ggplot(tb, aes(x=logp, y = parsed_name)) +
    geom_segment(aes(y=parsed_name,yend=parsed_name,x=0,xend=logp), lwd=1) + 
    geom_point(size = 4) +
    geom_vline(xintercept=-log10(0.05), linetype="dotted", lwd = 1, color = "red") +
    labs(x = "Log-10 P-Value", y = "") + 
    theme_bw() + 
    theme(plot.title = element_blank(),
          axis.text.y = element_text(size=22.5),
          axis.text.x = element_text(size=20),
          axis.title = element_text(size=25))
  return(vln/topic)
}

wrap_plots(list(miraplotter(2,topic2),miraplotter(3,topic3),miraplotter(15,topic15)),ncol=3) %>% 
  ggsave("MIRA_plot.png",.,width=15.975,height=15.975,dpi=600)

# SCENIC+ plots ----
## Fig. 3D ----
eregulon_df = fread("eregulons_df.csv") #output from SCENIC+, see accompanying code
eregulon_df$index = factor(eregulon_df$index, levels = levels(mglia$mglia_ident))
regulon_order = fread("regulon_order.csv",skip=1)$V2
eregulon_df$eRegulon_name = factor(eregulon_df$eRegulon_name, levels = regulon_order)

(ggplot(eregulon_df, aes(index, eRegulon_name)) +
    geom_tile(aes(fill = color_val)) +
    scale_fill_distiller(type = "div", palette = "RdYlBu") + 
    geom_point(aes(size = size_val), color = "black") +
    scale_size(range = c(3,10)) +
    theme_bw() +
    theme(axis.line=element_blank(),panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(),axis.text.x = element_text(angle = 45, hjust=1, size = 30),
          axis.text.y = element_text(size=25), legend.title = element_text(size=30), legend.text = element_text(size=25),
          axis.title.y = element_text(size=30), legend.spacing.y = unit(0.1, 'in')) +
    guides(size = guide_legend(byrow = TRUE), fill = guide_colorbar(barheight = 12)) +
    labs(x = "",y="eRegulon",fill="Expression",size="eRegulon Enrichment")) %>% 
  ggsave("scenicplus_eregulonplot.pdf",.,width=18.56,height=19,dpi=600)

## Fig. 3E ----
# made entirely in SCENIC+, see accompanying code

# CellPhoneDB ----
library(Seurat)
library(Signac)
library(dplyr)
library(data.table)
library(org.Hs.eg.db)

## Prep files ----
adapoe = readRDS("multiome.rds")
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
rm.these = which(is.na(rownames(counts)))
counts = counts[-rm.these,]
counts_dt = as.data.table(counts, keep.rownames = TRUE)
fwrite(counts_dt, "cpdb_counts_matrix.txt", sep = "\t", 
       row.names = FALSE, col.names = TRUE, quote = FALSE)

# also prep metadata file
Idents(adapoe) = adapoe$subclust_id
md = Idents(adapoe) %>% as.data.table(keep.rownames = TRUE)
colnames(md) = c("Cell","cell_type")
fwrite(md, "cpdb_metadata.txt", sep = "\t", 
       row.names = FALSE, col.names = TRUE, quote = FALSE)

## Plot after run ----
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
library(Cairo)
## get heatmap order from seurat idents
apoe = readRDS("multiome.rds")
cluster_order = levels(apoe$subclust_id)

## rename clusters and pairs
# outputs from CellPhoneDB
pvals = fread("pvalues.txt") 
means = fread("means.txt")

setnames(pvals, colnames(pvals[,12:540]), 
         sapply(colnames(pvals[,12:540]), function(x) str_replace(x, fixed("|"), paste0(" ",sprintf("\U2192")," "))) %>% as.character) #these are arrows
setnames(means, colnames(means[,12:540]), 
         sapply(colnames(means[,12:540]), function(x) str_replace(x, fixed("|"), paste0(" ",sprintf("\U2192")," "))) %>% as.character) # these are arrows
pvals$interacting_pair = sapply(pvals$interacting_pair, function(x) str_replace(x, fixed("_"), paste0(" ",sprintf("\U00D7")," "))) # these are x signs
means$interacting_pair = sapply(means$interacting_pair, function(x) str_replace(x, fixed("_"), paste0(" ",sprintf("\U00D7")," "))) # these are x signs

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
              pvalues_file = "out/pvalues.txt",
              count_filename = "heatmap_count.png",
              log_filename = "heatmap_log_count.png",
              count_network_filename = "out/network.txt",
              interaction_count_filename = "out/interaction_count.txt",
              count_network_separator = "\t",
              interaction_count_separator = "\t")

## Fig. 3F ----
dat = fread("network.txt") # generated in heatmaps_plot
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

# struct for subclust_id ordering
group = structure(c(rep("HMs",4),rep("IMs",2),rep("TIMs",1),rep("HMs",1),rep("Other Immune",15)),
                  names = rownames(dat.m))
grid.col = structure(c(rep(2,4),rep(3,2),rep(4,1),rep(2,1),rep(6,15)),
                     names = rownames(dat.m))
circos.clear()
chordDiagram(dat.m_filt, group = group, grid.col = grid.col, directional = 1,
             direction.type = c("arrows","diffHeight"), link.arr.type = "big.arrow",
             diffHeight = mm_h(3), target.prop.height = mm_h(2))

#with rotated labels
png("circos.png", height = 9, width = 9, units = "in", res = 600)
par(bg=NA)
circos.par(circle.margin = c(0.25,0.2,0.05,0.05))
chordDiagram(dat.m_filt, group = group, grid.col = grid.col, directional = 1,
             direction.type = c("arrows","diffHeight"), link.arr.type = "big.arrow",
             diffHeight = mm_h(3), target.prop.height = mm_h(2), 
             annotationTrack = "grid", preAllocateTracks = 1.1)
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)
dev.off()

## Fig. S3F ----
dat = fread("deconvoluted.txt") #output from cellphonedb
dat = subset(dat, !duplicated(dat[,complex_name]))
dat[complex_name == "", complex_name := protein_name]
dat = melt.data.table(dat, measure.vars = colnames(dat)[7:(length(cluster_order)+6)], variable.name = "cluster", value.name = "deconv") %>%
  setorder(-deconv)
dat$cluster = factor(dat$cluster, levels = levels(apoe$subclust_id))
dat$rank = seq_len(nrow(dat))
dat[, stat := (deconv - mean(deconv))/mean(deconv), by = cluster]
dat.mglia = dat[cluster %in% cluster_order[1:8]]
highlight.these = c("Dehydroepiandrosterone_bySTS","2arachidonoylglycerol_byDAGLB","TGFbeta_receptor1",
                    "TREM2_receptor","Histamine_byHDC","TGFbeta_receptor2")
dat_filt = dat.mglia[complex_name %in% dat.mglia[stat>=2, complex_name]]
dat_filt$complex_name = paste("<span style = 'color: ",
                              ifelse(dat_filt$complex_name %in% highlight.these,"blue","black"),
                              ";'>",
                              ifelse(dat_filt$complex_name %in% highlight.these, "<b>",""),
                              dat_filt$complex_name,
                              ifelse(dat_filt$complex_name %in% highlight.these, "</b>",""),
                              "</span>",sep = "")
(ggplot(dat_filt, aes(y = reorder(complex_name,-stat), x = stat)) +
    geom_boxplot(size=1.5) + 
    labs(y = "", x = "Mean Fold-Change\nover Cluster Mean") +
    coord_flip() +
    theme_bw() +
    theme(axis.text.x = ggtext::element_markdown(size=30, angle = 45, hjust = 1), axis.text.y = element_text(size=25),
          axis.title = element_text(size=35), plot.margin = margin(5,5,5,200))) %>% 
  ggsave("mglia_top_complexes.png",.,width=45.5,height=12,dpi=600)

#comparing TIMs to hom or DAM complexes
## Fig. 3G ----
dat.mglia[cluster %in% cluster_order[1:4], type := "Homeostatic"]
dat.mglia[cluster %in% cluster_order[7], type := "TIMs"]
dat.mglia[cluster %in% cluster_order[5:6], type := "DAMs"]

comp = dat.mglia[, mean(stat), by = .(complex_name,type)]
comp = comp[!is.na(type)]
comp$V1 = log(1.25+comp$V1) #log-norm to get rid of infs/nas
comp = dcast(comp, complex_name ~ type, value.var = "V1")
# for hom:
comp$score = comp$TIMs - comp$Homeostatic
comp$label = ifelse(comp$score > 0.4 | comp$score < -0.4, comp$complex_name, "")
(ggplot(comp, aes(x = reorder(complex_name,score), y = score, label = label)) + 
    geom_point() + 
    geom_text_repel(max.overlaps = Inf, box.padding = 0.5, force = 20) +
    labs(x = "Complex", y = "TIM Score - Homeostatic Score") +
    theme_bw() + 
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())) %>% 
  ggsave("tim_vs_hom_complexes.png",.,width=10,height=7,dpi=600)

int_comp = fread("tim_hom_comp.csv") # generated when making Figure 2
comp$concordance = "No"
comp[score > 0, concordance := ifelse(complex_name %in% int_comp[score>0,complex_name], "Yes", "No")]
comp[score < 0, concordance := ifelse(complex_name %in% int_comp[score<0,complex_name], "Yes", "No")]
comp$concordance = factor(comp$concordance, levels = c("Yes", "No"))

comp %>% 
  mutate(facet = ifelse(score >= 0, "Up in TIMs", "Up in Homeostatic")) %>% 
  mutate(facet = factor(facet, levels = c("Up in TIMs", "Up in Homeostatic"))) %>%
  arrange(desc(score)) %>% 
  filter((row_number() > max(row_number()) - 15 & score<0) | (row_number() <= 15 & score>0)) %>%
  mutate(complex_name = factor(complex_name, levels = c(rev(complex_name[1:15]),complex_name[16:30]))) %>%
  {ggplot(.,aes(x = score, y = complex_name)) + 
      geom_segment(aes(y=complex_name,yend=complex_name,x=0,xend=score), lwd=1.25) +
      geom_point(aes(color = concordance)) + 
      scale_color_manual(values = c("Yes" = "red", "No" = "blue")) + 
      labs(x = "TIM Score - Homeostatic Score", y = "", color = "Concordant with\nscRNAseq Dataset?") +
      theme_bw() + 
      facet_wrap(~facet,scales="free",ncol=1)} %>% 
  ggsave("tim_hom_complexes_lollipop.png",.,width=6,height=8,dpi=600)

# for dam:
comp$score = comp$TIMs - comp$DAMs
comp$label = ifelse(comp$score > 0.4 | comp$score < -0.4, comp$complex_name, "")
(ggplot(comp, aes(x = reorder(complex_name,score), y = score, label = label)) + 
    geom_point() + 
    geom_text_repel(max.overlaps = Inf, box.padding = 0.5, force = 20) +
    labs(x = "Complex", y = "TIM Score - DAM Score") +
    theme_bw() + 
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())) %>% 
  ggsave("tim_vs_dam_complexes.png",.,width=10,height=7,dpi=600)

int_comp = fread("tim_dam_comp.csv") # generated when making Figure 2
comp$concordance = "No"
comp[score > 0, concordance := ifelse(complex_name %in% int_comp[score>0,complex_name], "Yes", "No")]
comp[score < 0, concordance := ifelse(complex_name %in% int_comp[score<0,complex_name], "Yes", "No")]
comp$concordance = factor(comp$concordance, levels = c("Yes", "No"))

comp %>% 
  mutate(facet = ifelse(score >= 0, "Up in TIMs", "Up in DAMs")) %>% 
  mutate(facet = factor(facet, levels = c("Up in TIMs", "Up in DAMs"))) %>%
  arrange(desc(score)) %>% 
  filter((row_number() > max(row_number()) - 15 & score<0) | (row_number() <= 15 & score>0)) %>%
  mutate(complex_name = factor(complex_name, levels = c(rev(complex_name[1:15]),complex_name[16:30]))) %>%
  {ggplot(.,aes(x = score, y = complex_name)) + 
      geom_segment(aes(y=complex_name,yend=complex_name,x=0,xend=score), lwd=1.5) +
      geom_point(aes(color = concordance), size = 5) + 
      scale_color_manual(values = c("Yes" = "red", "No" = "blue")) + 
      labs(x = "TIM Score - DAM Score", y = "", color = "Concordant with\nscRNAseq Dataset?") +
      theme_bw() + 
      theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=25),
            axis.title = element_text(size=25), strip.text = element_text(size=27.5),
            legend.title = element_text(size=25), legend.text = element_text(size=22.5),
            legend.spacing.y = unit(0.2, 'in')) +
      guides(color=guide_legend(ncol=1,override.aes = list(size = 6),byrow = TRUE)) +
      facet_wrap(~facet,scales="free",ncol=1)} %>% 
  ggsave("tim_dam_complexes_lollipop.png",.,width=12.5,height=16.667,dpi=600)

# Compass ----
## Prep data ----
library(Seurat)
library(Signac)
library(dplyr)
library(data.table)

mglia = readRDS("multiome_mglia.rds")

md = mglia@meta.data %>% as.data.table
md$samp = 1:nrow(md)
md$agg_id = md$mglia_ident
md = split(md, md$agg_id)
for(i in 1:length(md)){
  num_to_pseudobulk = ceiling(nrow(md[[i]])/10) # pseudobulk into 10 cells per "cell"
  md[[i]]$pseudobulk_id = paste0(md[[i]]$agg_id,"_",rep(1:num_to_pseudobulk))
}
md = do.call("rbind",md)
setorder(md, "samp")
mglia$pseudobulk_id = md$pseudobulk_id
out = AverageExpression(mglia, group.by = "pseudobulk_id", assays = "RNA", slot = "counts")$RNA
write.table(out, file = "dat_matrix.tsv", sep = "\t", quote = FALSE)

type = unname(colnames(out))
type = sapply(strsplit(type, "_"), function(x) paste(x[1], collapse = "-"))
metadata = data.frame("type" = type)
rownames(metadata) = unname(colnames(out))
write.table(metadata, file = "cell_metadata.csv", sep = "\t", quote = FALSE)

## Analyze after run ----
library(data.table)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(matrixStats)

setwd("~/scratch/R_dir/1yr_multiome/Compass/Normed")

dat = fread("reactions.tsv")
penalties = as.matrix(dat[,-1])
rownames(penalties) = dat$V1
cell_metadata = fread("~/scratch/R_dir/1yr_multiome/Compass/cell_metadata.csv")
grps = cell_metadata$type
rxn_metadata = fread("~/scratch/R_dir/compass_reaction_metadata.csv")

reaction.consistencies = -log1p(penalties)
reaction.consistencies = reaction.consistencies[which(apply(reaction.consistencies, 1, max) - apply(reaction.consistencies, 1, min) > 1e-3),]
reaction.consistencies = reaction.consistencies - min(reaction.consistencies)
reaction.consistencies = scale(reaction.consistencies)

(rowRanges(reaction.consistencies)[,2] - rowRanges(reaction.consistencies)[,1]) %>% summary
filt = reaction.consistencies[which((rowRanges(reaction.consistencies)[,2] - rowRanges(reaction.consistencies)[,1]) > 0.4),]

filt = t(filt)
rownames(filt) = NULL

# rename reactions
rxns = strsplit(colnames(filt),"_(?!.*_)",perl = TRUE) %>% sapply("[[",1) %>% as.data.table
rxns$dir = strsplit(colnames(filt),"_(?!.*_)",perl = TRUE) %>% sapply("[[",2)
rxns$dir = ifelse(rxns$dir == "pos", "(+)","(-)")
rxns = merge(rxns, rxn_metadata, by.x = ".", by.y = "reaction_no_direction")
rxns[, new_name := paste0(stringr::str_trunc(reaction_name,43)," ",dir)]
colnames(filt) = rxns$new_name
# prep for annotations
col_fun = randomcoloR::distinctColorPalette(length(unique(rxns$subsystem)))
names(col_fun) = unique(rxns$subsystem)
confs = ifelse(rxns$confidence %in% c(0,4), "High", "Low")

ht = Heatmap(filt, name = "COMPASS\nReaction\nScores", cluster_rows = cluster_within_group(t(filt),grps), 
             left_annotation = rowAnnotation(Cluster=anno_block(gp = gpar(fill=randomcoloR::distinctColorPalette(19)))), row_split = 19,
             bottom_annotation = HeatmapAnnotation(Subsystem = rxns$subsystem, Confidence = confs, 
                                                   col = list(Subsystem = col_fun, Confidence = c("High" = "lightblue","Low" = "grey")), 
                                                   annotation_name_side = "left"),
             top_annotation = HeatmapAnnotation(Density = anno_boxplot(filt, height = unit(4, "cm")), 
                                                annotation_name_side = "left", annotation_name_rot = 0),
             row_title_rot = 0, column_names_max_height = unit(35,"cm"))

png("hmap.png",width=40,height=35,units="in",res=600)
draw(ht)
dev.off()

## Fig. 3H ----
library(ggblend)
library(ggtext)
library(ggrepel)
rc = t(reaction.consistencies) %>% as.data.table(keep.rownames = TRUE)
rc$rn = strsplit(rc$rn, "_") %>% sapply(.,"[[",1)
means = rc[rn %in% c("TIMs","DAM-2"), lapply(.SD,mean), by = rn][,-1] %>% as.matrix
mean_diffs = colDiffs(means) %>% as.vector # this is TIM - DAM, as we want
sds = (rc[rn %in% c("DAM-2","TIMs"), lapply(.SD, sd), .SDcols = colnames(rc)[-1]][1] %>% as.matrix)[1,]
cohends = mean_diffs/sds
cohends = data.table(rxn = names(cohends), cohend = cohends)
cohends$rxnrename = strsplit(cohends$rxn,"_(?!.*_)",perl = TRUE) %>% sapply("[[",1) %>% as.data.table
cohends$dir = strsplit(cohends$rxn,"_(?!.*_)",perl = TRUE) %>% sapply("[[",2)
cohends$dir = ifelse(cohends$dir == "pos", "(+)","(-)")
cohends = merge(cohends, rxn_metadata, by.x = "rxnrename", by.y = "reaction_no_direction")
cohends[, new_name := paste0(stringr::str_trunc(rxnrename,43)," ",dir)]
# pull in padj values with presto
presto_res = presto::wilcoxauc(X = t(rc[rn %in% c("TIMs","DAM-2"),-1]), y = rc[rn %in% c("TIMs","DAM-2"),rn]) %>% as.data.table
presto_res = presto_res[group == "TIMs"]
cohends$padj = presto_res$padj
# order and color
cohends[subsystem == "Arginine and Proline Metabolism", subsystem := "Arginine and proline metabolism"]
highlight.these.low = c("Glycolysis/gluconeogenesis","Pentose phosphate pathway","Galactose metabolism","Fructose and mannose metabolism")
highlight.these.high = c("Taurine and hypotaurine metabolism","Histidine metabolism","Arginine and proline metabolism","Phenylalanine metabolism",
                         "Glutamate metabolism","Lysine metabolism","Cysteine metabolism","Valine, leucine, and isoleucine metabolism")
cohends$subsystem = paste("<span style = 'color: ",
                          ifelse(cohends$subsystem %in% highlight.these.low,"blue",
                                 ifelse(cohends$subsystem %in% highlight.these.high, "red","black")),
                          ";'>",
                          ifelse(cohends$subsystem %in% c(highlight.these.low,highlight.these.high), "<b>",""),
                          cohends$subsystem,
                          ifelse(cohends$subsystem %in% c(highlight.these.low,highlight.these.high), "</b>",""),
                          "</span>",sep = "")
cohends$subsystem = factor(cohends$subsystem, levels = cohends[, median(cohend), by = subsystem][order(V1)]$subsystem)
cohends[, color := ifelse(cohend >= 0, "TIMs", "DAM-2")]
cohends$color = factor(cohends$color, levels = c("TIMs","DAM-2"))
cohends[, shading := padj < 0.1]
cohends$shading = factor(cohends$shading, levels = c(TRUE,FALSE))

these.subsystems = c(head(levels(cohends$subsystem),n=30),tail(levels(cohends$subsystem),n=30))
cohends[subsystem %in% these.subsystems[1:30], facet := "Down in TIMs"]
cohends[subsystem %in% these.subsystems[31:60], facet := "Up in TIMs"]
cohends$facet = factor(cohends$facet, levels = c("Up in TIMs","Down in TIMs"))

(ggplot(cohends[subsystem %in% these.subsystems], aes(x=cohend,y=subsystem)) + 
    geom_point(aes(fill=color, alpha = shading),size=5,shape=21,color="black",lwd=0.25) %>% blend("source") + 
    scale_fill_manual(values = c("TIMs" = "red", "DAM-2" = "blue")) + 
    scale_alpha_discrete(range=c(1,0.4)) +
    theme_bw() + 
    labs(x = "Cohen's d", y= "", fill = "Enriched in:", alpha = "Adjusted P-Value\n< 0.1?") +
    theme(strip.text = ggtext::element_markdown(size=25), axis.text.y = ggtext::element_markdown(size=22.5),
          axis.text.x = element_text(size=17.5),
          axis.title = element_text(size=25),
          legend.title = element_text(size=25),
          legend.text = element_text(size=22.5),
          legend.spacing.y = unit(0.2, 'in')) + 
    guides(fill=guide_legend(override.aes = list(size = 6),byrow=TRUE),
           alpha=guide_legend(override.aes = list(size = 6, fill = "black"),byrow=TRUE)) +
    facet_wrap(~facet, scales = "free")
) %>% 
  ggsave("cohen_plot_vs_dams.png",.,height=16.856,width=35.116,dpi=600)