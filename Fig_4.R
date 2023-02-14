# Prep ROSMAP_wholebrain ----

setwd("ROSMAP_WholeBrain")

library(Seurat)
library(data.table)
library(dplyr)
library(Matrix)
library(ggplot2)

counts = fread("ROSMAP_Brain.snRNAseq_counts_sparse_format_20201107.csv")
i = as.vector(counts$i)
j = as.vector(counts$j)
x = as.vector(counts$x)
A = sparseMatrix(i, j, x = x)
cells = read.csv("ROSMAP_Brain.snRNAseq_metadata_cells_20201107.csv", header=T, stringsAsFactors = F)
genes = read.csv("ROSMAP_Brain.snRNAseq_metadata_genes_20201107.csv", header=T, stringsAsFactors = F)
rownames(A)=genes$x
colnames(A)=cells$cell_name
# convert gene names to mouse homologs for mapping
homologs = fread("~/HMD_HumanPhenotype.rpt.txt") # from MGI
homologs$iterator = 1:nrow(homologs)
homologs = homologs[, .SD[which.min(iterator)], by = V1] # we do this to get only one instance of each human gene
keep.these = genes$x %in% homologs$V1 %>% which
my_conv = data.table(num = 1:length(keep.these), orig = genes$x[keep.these])
my_conv = merge.data.table(my_conv, homologs, by.x = "orig", by.y = "V1")
my_conv = setkey(my_conv, num)
# filter and rename
A = A[keep.these,]
rownames(A) = my_conv$V3
# create seurat object
cells_asmd = cells[,2:4]
rownames(cells_asmd) = cells$cell_name
rosmap = CreateSeuratObject(counts = A, meta.data = cells_asmd)
rosmap = subset(rosmap, subset = broad_class == "Micr")
# prep for query mapping
rosmap = SCTransform(rosmap)
rosmap = RunPCA(rosmap)
rosmap = RunUMAP(rosmap, reduction = "pca", dims = 1:30, verbose = FALSE)
(DimPlot(rosmap, group.by = "subtype", label = TRUE, repel = TRUE) + NoLegend()) %>% ggsave("Rplots.pdf",.)

library(Signac)
library(SeuratWrappers)
library(patchwork)
seu = readRDS("~/multiome_mglia.rds")
DefaultAssay(seu) = "SCT"
anchors = FindTransferAnchors(reference = seu, query = rosmap, dims = 1:30, reference.reduction = "pca", normalization.method = "SCT")
predictions = TransferData(anchorset = anchors, refdata = seu$mglia_ident, dims = 1:30)

rosmap = MapQuery(anchorset = anchors, reference = seu, query = rosmap, refdata = "mglia_ident", 
                  reference.reduction = "pca", reduction.model = "umap")
rosmap$predicted.id = factor(rosmap$predicted.id, levels = levels(seu$mglia_ident))
p1 = DimPlot(seu, reduction = "umap", group.by = "mglia_ident", label = TRUE, label.size = 3,
             repel = TRUE) + ggtitle("Reference annotations")
p2 = DimPlot(rosmap, reduction = "ref.umap", group.by = "predicted.id", label = TRUE,
             label.size = 3, repel = TRUE) + ggtitle("Query transferred labels")
(p1 + p2) %>% ggsave("Rplots.tiff",.,width=14,height=7)
saveRDS(rosmap, "rosmap_projected.rds")
# save our metadata separately for easy access
fwrite(as.data.table(rosmap@meta.data),"rosmap_wholebrain_md.csv")
fwrite(as.data.table(seu@meta.data),"~/multiome_md.csv")

# Prep BA10 ----
setwd("snRNAseqPFC_BA10")

library(Seurat)
library(data.table)
library(dplyr)
library(Matrix)
library(ggplot2)

cells = fread("filtered_column_metadata.txt")
genes = read.csv("filtered_gene_row_names.txt", header = F, stringsAsFactors = F)
A = readMM("filtered_count_matrix.mtx")
rownames(A)=genes$V1
colnames(A)=cells$TAG

# convert gene names to mouse homologs for mapping
homologs = fread("~/HMD_HumanPhenotype.rpt.txt") # from MGI
homologs$iterator = 1:nrow(homologs)
homologs = homologs[, .SD[which.min(iterator)], by = V1] # we do this to get only one instance of each human gene
keep.these = genes$V1 %in% homologs$V1 %>% which
my_conv = data.table(num = 1:length(keep.these), orig = genes$V1[keep.these])
my_conv = merge.data.table(my_conv, homologs, by.x = "orig", by.y = "V1")
my_conv = setkey(my_conv, num)
# filter and rename
A = A[keep.these,]
rownames(A) = my_conv$V3
# create seurat object
cells_asmd = cells[,-1]
rownames(cells_asmd) = cells$TAG
BA10 = CreateSeuratObject(counts = A, meta.data = cells_asmd)
BA10 = subset(BA10, subset = broad.cell.type == "Mic")

# prep for query mapping
BA10 = SCTransform(BA10)
BA10 = RunPCA(BA10)
BA10 = RunUMAP(BA10, reduction = "pca", dims = 1:30, verbose = FALSE)
(DimPlot(BA10, group.by = "Subcluster", label = TRUE, repel = TRUE) + NoLegend()) %>% ggsave("Rplots.pdf",.)

# now we get our own multiome data so we can project sn -> sn
library(Signac)
library(SeuratWrappers)
library(patchwork)
seu = readRDS("~/multiome_mglia.rds")
DefaultAssay(seu) = "SCT"
# transfer.
anchors = FindTransferAnchors(reference = seu, query = BA10, dims = 1:30, reference.reduction = "pca", normalization.method = "SCT")
predictions = TransferData(anchorset = anchors, refdata = seu$mglia_ident, dims = 1:30)

BA10 = MapQuery(anchorset = anchors, reference = seu, query = BA10, refdata = "mglia_ident", 
                reference.reduction = "pca", reduction.model = "umap")
BA10$predicted.id = factor(BA10$predicted.id, levels = levels(seu$mglia_ident))
p1 = DimPlot(seu, reduction = "umap", group.by = "mglia_ident", label = TRUE, label.size = 3,
             repel = TRUE) + ggtitle("Reference annotations")
p2 = DimPlot(BA10, reduction = "ref.umap", group.by = "predicted.id", label = TRUE,
             label.size = 3, repel = TRUE) + ggtitle("Query transferred labels")
(p1 + p2) %>% ggsave("Rplots.pdf",.,width=14,height=7)

BA10_md = BA10@meta.data %>% as.data.table
seu_md = seu@meta.data %>% as.data.table

freqs_BA10 = BA10_md[, .N/nrow(BA10_md), by = "predicted.id"]
freqs_seu = seu_md[, .N/nrow(seu_md), by = "mglia_ident"]
freqs = merge.data.table(x = freqs_BA10, y = freqs_seu, by.x = "predicted.id", by.y = "mglia_ident") %>% 
  melt.data.table(., id.vars = "predicted.id", measure.vars = c("V1.x", "V1.y"))

saveRDS(BA10, "BA10_projected.rds")
# save our metadata separately for easy access
fwrite(as.data.table(BA10@meta.data),"BA10_md.csv")

# Prep SEA-AD ----
setwd("SEA-AD")

library(Seurat)
library(data.table)
library(dplyr)
library(org.Hs.eg.db)
library(ggplot2)
library(rlang)

# microglia only from the Allen dataset:
# https://cellxgene.cziscience.com/collections/1ca90a2d-2943-483d-b678-b809bf464c30

sea = readRDS("allen_mglia.rds")

# extract counts matrix and overwrite with mouse homologs
counts = GetAssayData(sea, slot = "counts")
convs = select(org.Hs.eg.db, keys = rownames(counts), columns = "SYMBOL", keytype = "ENSEMBL") %>% as.data.table
convs$iterator = 1:nrow(convs)
convs = convs[, .SD[which.min(iterator)], by = ENSEMBL] # we do this to get only one instance of each human gene
homologs = fread("~/HMD_HumanPhenotype.rpt.txt") # from MGI
homologs$iterator = 1:nrow(homologs)
homologs = homologs[, .SD[which.min(iterator)], by = V1] # we do this to get only one instance of each human gene
convs$iterator = 1:nrow(convs) #overwrite iterator to set order
convs = merge.data.table(convs,homologs,by.x="SYMBOL",by.y="V1")
setkey(convs,"iterator.x")
counts = counts[convs$iterator.x,]
rownames(counts) = convs$V3

# prep seurat object and grab md from original rds
rm(sea)
sea_mouse = CreateSeuratObject(counts, meta.data = sea@meta.data)
sea_mouse = SCTransform(sea_mouse, method = "glmGamPoi", vst.flavor = "v2")
sea_mouse = RunPCA(sea_mouse, verbose = FALSE)
sea_mouse = RunUMAP(sea_mouse, reduction = "pca", dims = 1:30, verbose = FALSE)
sea_mouse = FindNeighbors(sea_mouse, reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.5, verbose = FALSE)
(DimPlot(sea_mouse, group.by = "Supertype", label = TRUE, repel = TRUE) + NoLegend()) %>% ggsave("Rplots.pdf",.)
(DimPlot(sea_mouse, group.by = "seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()) %>% ggsave("Rplots.pdf",.)
(DimPlot(sea_mouse, group.by = "Specimen.ID") + NoLegend()) %>% ggsave("Rplots.pdf",.)

# project onto my reference
library(Signac)
library(SeuratWrappers)
library(patchwork)
seu = readRDS("~/multiome_mglia.rds")
DefaultAssay(seu) = "SCT"
# transfer.
anchors = FindTransferAnchors(reference = seu, query = sea_mouse, dims = 1:30, reference.reduction = "pca", normalization.method = "SCT")
predictions = TransferData(anchorset = anchors, refdata = seu$mglia_ident, dims = 1:30)

sea_mouse = MapQuery(anchorset = anchors, reference = seu, query = sea_mouse, refdata = "mglia_ident", 
                     reference.reduction = "pca", reduction.model = "umap")
sea_mouse$predicted.id = factor(sea_mouse$predicted.id, levels = levels(seu$mglia_ident))
p1 = DimPlot(seu, reduction = "umap", group.by = "mglia_ident", label = TRUE, label.size = 3,
             repel = TRUE) + ggtitle("Reference annotations")
p2 = DimPlot(sea_mouse, reduction = "ref.umap", group.by = "predicted.id", label = TRUE,
             label.size = 3, repel = TRUE) + ggtitle("Query transferred labels")
(p1 + p2) %>% ggsave("Rplots.pdf",.,width=14,height=7)

# save our metadata separately for easy access
fwrite(as.data.table(sea_mouse@meta.data),"sea_md.csv")
saveRDS(sea_mouse, "sea_projected.rds")

# Prep MIND ----
setwd("swarup_MIND")

library(Seurat)
library(data.table)
library(dplyr)
library(Matrix)
library(ggplot2)

mat = ReadMtx("snRNA_counts.mtx",cells="barcodes_rna.csv",features ="genes.csv",feature.column = 1,cell.column = 1)
# convert gene names to mouse homologs for mapping
homologs = fread("~/HMD_HumanPhenotype.rpt.txt") # from MGI
homologs$iterator = 1:nrow(homologs)
homologs = homologs[, .SD[which.min(iterator)], by = V1] # we do this to get only one instance of each human gene
keep.these = rownames(mat) %in% homologs$V1 %>% which
my_conv = data.table(num = 1:length(keep.these), orig = rownames(mat)[keep.these])
my_conv = merge.data.table(my_conv, homologs, by.x = "orig", by.y = "V1")
my_conv = setkey(my_conv, num)
# filter and rename
mat = mat[keep.these,]
rownames(mat) = my_conv$V3
# create seurat object
md = fread("snRNA_metadta.csv")
md_toimport = md[,-1]
rownames(md_toimport) = md$V1
mind = CreateSeuratObject(counts = mat, meta.data = md_toimport)
mind = subset(mind, subset = celltype == "MG")

# prep for query mapping
mind = SCTransform(mind, vst.flavor = "v2")
mind = RunPCA(mind)
mind = RunUMAP(mind, reduction = "pca", dims = 1:30, verbose = FALSE)
(DimPlot(mind, group.by = "cluster", label = TRUE, repel = TRUE) + NoLegend()) %>% ggsave("Rplots.pdf",.)

#transfer
library(Signac)
library(SeuratWrappers)
library(patchwork)
seu = readRDS("~/multiome_mglia.rds")

anchors = FindTransferAnchors(reference = seu, query = mind, dims = 1:30, reference.reduction = "pca", normalization.method = "SCT")
predictions = TransferData(anchorset = anchors, refdata = seu$mglia_ident, dims = 1:30)

mind = MapQuery(anchorset = anchors, reference = seu, query = mind, refdata = "mglia_ident", 
                reference.reduction = "pca", reduction.model = "umap")
mind$predicted.id = factor(mind$predicted.id, levels = levels(seu$mglia_ident))
p1 = DimPlot(seu, reduction = "umap", group.by = "mglia_ident", label = TRUE, label.size = 3,
             repel = TRUE) + ggtitle("Reference annotations")
p2 = DimPlot(mind, reduction = "ref.umap", group.by = "predicted.id", label = TRUE,
             label.size = 3, repel = TRUE) + ggtitle("Query transferred labels")
(p1 + p2) %>% ggsave("Rplots.pdf",.,width=14,height=7)

saveRDS(mind, "mind_projected.rds")
# save our metadata separately for easy access
fwrite(as.data.table(mind@meta.data),"mind_md.csv")
# Prep GSE157827 ----
setwd("GSE157827")

library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(pbapply)
library(parallel)

files = list.files("raw_data/")
prefixes = str_split(files,"_")
prefixes = lapply(prefixes, function(x){paste(x[1],x[2],sep="_")}) %>% unlist() %>% unique()
homologs = fread("~/HMD_HumanPhenotype.rpt.txt") # from MGI
homologs$iterator = 1:nrow(homologs)
homologs = homologs[, .SD[which.min(iterator)], by = V1] # we do this to get only one instance of each human gene

# write a function ro parse each sample's mtx and convert to mouse genes
seu_prepper = function(this_prefix) {
  seu = ReadMtx(mtx = paste0("raw_data/",this_prefix,"_matrix.mtx.gz"),
                cells = paste0("raw_data/",this_prefix,"_barcodes.tsv.gz"),
                features = paste0("raw_data/",this_prefix,"_features.tsv.gz"))
  keep.these = rownames(seu) %in% homologs$V1 %>% which
  my_conv = data.table(num = 1:length(keep.these), orig = rownames(seu)[keep.these])
  my_conv = merge.data.table(my_conv, homologs, by.x = "orig", by.y = "V1")
  my_conv = setkey(my_conv, num)
  # filter and rename
  seu = seu[keep.these,]
  rownames(seu) = my_conv$V3
  colnames(seu) = paste(str_split(this_prefix,"_")[[1]][2], colnames(seu), sep="_")
  return(seu)
}

cl = makeCluster(7)
clusterEvalQ(cl, lapply(c("Seurat","data.table",
                          "dplyr","stringr"),require,character.only=TRUE)) # load libraries onto each node, takes ~15s
clusterExport(cl, "homologs") # push homologs table to all nodes. ~instant
seu_list = pblapply(prefixes, seu_prepper, cl = cl)
stopCluster(cl); rm(cl)

seu_bound = do.call(cbind,seu_list)
seu = CreateSeuratObject(seu_bound)
seu = SCTransform(seu, method = "GlmGamPoi", vst.flavor = "v2")
seu = RunPCA(seu, verbose = FALSE)
seu = RunUMAP(seu, reduction = "pca", dims = 1:30, verbose = TRUE)
seu = FindNeighbors(seu, dims = 1:30)
seu = FindClusters(seu, resolution = 0.3)
(DimPlot(seu, label = TRUE, repel = TRUE) + NoLegend()) %>% ggsave("Rplots.pdf",.)
saveRDS(seu, "GSE157827_all_cells.rds")

# subset microglia
mg = subset(seu, subset = seurat_clusters == 9) # only cluster 9 has appreciable expression of microglial markers (DOCK8, P2RY12, PTPRC)
DefaultAssay(mg) = "RNA"
mg$orig.ident = str_split(Cells(mg),"_") %>% sapply(.,"[[",1)

# transfer!
library(Signac)
library(SeuratWrappers)
library(patchwork)
seu = readRDS("~/multiome_mglia.rds")
DefaultAssay(seu) = "SCT"
# transfer.
anchors = FindTransferAnchors(reference = seu, query = mg, dims = 1:30, reference.reduction = "pca", normalization.method = "SCT")
predictions = TransferData(anchorset = anchors, refdata = seu$mglia_ident, dims = 1:30)

mg = MapQuery(anchorset = anchors, reference = seu, query = mg, refdata = "mglia_ident", 
              reference.reduction = "pca", reduction.model = "umap")
mg$predicted.id = factor(mg$predicted.id, levels = levels(seu$mglia_ident))
mg$disease = substr(mg$orig.ident, 1, 2)
p1 = DimPlot(seu, reduction = "umap", group.by = "mglia_ident", label = TRUE, label.size = 3,
             repel = TRUE) + ggtitle("Reference annotations")
p2 = DimPlot(mg, reduction = "ref.umap", group.by = "predicted.id", label = TRUE,
             label.size = 3, repel = TRUE) + ggtitle("Query transferred labels")
(p1 + p2) %>% ggsave("Rplots.pdf",.,width=14,height=7)

saveRDS(mg, "mg_projected.rds")
# save our metadata separately for easy access
fwrite(as.data.table(mg@meta.data),"mg_md.csv")

# Prep Sel Vul ----
setwd("sel_vul")

library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)

sfg = readRDS("sce.SFG.Micro.scAlign.rds")
ec = readRDS("sce.EC.Micro.scAlign.rds")

sfg = as.Seurat(sfg, counts = "counts", data = "logcounts")
ec = as.Seurat(ec, counts = "counts", data = "logcounts")
selvul = merge(sfg,ec, add.cell.ids = c("sfg","ec"), project = "selvul")

# convert to mouse.
counts = selvul@assays$originalexp@counts
homologs = fread("~/HMD_HumanPhenotype.rpt.txt") # from MGI
homologs$iterator = 1:nrow(homologs)
homologs = homologs[, .SD[which.min(iterator)], by = V1] # we do this to get only one instance of each human gene
keep.these = rownames(counts) %in% homologs$V1 %>% which
my_conv = data.table(num = 1:length(keep.these), orig = rownames(counts)[keep.these])
my_conv = merge.data.table(my_conv, homologs, by.x = "orig", by.y = "V1")
my_conv = setkey(my_conv, num)
# filter and rename
counts = counts[keep.these,]
rownames(counts) = my_conv$V3
# make new seurat object
selvul_conv = CreateSeuratObject(counts, meta.data = selvul@meta.data)

# time to project.
library(Signac)
library(SeuratWrappers)
library(patchwork)
seu = readRDS("~/multiome_mglia.rds")
DefaultAssay(seu) = "SCT"
# transfer.
anchors = FindTransferAnchors(reference = seu, query = selvul_conv, dims = 1:30, reference.reduction = "pca", normalization.method = "SCT")
predictions = TransferData(anchorset = anchors, refdata = seu$mglia_ident, dims = 1:30)

selvul_conv = MapQuery(anchorset = anchors, reference = seu, query = selvul_conv, refdata = "mglia_ident", 
                       reference.reduction = "pca", reduction.model = "umap")
selvul_conv$predicted.id = factor(selvul_conv$predicted.id, levels = levels(seu$mglia_ident))
p1 = DimPlot(seu, reduction = "umap", group.by = "mglia_ident", label = TRUE, label.size = 3,
             repel = TRUE) + ggtitle("Reference annotations")
p2 = DimPlot(selvul_conv, reduction = "ref.umap", group.by = "predicted.id", label = TRUE,
             label.size = 3, repel = TRUE) + ggtitle("Query transferred labels")
(p1 + p2) %>% ggsave("Rplots.pdf",.,width=14,height=7)

saveRDS(selvul_conv, "selvul_projected.rds")
# save our metadata separately for easy access
fwrite(as.data.table(selvul_conv@meta.data),"selvul_md.csv")

# Prep Tsai ---- 
setwd("tsai_myelin")

library(Seurat)
library(data.table)
library(dplyr)
library(Matrix)
library(ggplot2)

cells = fread("qc_column_metadata.csv")
genes = fread("qc_gene_names.txt", header = F)
A = readMM("qc_counts.mtx")
rownames(A)=genes$V1
colnames(A)=cells$V1

# convert gene names to mouse homologs for mapping
homologs = fread("~/HMD_HumanPhenotype.rpt.txt") # from MGI
homologs$iterator = 1:nrow(homologs)
homologs = homologs[, .SD[which.min(iterator)], by = V1] # we do this to get only one instance of each human gene
keep.these = genes$V1 %in% homologs$V1 %>% which
my_conv = data.table(num = 1:length(keep.these), orig = genes$V1[keep.these])
my_conv = merge.data.table(my_conv, homologs, by.x = "orig", by.y = "V1")
my_conv = setkey(my_conv, num)
# filter and rename
A = A[keep.these,]
rownames(A) = my_conv$V3
# let's also prefilter the sparsematrix to only microglia
mglia_indices = which(cells$cell.type == "Mic")
A_mg = A[,mglia_indices]
# create seurat object
cells_asmd = cells[mglia_indices,-1]
rownames(cells_asmd) = cells[mglia_indices,V1]
tsai = CreateSeuratObject(counts = A_mg, meta.data = cells_asmd)
rm(A,A_mg,cells,cells_asmd,genes,homologs,my_conv)

# prep for query mapping
tsai = SCTransform(tsai, method = "glmGamPoi", vst.flavor = "v2")
tsai = RunPCA(tsai, verbose = FALSE)
tsai = RunUMAP(tsai, reduction = "pca", dims = 1:30, verbose = TRUE)
(DimPlot(tsai, group.by = "apoe_genotype", label = TRUE, repel = TRUE)) %>% ggsave("Rplots.pdf",.)

# now we get our own multiome data so we can project sn -> sn
library(Signac)
library(SeuratWrappers)
library(patchwork)
seu = readRDS("~/multiome_mglia.rds")
DefaultAssay(seu) = "SCT"
# transfer.
anchors = FindTransferAnchors(reference = seu, query = tsai, dims = 1:30, reference.reduction = "pca", normalization.method = "SCT")
predictions = TransferData(anchorset = anchors, refdata = seu$mglia_ident, dims = 1:30)

tsai = MapQuery(anchorset = anchors, reference = seu, query = tsai, refdata = "mglia_ident", 
                reference.reduction = "pca", reduction.model = "umap")
tsai$predicted.id = factor(tsai$predicted.id, levels = levels(seu$mglia_ident))
p1 = DimPlot(seu, reduction = "umap", group.by = "mglia_ident", label = TRUE, label.size = 3,
             repel = TRUE) + ggtitle("Reference annotations")
p2 = DimPlot(tsai, reduction = "ref.umap", group.by = "predicted.id", label = TRUE,
             label.size = 3, repel = TRUE) + ggtitle("Query transferred labels")
(p1 + p2) %>% ggsave("Rplots.pdf",.,width=14,height=7)

saveRDS(tsai, "tsai_projected.rds")
# save our metadata separately for easy access
fwrite(as.data.table(tsai@meta.data),"tsai_md.csv")

# Prep Ento ----
library(Seurat)
library(data.table)
library(dplyr)
library(Matrix)
library(ggplot2)

setwd("entorhinal")

dat = fread("scRNA_rawCounts.tsv")
mat = dat[,-1] %>% as.matrix
rownames(mat) = dat$geneName
# convert gene names to mouse homologs for mapping
homologs = fread("~/HMD_HumanPhenotype.rpt.txt") # from MGI
homologs$iterator = 1:nrow(homologs)
homologs = homologs[, .SD[which.min(iterator)], by = V1] # we do this to get only one instance of each human gene
keep.these = rownames(mat) %in% homologs$V1 %>% which
my_conv = data.table(num = 1:length(keep.these), orig = rownames(mat)[keep.these])
my_conv = merge.data.table(my_conv, homologs, by.x = "orig", by.y = "V1")
my_conv = setkey(my_conv, num)
# filter and rename
mat = mat[keep.these,]
rownames(mat) = my_conv$V3
# create seurat object
md = fread("scRNA_metadata.tsv")
md_toimport = md[,-1]
rownames(md_toimport) = md$sampleID
ento = CreateSeuratObject(counts = mat, meta.data = md_toimport)
ento = subset(ento, subset = cellType == "mg")
ento = SCTransform(ento, method = "glmGamPoi", vst.flavor = "v2")

#transfer
library(Signac)
library(SeuratWrappers)
library(patchwork)
seu = readRDS("~/multiome_mglia.rds")

anchors = FindTransferAnchors(reference = seu, query = ento, dims = 1:30, reference.reduction = "pca", normalization.method = "SCT")
predictions = TransferData(anchorset = anchors, refdata = seu$mglia_ident, dims = 1:30)

ento = MapQuery(anchorset = anchors, reference = seu, query = ento, refdata = "mglia_ident", 
                reference.reduction = "pca", reduction.model = "umap")
ento$predicted.id = factor(ento$predicted.id, levels = levels(seu$mglia_ident))
p1 = DimPlot(seu, reduction = "umap", group.by = "mglia_ident", label = TRUE, label.size = 3,
             repel = TRUE) + ggtitle("Reference annotations")
p2 = DimPlot(ento, reduction = "ref.umap", group.by = "predicted.id", label = TRUE,
             label.size = 3, repel = TRUE) + ggtitle("Query transferred labels")
(p1 + p2) %>% ggsave("Rplots.pdf",.,width=14,height=7)

saveRDS(ento, "ento_projected.rds")
# save our metadata separately for easy access
fwrite(as.data.table(ento@meta.data),"ento_md.csv")



# Generate plots ----
library(data.table)
library(dplyr)
library(ggplot2)

mglia_levels = c("Homeostatic Microglia","Selplg-lo Microglia","mt-Enriched Microglia","mt-Depleted Microglia",
                 "DAM-1","DAM-2","TIMs","Siglech-hi Microglia")

seu_md = fread("multiome_md.csv")
rosmap_md = fread("ROSMAP_WholeBrain/rosmap_wholebrain_md.csv")
BA10_md = fread("snRNAseqPFC_BA10/BA10_md.csv")
sea_md = fread("SEA-AD/sea_md.csv")
mind_md = fread("swarup_MIND/mind_md.csv")
gse157827_md = fread("GSE157827/mg_md.csv")
selvul_md = fread("sel_vul/selvul_md.csv")
tsai_md = fread("tsai_myelin/tsai_md.csv")
ento_md = fread("entorhinal/ento_md.csv")

freqs_seu = seu_md[, .N/nrow(seu_md), by = "mglia_ident"]
freqs_rosmap = rosmap_md[, .N/nrow(rosmap_md), by = "predicted.id"]
freqs_ba10 = BA10_md[, .N/nrow(BA10_md), by = "predicted.id"]
freqs_sea = sea_md[, .N/nrow(sea_md), by = "predicted.id"]
freqs_mind = mind_md[, .N/nrow(mind_md), by = "predicted.id"]
freqs_gse157827 = gse157827_md[, .N/nrow(gse157827_md), by = "predicted.id"]
freqs_selvul = selvul_md[, .N/nrow(selvul_md), by = "predicted.id"]
freqs_tsai = tsai_md[, .N/nrow(tsai_md), by = "predicted.id"]
freqs_ento = ento_md[, .N/nrow(ento_md), by = "predicted.id"]

colnames(freqs_seu) = c("mglia_ident","Multiome Reference\n(This Study)")
colnames(freqs_rosmap) = c("mglia_ident","ROSMAP DLPFC") #https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyDetails?Study=syn3219045, method: https://protocolexchange.researchsquare.com/article/nprot-6163/v1
colnames(freqs_ba10) = c("mglia_ident","Mathys, 2019") #https://www.nature.com/articles/s41586-019-1195-2
colnames(freqs_sea) = c("mglia_ident","SEA-AD") #https://portal.brain-map.org/explore/seattle-alzheimers-disease
colnames(freqs_mind) = c("mglia_ident","Morabito, 2021") # https://www.nature.com/articles/s41588-021-00894-z
colnames(freqs_gse157827) = c("mglia_ident","Lau, 2020") # https://www.pnas.org/doi/abs/10.1073/pnas.2008762117
colnames(freqs_selvul) = c("mglia_ident","Leng, 2021") # https://www.nature.com/articles/s41593-020-00764-7
colnames(freqs_tsai) = c("mglia_ident","Blanchard, 2022") # https://www.nature.com/articles/s41586-022-05439-w
colnames(freqs_ento) = c("mglia_ident","Grubman, 2019") #https://www.nature.com/articles/s41593-019-0539-4

freqs = Reduce(merge,list(freqs_seu,freqs_rosmap,freqs_ba10,freqs_sea,freqs_mind,freqs_gse157827,freqs_selvul,freqs_tsai,freqs_ento)) %>% 
  melt.data.table(., id.vars = "mglia_ident") %>% 
  mutate(variable = factor(variable,levels = c("Multiome Reference\n(This Study)","ROSMAP DLPFC","Mathys, 2019",
                                               "SEA-AD","Morabito, 2021","Lau, 2020","Leng, 2021","Blanchard, 2022",
                                               "Grubman, 2019"))) %>% 
  mutate(mglia_ident = factor(mglia_ident, levels = mglia_levels))

(ggplot(freqs, aes(x = mglia_ident, y = 100*value, fill = variable)) +
    geom_bar(position = "dodge", stat = "identity") + 
    labs(x = "Cluster", y = "Frequency (%)", fill = "Data Source") + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    theme(plot.margin = margin(5,5,5,20))) %>% 
  ggsave("transfer_freqs.png",.,height=7,width=10,dpi=600)

timsonly = freqs[mglia_ident == "TIMs"]
timsonly$method = factor(c("Enzymatic","Cold Dounce","Cold Dounce",
                           "Enzymatic","Cold Dounce","Cold Dounce",
                           "Cold Dounce","Cold Dounce","Cold Dounce"), levels = c("Enzymatic","Cold Dounce"))
## Fig. 4B ----
plot = timsonly %>% 
  mutate(variable = forcats::fct_reorder(variable,-value)) %>%
  ggplot(aes(x = variable, y = 100*value, fill = variable, pattern = method)) +
  geom_bar(stat = "identity", color = "black") + 
  labs(x = "", y = "Frequency of TIMs (%)") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  theme(plot.margin = margin(5,5,5,20)) + 
  facet_grid(~method, scales = "free_x", space = "free") +
  theme(legend.position = "none")

ggsave("transfer_freqs_timsonly.png",plot,height=4.5,width=4.5,dpi=600)

## now let's make all our projections
library(Seurat)
library(Signac)
library(patchwork)

seu = readRDS("~/scratch/R_dir/1yr_multiome/multiome_mglia.rds")
rosmap = readRDS("ROSMAP_WholeBrain/rosmap_projected.rds")
BA10 = readRDS("snRNAseqPFC_BA10/BA10_projected.rds")
sea_mouse = readRDS("SEA-AD/sea_projected.rds")
mind = readRDS("swarup_MIND/mind_projected.rds")
gse157827 = readRDS("GSE157827/mg_projected.rds")
selvul = readRDS("sel_vul/selvul_projected.rds")
tsai = readRDS("tsai_myelin/tsai_projected.rds")
ento = readRDS("entorhinal/ento_projected.rds")

p1 = DimPlot(seu, reduction = "umap", group.by = "mglia_ident") + ggtitle("Multiome Reference\n(This Study)") + xlim(-6,6) + ylim(-6,6)
p2 = DimPlot(rosmap, reduction = "ref.umap", group.by = "predicted.id") + ggtitle("ROSMAP DLPFC") + xlim(-6,6) + ylim(-6,6)
p3 = DimPlot(BA10, reduction = "ref.umap", group.by = "predicted.id") + ggtitle("Mathys, 2019") + xlim(-6,6) + ylim(-6,6)
p4 = DimPlot(sea_mouse, reduction = "ref.umap", group.by = "predicted.id") + ggtitle("SEA-AD") + xlim(-6,6) + ylim(-6,6)
p5 = DimPlot(mind, reduction = "ref.umap", group.by = "predicted.id") + ggtitle("Morabito, 2021") + xlim(-6,6) + ylim(-6,6)
p6 = DimPlot(gse157827, reduction = "ref.umap", group.by = "predicted.id") + ggtitle("Lau, 2020") + xlim(-6,6) + ylim(-6,6)
p7 = DimPlot(selvul, reduction = "ref.umap", group.by = "predicted.id") + ggtitle("Leng, 2021") + xlim(-6,6) + ylim(-6,6)
p8 = DimPlot(tsai, reduction = "ref.umap", group.by = "predicted.id") + ggtitle("Blanchard, 2022") + xlim(-6,6) + ylim(-6,6)
p9 = DimPlot(ento, reduction = "ref.umap", group.by = "predicted.id") + ggtitle("Grubman, 2019") + xlim(-6,6) + ylim(-6,6)

## Fig. 4A ----
wrap_plots(list(p1,p4,p9,p7,p2,p6,p3,p5,p8)) %>% ggsave("projections.png",.,width=14,height=10,dpi=600)

# specific metadata plots for supplement
# from Lau, 2020
## Fig. 4C ----
setkey(gse157827_md, "predicted.id","orig.ident") 
dat = gse157827_md[CJ(predicted.id, orig.ident, unique = TRUE),.N,by=.EACHI]
dat[, frac := N/sum(N) * 100, by = orig.ident]
dat = dat[predicted.id == "TIMs" & frac < 7]
dat$disease = substr(dat$orig.ident, 1, 2)
t.test(dat$frac ~ dat$disease)
dat$disease = ifelse(dat$disease == "AD", "Alzheimer's","Normal Control")
(ggplot(dat, aes(x = disease, y = frac, fill = disease)) +
    geom_boxplot() + 
    labs(x = "", y = "Frequency of TIMs (%)") + 
    theme_bw() + 
    #theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    NoLegend() + 
    ggtitle("Lau, 2020: TIM Frequency\nby Disease State")) %>% 
  ggsave("lau_2020_disease_state.png",.,height=4,width=4,dpi=600)

# from sea-ad
## Fig. 4F ----
setkey(sea_md, "predicted.id","Specimen.ID")
dat = sea_md[CJ(predicted.id, Specimen.ID, unique = TRUE),.N,by=.EACHI]
dat[, frac := N/sum(N) * 100, by = Specimen.ID]
col = "PMI"
md_to_add = sea_md[,.N,by=c(col,"Specimen.ID")]
dat = merge.data.table(dat,md_to_add,by = "Specimen.ID") #N.x is number of cells in predicted.id, N.y is total number of cells
toplot = dat[,mean(frac),by=c(col,"predicted.id")]
count_per = sea_md[,length(unique(.SD[["Specimen.ID"]])),by=col]
toplot = merge.data.table(toplot, count_per, by = col)
toplot = toplot %>%
  mutate(UQ(rlang::sym(col)) := factor(UQ(rlang::sym(col)), levels = levels(sea_mouse[[col]][,1])))
toplot = toplot[predicted.id == "TIMs"][order(match(get(col),levels(sea_mouse[[col]][,1]))),]
toplot$V1.y = paste0("N = ",toplot$V1.y)
(ggplot(toplot,aes(x = UQ(rlang::sym(col)), y = V1.x, fill = UQ(rlang::sym(col)))) + 
    geom_bar(stat = "identity", color = "black") +
    geom_text(aes(label = V1.y), vjust = -1) + 
    theme_bw() + 
    theme(axis.title.x = element_blank()) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(stringr::str_replace(x,"-"," "), width = 10)) +
    labs(y = "Frequency of TIMs (%)", title = "SEA-AD: TIM Frequency by PMI") + 
    ylim(0,2.5) +
    NoLegend()) %>% 
  ggsave("SEA_tims_by_pmi.png",.,height = 4, width = 4, dpi = 600)

# from Leng, 2021
## Fig. 4D ----
setkey(selvul_md, "predicted.id","SampleID") 
dat = selvul_md[CJ(predicted.id, SampleID, unique = TRUE),.N,by=.EACHI]
dat[, frac := N/sum(N) * 100, by = SampleID]
dat = dat[predicted.id == "TIMs" & frac < 10]
braaks = selvul_md[,c("BraakStage","SampleID")] %>% unique
dat = merge.data.table(dat, braaks, by = "SampleID")
dat$disease = ifelse(dat$BraakStage > 0, "Braak High","Braak Low")
t.test(dat$frac ~ dat$disease)

(ggplot(dat, aes(x = disease, y = frac, fill = disease)) +
    geom_boxplot() + 
    labs(x = "", y = "Frequency of TIMs (%)") + 
    theme_bw() + 
    #theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    NoLegend() + 
    ggtitle("Leng, 2021: TIM Frequency\nby Braak Stage")) %>% 
  ggsave("leng_2021_tims_by_braak.png",.,height = 4, width = 4, dpi = 600)

# from blanchard, 2022
## Fig. 4F ----
tsai_md$is_tim = ifelse(tsai_md$predicted.id == "TIMs", TRUE, FALSE)
tsai_md$amyloid_hi = ifelse(tsai_md$amyloid > 10, "yes", "no")

dat_amyloidhi = tsai_md[, .N, by = .(is_tim,amyloid_hi)]
dat_amyloidhi = dcast.data.table(dat_amyloidhi, amyloid_hi ~ is_tim, value.var = "N")
setnafill(dat_amyloidhi, fill = 0)
dat_amyloidhi$frac = (100*dat_amyloidhi$`TRUE`/(dat_amyloidhi$`TRUE`+dat_amyloidhi$`FALSE`))

dat_gen = tsai_md[, .N, by = .(is_tim,APOE4)]
dat_gen = dcast.data.table(dat_gen, APOE4 ~ is_tim, value.var = "N")
setnafill(dat_gen, fill = 0)
dat_gen$frac = (100*dat_gen$`TRUE`/(dat_gen$`TRUE`+dat_gen$`FALSE`))

toplot = data.table(lab = c("Amyloid-hi?","Amyloid-hi?","APOE4?","APOE4?"),
                    case = stringr::str_to_title(c(dat_amyloidhi$amyloid_hi,dat_gen$APOE4)),
                    val = c(dat_amyloidhi$frac,dat_gen$frac))

(ggplot(toplot, aes(x = lab, y = val, fill = case)) +
    geom_bar(position = "dodge", stat = "identity", color = "black") + 
    labs(x = "", y = "Frequency of TIMs (%)", fill = "Status") + 
    theme_bw() + 
    #theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    ggtitle("Blanchard, 2022: TIM Frequency by\nAmyloid Burden and APOE Isoform")) %>% 
  ggsave("blanchard_2022_tims_by_amyloid_and_apoe.png",.,height = 4, width = 4, dpi = 600)

# Tabula muris senis ----
library(Seurat)
library(dplyr)
library(data.table)
library(ggplot2)
library(org.Mm.eg.db)

int_mglia = readRDS("~/mglia_only.rds")
tabula = readRDS("brain_myeloid.rds")
tabula = subset(tabula, subset = cell_type == "microglial cell")
tabula_counts = tabula@assays$RNA@counts
symbols = mapIds(org.Mm.eg.db,keys=rownames(tabula_counts),column="SYMBOL", keytype="ENSEMBL", multiVals="first")
rownames(tabula_counts) = symbols
tabula_counts = tabula_counts[-which(is.na(symbols)),]
tabula_rename = CreateSeuratObject(tabula_counts, meta.data = tabula@meta.data)
tabula_rename = SCTransform(tabula_rename, method = "glmGamPoi", vst.flavor = "v2")

DefaultAssay(int_mglia) = "SCT"
int_mglia = RunUMAP(int_mglia, dims = 1:46, verbose = TRUE, return.model = TRUE) #rerun to store model

anchors = FindTransferAnchors(reference = int_mglia, query = tabula_rename, dims = 1:30, reference.reduction = "pca", normalization.method = "SCT")
predictions = TransferData(anchorset = anchors, refdata = int_mglia$mglia_ident, dims = 1:30)

tabula_rename = MapQuery(anchorset = anchors, reference = int_mglia, query = tabula_rename, refdata = "mglia_ident", 
                         reference.reduction = "pca", reduction.model = "umap")
tabula_rename$predicted.id = factor(tabula_rename$predicted.id, levels = levels(int_mglia$mglia_ident))
p1 = DimPlot(int_mglia, reduction = "umap", group.by = "mglia_ident", label = TRUE, label.size = 3,
             repel = TRUE) + ggtitle("Reference annotations")
p2 = DimPlot(tabula_rename, reduction = "ref.umap", group.by = "predicted.id", label = TRUE,
             label.size = 3, repel = TRUE) + ggtitle("Query transferred labels")
(p1 + p2) %>% ggsave("Rplots.pdf",.,width=20,height=7)

## Fig. 4G ----
md = tabula_rename@meta.data %>% as.data.table
md[, istim := grepl("TIM",predicted.id,fixed=TRUE)]
md[, regions_of_interest := grepl("Cortex|Hippocampus",subtissue)]
plot = md[regions_of_interest == TRUE, .N, by = .(age,istim)] %>%
  dcast.data.table(age ~ istim, value.var = "N") %>%
  mutate(frac = `TRUE`/(`TRUE`+`FALSE`)) %>% 
  mutate(frac = 100*frac) %>%
  mutate(age = as.character(age)) %>% 
  mutate(age = paste0(age,"o")) %>% 
  mutate(age = factor(age, levels = c("3mo","18mo","24mo"))) %>% 
  ggplot(aes(x=age,y=frac,fill=age)) + 
  geom_bar(stat = "identity",color="black") + 
  theme_bw() + 
  theme(legend.position = "none") +
  labs(x="",y="Frequency of TIMs (%)") + 
  ggtitle("Tabula Muris Senis, 2020")
ggsave("tim_freq_by_age.png",plot,dpi=600,width=4,height=4)

saveRDS(tabula_rename, "tabula_projected.rds")

