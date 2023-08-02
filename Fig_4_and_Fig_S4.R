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

# Prep Prater ----
setwd("PU1")

library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)

load("Prater_Green_PU1_MGsubset_10clusters_DeID.rdata")
pu1 = ss_data_norm
rm(ss_data_norm)

# convert to mouse.
counts = pu1@assays$RNA@counts
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
pu1_conv = CreateSeuratObject(counts, meta.data = pu1@meta.data)

# time to project.
library(Signac)
library(SeuratWrappers)
library(patchwork)
seu = readRDS("~/multiome_mglia.rds")
DefaultAssay(seu) = "SCT"
# transfer.
anchors = FindTransferAnchors(reference = seu, query = pu1_conv, dims = 1:30, reference.reduction = "pca", normalization.method = "SCT")
predictions = TransferData(anchorset = anchors, refdata = seu$mglia_ident, dims = 1:30)

pu1_conv = MapQuery(anchorset = anchors, reference = seu, query = pu1_conv, refdata = "mglia_ident", 
                    reference.reduction = "pca", reduction.model = "umap")
pu1_conv$predicted.id = factor(pu1_conv$predicted.id, levels = levels(seu$mglia_ident))
p1 = DimPlot(seu, reduction = "umap", group.by = "mglia_ident", label = TRUE, label.size = 3,
             repel = TRUE) + ggtitle("Reference annotations")
p2 = DimPlot(pu1_conv, reduction = "ref.umap", group.by = "predicted.id", label = TRUE,
             label.size = 3, repel = TRUE) + ggtitle("Query transferred labels")
(p1 + p2) %>% ggsave("transfer.png",.,width=14,height=7)

pu1_md = pu1_conv@meta.data %>% as.data.table
seu_md = seu@meta.data %>% as.data.table

freqs_pu1 = pu1_md[, .N/nrow(pu1_md), by = "predicted.id"]
freqs_seu = seu_md[, .N/nrow(seu_md), by = "mglia_ident"]
freqs = merge.data.table(x = freqs_pu1, y = freqs_seu, by.x = "predicted.id", by.y = "mglia_ident") %>% 
  melt.data.table(., id.vars = "predicted.id", measure.vars = c("V1.x", "V1.y"))

(ggplot(freqs, aes(x = predicted.id, y = value, fill = variable)) +
    geom_bar(position = "dodge", stat = "identity") + 
    scale_fill_discrete("Data Source",labels = c("PU1","1yr_Multiome")) + 
    labs(x = "Cluster", y = "Frequency") + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))) %>% 
  ggsave("transferfreqs.png",.)

saveRDS(pu1_conv, "pu1_projected.rds")
# save our metadata separately for easy access
fwrite(as.data.table(pu1_conv@meta.data),"pu1_md.csv")

# Prep ROSMAP DLPFC-2 ----
library(Seurat)
library(dplyr)
library(data.table)

setwd("DLPFC_2")

prefixes = list.files("mtx_files") %>% strsplit(.,".",fixed=TRUE) %>% sapply(.,"[[",1) %>% unique
dir = paste0(getwd(),"/mtx_files/")

md = fread("cell-annotation.csv")

for (i in 101:127) { # run separately in intervals of 20
  tmp = ReadMtx(mtx = paste0(dir,prefixes[i],".matrix.mtx.gz"),
                cells = paste0(dir,prefixes[i],".barcodes.tsv.gz"),
                features = paste0(dir,prefixes[i],".features.tsv.gz"))
  tmp = CreateSeuratObject(tmp)
  tmp$identifier = prefixes[i]
  tmp = subset(tmp, cells = md[libraryBatch == prefixes[i],cellBarcode])
  tmp = RenameCells(tmp, add.cell.id = prefixes[i])
  if(i == 101){
    merge = tmp
    print(paste0(i,"/127 Done!"))
  } else {
    merge = merge(x = merge, y = tmp)
    print(paste0(i,"/127 Done!"))
  }
}
saveRDS(merge, "tmpmerge_101-127.rds")

# after importing:
rownames(md) = paste0(md$libraryBatch,"_",md$cellBarcode)

seu1 = readRDS("tmpmerge_1-20.rds")
seu1 = AddMetaData(seu1, md)
seu1 = subset(seu1, subset = cell.type == "Microglia")

seu2 = readRDS("tmpmerge_21-40.rds")
seu2 = AddMetaData(seu2, md)
seu2 = subset(seu2, subset = cell.type == "Microglia")
seu = merge(seu1,seu2)

seu3 = readRDS("tmpmerge_41-60.rds")
seu3 = AddMetaData(seu3, md)
seu3 = subset(seu3, subset = cell.type == "Microglia")
seu = merge(seu,seu3)

seu4 = readRDS("tmpmerge_61-80.rds")
seu4 = AddMetaData(seu4, md)
seu4 = subset(seu4, subset = cell.type == "Microglia")
seu = merge(seu,seu4)

seu5 = readRDS("tmpmerge_81-100.rds")
seu5 = AddMetaData(seu5, md)
seu5 = subset(seu5, subset = cell.type == "Microglia")
seu = merge(seu,seu5)

seu6 = readRDS("tmpmerge_101-127.rds")
seu6 = AddMetaData(seu6, md)
seu6 = subset(seu6, subset = cell.type == "Microglia")
seu = merge(seu,seu6)

rm(seu1,seu2,seu3,seu4,seu5,seu6)
saveRDS(seu, "dlpfc_mglia.rds")

# convert to homologs
counts = GetAssayData(seu, slot = "counts")
convs = data.table(human = rownames(counts))
convs$iterator = 1:nrow(convs)
convs = convs[, .SD[which.min(iterator)], by = human] # we do this to get only one instance of each human gene
homologs = fread("~/HMD_HumanPhenotype.rpt.txt") # from MGI
homologs$iterator = 1:nrow(homologs)
homologs = homologs[, .SD[which.min(iterator)], by = V1] # we do this to get only one instance of each human gene
convs$iterator = 1:nrow(convs) #overwrite iterator to set order
convs = merge.data.table(convs,homologs,by.x="human",by.y="V1")
convs = convs[, .SD[which.min(iterator.x)], by = V3] # we do this to get only one instance of each mouse gene
setkey(convs,"iterator.x")
counts = counts[convs$iterator.x,]
rownames(counts) = convs$V3

# prep seurat object and grab md from original rds
library(ggplot2)
seu_mouse = CreateSeuratObject(counts, meta.data = seu@meta.data)
seu_mouse = SCTransform(seu_mouse, method = "glmGamPoi", vst.flavor = "v2")
seu_mouse = RunPCA(seu_mouse, verbose = FALSE)
seu_mouse = RunUMAP(seu_mouse, reduction = "pca", dims = 1:30, verbose = FALSE)
seu_mouse = FindNeighbors(seu_mouse, reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.5, verbose = FALSE)
(DimPlot(seu_mouse, group.by = "state", label = TRUE, repel = TRUE) + NoLegend()) %>% ggsave("Rplots.pdf",.)
(DimPlot(seu_mouse, group.by = "seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()) %>% ggsave("Rplots.pdf",.)
(DimPlot(seu_mouse, group.by = "individualID") + NoLegend()) %>% ggsave("Rplots.pdf",.)
(DimPlot(seu_mouse, group.by = "libraryBatch") + NoLegend()) %>% ggsave("Rplots.pdf",.)
# looks well integrated!

# transfer
library(Signac)
library(SeuratWrappers)
library(patchwork)
seu = readRDS("~/multiome_mglia.rds")

anchors = FindTransferAnchors(reference = seu, query = seu_mouse, dims = 1:30, reference.reduction = "pca", normalization.method = "SCT")
predictions = TransferData(anchorset = anchors, refdata = seu$mglia_ident, dims = 1:30)

seu_mouse = MapQuery(anchorset = anchors, reference = seu, query = seu_mouse, refdata = "mglia_ident", 
                     reference.reduction = "pca", reduction.model = "umap")
seu_mouse$predicted.id = factor(seu_mouse$predicted.id, levels = levels(seu$mglia_ident))
p1 = DimPlot(seu, reduction = "umap", group.by = "mglia_ident", label = TRUE, label.size = 3,
             repel = TRUE) + ggtitle("Reference annotations")
p2 = DimPlot(seu_mouse, reduction = "ref.umap", group.by = "predicted.id", label = TRUE,
             label.size = 3, repel = TRUE) + ggtitle("Query transferred labels")
(p1 + p2) %>% ggsave("Rplots.pdf",.,width=14,height=7)

dlpfc_md = seu_mouse@meta.data %>% as.data.table
seu_md = seu@meta.data %>% as.data.table

freqs_dlpfc = dlpfc_md[, .N/nrow(dlpfc_md), by = "predicted.id"]
freqs_seu = seu_md[, .N/nrow(seu_md), by = "mglia_ident"]
freqs = merge.data.table(x = freqs_dlpfc, y = freqs_seu, by.x = "predicted.id", by.y = "mglia_ident") %>% 
  melt.data.table(., id.vars = "predicted.id", measure.vars = c("V1.x", "V1.y"))

(ggplot(freqs, aes(x = predicted.id, y = value, fill = variable)) +
    geom_bar(position = "dodge", stat = "identity") + 
    scale_fill_discrete("Data Source",labels = c("DLPFC","1yr_Multiome")) + 
    labs(x = "Cluster", y = "Frequency") + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))) %>% 
  ggsave("Rplots.pdf",.)

saveRDS(seu_mouse, "dlpfc_mglia_projected.rds")
# save our metadata separately for easy access
fwrite(as.data.table(seu_mouse@meta.data),"dlpfc_md.csv")

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
dlpfc2_md = fread("DLPFC_2/dlpfc_md.csv")
pu1_md = fread("PU1/pu1_md.csv")

freqs_seu = seu_md[, .N/nrow(seu_md), by = "mglia_ident"]
freqs_rosmap = rosmap_md[, .N/nrow(rosmap_md), by = "predicted.id"]
freqs_ba10 = BA10_md[, .N/nrow(BA10_md), by = "predicted.id"]
freqs_sea = sea_md[, .N/nrow(sea_md), by = "predicted.id"]
freqs_mind = mind_md[, .N/nrow(mind_md), by = "predicted.id"]
freqs_gse157827 = gse157827_md[, .N/nrow(gse157827_md), by = "predicted.id"]
freqs_selvul = selvul_md[, .N/nrow(selvul_md), by = "predicted.id"]
freqs_tsai = tsai_md[, .N/nrow(tsai_md), by = "predicted.id"]
freqs_ento = ento_md[, .N/nrow(ento_md), by = "predicted.id"]
freqs_dlpfc2 = dlpfc2_md[, .N/nrow(dlpfc2_md), by = "predicted.id"]
freqs_pu1 = pu1_md[, .N/nrow(pu1_md), by = "predicted.id"]

colnames(freqs_seu) = c("mglia_ident","Multiome Reference\n(This Study)")
colnames(freqs_rosmap) = c("mglia_ident","ROSMAP DLPFC") #https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyDetails?Study=syn3219045, method: https://protocolexchange.researchsquare.com/article/nprot-6163/v1
colnames(freqs_ba10) = c("mglia_ident","Mathys, 2019") #https://www.nature.com/articles/s41586-019-1195-2
colnames(freqs_sea) = c("mglia_ident","SEA-AD") #https://portal.brain-map.org/explore/seattle-alzheimers-disease
colnames(freqs_mind) = c("mglia_ident","Morabito, 2021") # https://www.nature.com/articles/s41588-021-00894-z
colnames(freqs_gse157827) = c("mglia_ident","Lau, 2020") # https://www.pnas.org/doi/abs/10.1073/pnas.2008762117
colnames(freqs_selvul) = c("mglia_ident","Leng, 2021") # https://www.nature.com/articles/s41593-020-00764-7
colnames(freqs_tsai) = c("mglia_ident","Blanchard, 2022") # https://www.nature.com/articles/s41586-022-05439-w
colnames(freqs_ento) = c("mglia_ident","Grubman, 2019") #https://www.nature.com/articles/s41593-019-0539-4
colnames(freqs_dlpfc2) = c("mglia_ident","ROSMAP DLPFC-2") #https://www.synapse.org/#!Synapse:syn31512863
colnames(freqs_pu1) = c("mglia_ident","Prater, 2023") #https://www.nature.com/articles/s43587-023-00424-y

str_wrap_factor = function(x,width) {
  levels(x) = str_wrap(levels(x), width=width)
  x
}

freqs = Reduce(merge,list(freqs_seu,freqs_rosmap,freqs_ba10,freqs_sea,freqs_mind,freqs_gse157827,
                          freqs_selvul,freqs_tsai,freqs_ento,freqs_dlpfc2, freqs_pu1)) %>% 
  melt.data.table(., id.vars = "mglia_ident") %>% 
  mutate(variable = factor(variable,levels = c("Multiome Reference (This Study)","ROSMAP DLPFC-1","Mathys, 2019",
                                               "SEA-AD","Morabito, 2021","Lau, 2020","Leng, 2021","Blanchard, 2022",
                                               "Grubman, 2019","ROSMAP DLPFC-2","Prater, 2023"))) %>% 
  mutate(variable = str_wrap_factor(variable,width=5)) %>%
  mutate(mglia_ident = factor(mglia_ident, levels = mglia_levels))

(ggplot(freqs, aes(x = mglia_ident, y = 100*value, fill = variable)) +
    geom_bar(position = "dodge", stat = "identity") + 
    labs(x = "Cluster", y = "Frequency (%)", fill = "Data Source") + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    theme(plot.margin = margin(5,5,5,20))) %>% 
  ggsave("transfer_freqs.png",.,height=7,width=10,dpi=600)

## Fig. 4B ----
timsonly = freqs[mglia_ident == "TIMs"]
timsonly$method = factor(c("Enzymatic","Cold Dounce","Cold Dounce",
                           "Enzymatic","Cold Dounce","Cold Dounce",
                           "Cold Dounce","Cold Dounce","Cold Dounce",
                           "Cold Dounce", "Enzymatic"), levels = c("Enzymatic","Cold Dounce"))
plot = timsonly %>% 
  mutate(variable = forcats::fct_reorder(variable,value)) %>%
  ggplot(aes(x = variable, y = 100*value, fill = variable)) +
  geom_bar(stat = "identity", color = "black") + 
  labs(x = "", y = "Frequency of TIMs (%)") + 
  theme_bw(base_size = 16) + 
  facet_grid(rows = vars(method), scales = "free_y", space = "free") +
  coord_flip() +
  theme(legend.position = "none")

ggsave("transfer_freqs_timsonly.png",plot,height=12,width=5.75,dpi=600)


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
dlpfc2 = readRDS("DLPFC_2/dlpfc_mglia_projected.rds")
pu1 = readRDS("PU1/pu1_projected.rds")

p1 = DimPlot(seu, reduction = "umap", group.by = "mglia_ident") + ggtitle("Multiome Reference\n(This Study)") + xlim(-6,6) + ylim(-6,6)
p2 = DimPlot(rosmap, reduction = "ref.umap", group.by = "predicted.id") + ggtitle("ROSMAP DLPFC-1") + xlim(-6,6) + ylim(-6,6)
p3 = DimPlot(BA10, reduction = "ref.umap", group.by = "predicted.id") + ggtitle("Mathys, 2019") + xlim(-6,6) + ylim(-6,6)
p4 = DimPlot(sea_mouse, reduction = "ref.umap", group.by = "predicted.id") + ggtitle("SEA-AD") + xlim(-6,6) + ylim(-6,6)
p5 = DimPlot(mind, reduction = "ref.umap", group.by = "predicted.id") + ggtitle("Morabito, 2021") + xlim(-6,6) + ylim(-6,6)
p6 = DimPlot(gse157827, reduction = "ref.umap", group.by = "predicted.id") + ggtitle("Lau, 2020") + xlim(-6,6) + ylim(-6,6)
p7 = DimPlot(selvul, reduction = "ref.umap", group.by = "predicted.id") + ggtitle("Leng, 2021") + xlim(-6,6) + ylim(-6,6)
p8 = DimPlot(tsai, reduction = "ref.umap", group.by = "predicted.id") + ggtitle("Blanchard, 2022") + xlim(-6,6) + ylim(-6,6)
p9 = DimPlot(ento, reduction = "ref.umap", group.by = "predicted.id") + ggtitle("Grubman, 2019") + xlim(-6,6) + ylim(-6,6)
p10 = DimPlot(dlpfc2, reduction = "ref.umap", group.by = "predicted.id") + ggtitle("ROSMAP DLPFC-2") + xlim(-6,6) + ylim(-6,6)
p11 = DimPlot(pu1, reduction = "ref.umap", group.by = "predicted.id", raster = F) + ggtitle("Prater, 2023") + xlim(-6,6) + ylim(-6,6)

## Fig. 4A ----
(wrap_plots(list(p1,p11,p4,p9,p7,p2,p6,p3,p5,p8,p10)) + plot_layout(guides = "collect")) %>% 
  ggsave("projections.png",.,width=12,height=10,dpi=600)

# specific metadata plots for supplement
# from Lau, 2020
## Fig. 4C ----
setkey(gse157827_md, "predicted.id","orig.ident") 
dat = gse157827_md[CJ(predicted.id, orig.ident, unique = TRUE),.N,by=.EACHI]
dat[, frac := N/sum(N) * 100, by = orig.ident]
dat = dat[predicted.id == "TIMs" & frac < 7 & frac > 0] #filter out samples with too-high TIM fraction, as these are likely stress induced
dat$disease = substr(dat$orig.ident, 1, 2)
pval = compare_means(frac ~ disease, data = dat, method = "t.test")$p.adj
dat$disease = ifelse(dat$disease == "AD", "Alzheimer's","Normal Control")
dat$disease = factor(dat$disease, levels = c("Normal Control","Alzheimer's"))
(ggplot(dat, aes(x = disease, y = frac)) +
    geom_boxplot(aes(fill = disease), outlier.shape = NA) + 
    geom_point(size=0.75) +
    labs(x = "", y = "Frequency of TIMs (%)") + 
    theme_bw() + 
    ylim(NA,7.5) +
    theme(legend.position = "none") + 
    geom_bracket(xmin = "Alzheimer's", xmax = "Normal Control", y.position = 7, label = paste0("p = ", pval)) + 
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
dat = dat[predicted.id == "TIMs" & frac>0 & PMI != "Reference"]
dat = dat %>%
  mutate(UQ(rlang::sym(col)) := factor(UQ(rlang::sym(col)), levels = levels(sea_mouse[[col]][,1])))
counts = dat[,.N,by=PMI]
(ggplot() + 
    geom_boxplot(data = dat, aes(x = UQ(rlang::sym(col)), y = frac, fill = UQ(rlang::sym(col))),color = "black", outlier.shape = NA) +
    geom_text(data = counts, aes(x=PMI,y=0.1,label = paste0("N = ",N))) + 
    geom_point(data = dat, aes(x = UQ(rlang::sym(col)), y = frac), size=0.75) +
    stat_compare_means(data = dat,aes(x = UQ(rlang::sym(col)), y = frac),method="kruskal.test", label.x = 0.8, label.y = 5.9) +
    theme_bw() + 
    theme(axis.title.x = element_blank()) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(stringr::str_replace(x,"-"," "), width = 10)) +
    labs(y = "Frequency of TIMs (%)", title = "SEA-AD: TIM Frequency by PMI") + 
    coord_cartesian(ylim = c(0,6)) +
    theme(legend.position = "none")) %>% 
  ggsave("SEA_tims_by_pmi.png",.,height = 4, width = 4, dpi = 600)

# from Leng, 2021
## Fig. 4D ----
setkey(selvul_md, "predicted.id","SampleID") 
dat = selvul_md[CJ(predicted.id, SampleID, unique = TRUE),.N,by=.EACHI]
dat[, frac := N/sum(N) * 100, by = SampleID]
dat = dat[predicted.id == "TIMs" & frac < 7 & frac > 0] # filter out samples with too-high TIM fraction, as these are likely stress induced
braaks = selvul_md[,c("BraakStage","SampleID")] %>% unique
dat = merge.data.table(dat, braaks, by = "SampleID")
dat$disease = ifelse(dat$BraakStage > 0, "High","Low")
ri.stat = t.test(dat[disease == "High",frac],
                 mu=mean(dat[disease == "Low",frac]),
                 alternative="two.sided")$p.value %>% round(digits=4)

dat$disease = factor(dat$disease, levels = c("Low","High"))
(ggplot(dat, aes(x = disease, y = frac)) +
    geom_boxplot(aes(fill = disease),color="black") + 
    geom_point() + 
    labs(x = "", y = "Frequency of TIMs (%)",fill="Braak\nStage") + 
    theme_bw() + 
    coord_cartesian(ylim=c(NA,8)) +
    geom_bracket(xmin = 1, xmax = 2, y.position = 7.5, label = paste0("p = ", ri.stat), vjust = -0.3) + 
    ggtitle("Leng, 2021: TIM Frequency\nby Braak Stage") + 
    theme(legend.position = "none")) %>% 
  ggsave("leng_2021_tims_by_braak_alt.png",.,height = 4, width = 4, dpi = 600)

# from blanchard, 2022
## Fig. 4E ----
tsai_summary = tsai_md[,.N,by=.(pmi,amyloid_hi,APOE4,nft)]
tsai_summary$num_tim = tsai_md[,sum(is_tim),by=.(pmi,amyloid_hi,APOE4,nft)]$V1
tsai_summary$frac = tsai_summary$num_tim/tsai_summary$N
tsai_summary_filt = tsai_summary[N > 50 & nft > 0]

pval = compare_means(frac ~ APOE4, data = tsai_summary_filt, method = "t.test")$p.adj

(ggplot(tsai_summary_filt, aes(x=APOE4,y=100*frac)) + 
    geom_boxplot(aes(fill=APOE4),color = "black",outlier.shape = NA) + 
    geom_point(size=0.75) + 
    theme_bw() +
    labs(x = "", y = "Frequency of TIMs (%)") + 
    geom_bracket(xmin = "no", xmax = "yes", y.position = 3.9, label = paste0("p = ", pval), tip.length = 0.02) + 
    theme(legend.position = "none") + 
    coord_cartesian(ylim = c(NA,4)) +
    scale_x_discrete(labels = c("APOE3/APOE3","APOE4 Carrier")) + 
    ggtitle("Blanchard, 2022: TIM\nFrequency by APOE Isoform")) %>%
  ggsave("blanchard_with_ttest.png",.,height=4,width=4,dpi=600)

# from ROSMAP DPLFC-2
## Fig. 4G ----
rosmap = fread("DLPFC_2/ROSMAP_clinical.csv")
setkey(dlpfc2_md, "predicted.id","individualID") 
dat = dlpfc2_md[CJ(predicted.id, individualID, unique = TRUE),.N,by=.EACHI]
dat[, frac := N/sum(N) * 100, by = individualID]
dat[, total := sum(N), by = individualID]
dat = dat[!is.na(individualID) & predicted.id == "TIMs"]
dat = merge.data.table(dat, rosmap, by = "individualID")
dat$e4_count = lengths(regmatches(dat$apoe_genotype, gregexpr("4", dat$apoe_genotype))) %>% factor(ordered = T)
dat$sex = ifelse(dat$msex == 1, "Male","Female")
dat[, age_death_numeric := ifelse(age_death == "90+", 90, as.numeric(age_death))]

lmsum = lm(frac ~ 1 + e4_count + braaksc + ceradsc + cogdx + sex + pmi + educ + age_death_numeric + e4_count*sex + e4_count*braaksc, 
           data = dat %>% mutate(e4_count = as.numeric(e4_count))) %>% summary
lmsum_dat = data.table(lmsum$coefficients,keep.rownames = TRUE)[-c(1,10:11),]
lmsum_dat$rn = c("E4 Allele Count","Braak Score","CERAD Score","Cognitive Diagnosis","Male Sex","Postmortem Interval","Education","Age at Death")
lmsum_dat$rn = factor(lmsum_dat$rn) %>% forcats::fct_reorder(lmsum_dat$Estimate)
lmsum_dat$lower = lmsum_dat$Estimate - lmsum_dat$`Std. Error`
lmsum_dat$upper = lmsum_dat$Estimate + lmsum_dat$`Std. Error`
lmsum_dat$se = paste0("\u00B1",formatC(lmsum_dat$`Std. Error`,format="f",digits=2))
lmsum_dat$est_alt = formatC(lmsum_dat$Estimate,format="f",digits=2)
lmsum_dat$label_alt = ifelse(lmsum_dat$`Pr(>|t|)` < 0.05,
                             formatC(lmsum_dat$`Pr(>|t|)`, format = "e", digits = 2),
                             formatC(lmsum_dat$`Pr(>|t|)`, format = "f", digits = 2))
lmsum_dat = lmsum_dat[order(-rn)]

png("forest_plot_se.png",height=4.5,width=9.5,units="in",res = 600)
lmsum_dat %>% 
  forestplot(labeltext = c("rn","est_alt","se","label_alt"), 
             mean = "Estimate", graph.pos = 2, boxsize = 0.2, vertices = T,
             title = "ROSMAP DLPFC-2: Multiple Linear Regression of TIM Frequency") %>% 
  fp_add_lines(h_2 = gpar(col = "black"),
               h_5 = gpar(col = "black", lty = 2)) %>%
  fp_set_style(box = "royalblue4", line = "black",align = "lccc",
               txt_gp = fpTxtGp(ticks = gpar(cex = 1.1))) %>% 
  fp_add_header("rn" = c("") %>% fp_align_center(),
                "est_alt"=c("Estimate") %>% fp_align_center(),
                "se"=c("Std. Error") %>% fp_align_center(),
                "label_alt"=c("P-Value") %>% fp_align_center()) %>% 
  fp_set_zebra_style("#EFEFEF")
dev.off()

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

## Fig. S4 ----
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

