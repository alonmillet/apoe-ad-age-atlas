library(data.table)
library(dplyr)
library(EnsDb.Mmusculus.v79)
library(BisqueRNA)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(patchwork)

setwd("./salmon") # directory containing salmon quant results

# prep bulk data
dirs = list.files()
for(i in seq_along(dirs)){
  dat = fread(paste0(dirs[i],"/quant.sf"))
  dat$genotype = strsplit(dirs[i],"_") %>% sapply(.,"[[",1)
  dat$samplenum = strsplit(dirs[i],"_") %>% sapply(.,"[[",2)
  if(i == 1){
    merged = dat
  } else {
    merged = rbind(merged,dat)
  }
}

merged$Name = strsplit(merged$Name,".",fixed=T) %>% sapply(.,"[[",1)
merged$samp = paste0(merged$genotype, "-", merged$samplenum)
genenames = ensembldb::select(EnsDb.Mmusculus.v79, 
                              keys=merged$Name, 
                              keytype = "TXNAME", 
                              columns = c("SYMBOL","TXNAME"))

merged = merge.data.table(merged,genenames,by.x="Name",by.y="TXNAME")
agg = merged[,.(tpm = sum(TPM)), by = c("SYMBOL","samp")]
agg = dcast(agg, SYMBOL ~ samp, value.var = "tpm")[-1,]
agg.m = agg[,-1] %>% as.matrix
rownames(agg.m) = agg$SYMBOL
bulk.eset = ExpressionSet(agg.m)

# prep singlecell ref
seu = readRDS("adapoe.rds")
seu = RenameCells(seu, add.cell.id = paste0(Cells(seu) %>% strsplit(.,"_",fixed=T) %>% sapply(.,"[[",1),
                                            "-", 
                                            Cells(seu) %>% strsplit(.,"_",fixed=T) %>% sapply(.,"[[",2)))
sc.eset = SeuratToExpressionSet(seu, delimiter = "_", position = 1, version = "v3")

# run Bisque
res = ReferenceBasedDecomposition(bulk.eset, sc.eset, markers = NULL, use.overlap = F)
props = res$bulk.props %>% as.data.table(keep.rownames = T)
fwrite(props, "bisque_decomposition_results.txt")

props = fread("bisque_decomposition_results.txt")
timprops = props[rn %in% c("TIMs","Serpine1+ TIMs","Effector-hi TIMs"),] %>% as.data.table %>% data.table::transpose()
colnames(timprops) = timprops[1,] %>% unlist
colnames(timprops)[1] = "Effector-lo TIMs" #overwrite
timprops = timprops[-1,] 
timprops = timprops[, lapply(.SD, as.numeric)] * 100
timprops$genotype = c(rep("APOE2",5),rep("APOE3",5),rep("APOE4",5))
timprops = timprops[-c(6,13)] # remove outliers

p1 = ggplot(timprops, aes(x=genotype,y=`Effector-hi TIMs`)) + 
  geom_boxplot(aes(fill=genotype)) +
  geom_point() + 
  theme_bw(base_size = 25) + 
  labs(x = "", y = "Fraction of Cells (%)", title = "Effector-hi TIMs") + 
  stat_compare_means(method = "t.test", 
                     comparisons = list(c("APOE2","APOE3"),c("APOE3","APOE4"),c("APOE2","APOE4")),
                     aes(label = after_stat(p.signif)), size = 10, label.y = c(18,18,22)) + 
  theme(legend.position = "none")

p2 = ggplot(timprops, aes(x=genotype,y=`Effector-lo TIMs`)) + 
  geom_boxplot(aes(fill=genotype)) +
  geom_point() + 
  theme_bw(base_size = 25) + 
  labs(x = "", y = "Fraction of Cells (%)", title = "Effector-lo TIMs") + 
  stat_compare_means(method = "t.test", 
                     comparisons = list(c("APOE2","APOE3"),c("APOE3","APOE4"),c("APOE2","APOE4")),
                     aes(label = after_stat(p.signif)), size = 10, label.y = c(18,18,22)) + 
  theme(legend.position = "none")

((p2/p1) & ylim(0,26)) %>% ggsave("bisque_frequencies.png",.,dpi=600,width=11,height=10)