
.libPaths()
.libPaths(c('~/R/x86_64-pc-linux-gnu-library/4.4',
            '/refdir/Rlib',
            '/usr/local/lib/R/library'))

################################################################################
# scenic_R
rm(list=ls())
gc()
library(Seurat)
load("./rdata/scRNA_cellcall.rda")
table(scRNA_cellcall@meta.data$group)
scRNAsub = subset(scRNA_cellcall,group %in% c("Gpnmb-Mac",    "Gpnmb+Mac"))
set.seed(1314)
# setwd('./myelo_pseudotime/')
table(scRNAsub@meta.data$group)
# Gpnmb-Mac Gpnmb+Mac 
# 774       149 
table(scRNAsub@meta.data$condition)
# D   E   M   P 
# 267 266 169 221 

save(scRNAsub,file = "./rdata/scenic_mac.rda")
rm(list=ls())
gc()

# 有错记得本地安装

BiocManager::install("RcisTarget")
# devtools::install_local('./RcisTarget-master.zip')
devtools::install_github("aertslab/SCENIC")
# .libPaths(new = "~/R/R4_3/")
# .libPaths(new = "~/R/x86_64-pc-linux-gnu-library/4.2/")
# .libPaths()
install.packages(pkgs = "RSQLite", 
                 dependencies = c("Depends", "Imports"))
# To download the rankings:


# Sys.setenv(LIBARROW_MINIMAL = "false") (for all optional features, including 'zstd')
Sys.setenv(ARROW_WITH_ZSTD = "ON") # (for just 'zstd')
BiocManager::install(c("purrr"))
install.packages(c("assertthat", "bit64", "purrr", "tidyselect", "cpp11"))
install.packages("arrow")
library(RSQLite)
library(SCENIC)
library(Seurat)
library(RcisTarget)
## 内存不够控制细胞数目！！！不要超过1000细胞！！！！！！！
set.seed(123456)
load("./rdata/scenic_mac.rda")
#a=sample(1:ncol(scRNAsub),200)
#scRNAsub@reductions$nmf=NULL
#scRNAsub=scRNAsub[,a]

exprMat <- as.matrix(scRNAsub@assays$RNA@counts)
cellInfo <-scRNAsub@meta.data

################################################################################
# CELLCYCLE
# 进入子文件夹
setwd("./SCENIC_CYCLE") 
load("../rdata/scenic_mac.rda")
Idents(scRNAsub)=scRNAsub$condition

### Initialize settings
dir.create("int")
saveRDS(cellInfo, file="int/cellInfo.Rds")

Idents(scRNAsub) <- "condition"

# 人物种
org <- "mgi" # or hgnc, or dmel
dbDir="cisTarget_databases" # RcisTarget databases location, 我们要把数据自行防止
# 或者create a link in the current directory to avoid duplicating them
# like system("ln -s ~/path/to/dbs databases")
myDatasetTitle="SCENIC Challenge" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
dbs
# data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
# motifAnnotations_hgnc <- motifAnnotations_hgnc_v9
# motifAnnotations_mgi
data(list="motifAnnotations_mgi_v9", package="RcisTarget")
motifAnnotations_mgi <- motifAnnotations_mgi_v9
scenicOptions <- initializeScenic(org=org,
                                  dbDir=dbDir, 
                                  dbs = dbs,
                                  datasetTitle=myDatasetTitle, 
                                  nCores=32)

scenicOptions@inputDatasetInfo$cellInfo <- "cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "colVars.Rds"

# 节省内存，我们只取一个基因组文件
# 人：用下面！！！！！！！！！
# scenicOptions@settings$dbs <- c("mm9-5kb-mc8nr" = "hg19-tss-centered-10kb-7species.mc9nr.feather")
# scenicOptions@settings$db_mcVersion <- "v8"
# saveRDS(scenicOptions, file="scenicOptions.Rds") 

#Gene filter: Filter by the number of cells in which the gene is detected 
# (minCountsPerGene, by default 6 UMI counts across all samples)
# and by the number of cells in which the gene is detected 
# (minSamples, by default 1% of the cells)
# 尽可能减少基因数目
# 如果需要更多基因，0.1调成0.05或0.01
genesKept <- geneFiltering(exprMat, 
                           scenicOptions=scenicOptions,
                           minCountsPerGene=3*.1*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept,]
dim(exprMat_filtered)
# [1] 1753  923
#calculating correlation
runCorrelation(exprMat_filtered, scenicOptions)

#GENIE3: infer potential transcription factor targets based on the expression data
exprMat_filtered <- log2(exprMat_filtered+1) 


# 等待较久的一步，纳入基因和细胞越多，等待越久，此处只用了1000多基因，200个细胞
# 对应的TF也只有100多个了，其实有1000多个TF，
# 如果是性能好的计算机可以考虑前面0.1调成0.05或0.01
runGenie3(exprMat_filtered, scenicOptions)


scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 16
scenicOptions@settings$seed <- 123


exprMat_log <- log2(exprMat+1)
dim(exprMat)
# [1] 23657   923
#scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # For toy run
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)

# step2容易崩，节省内存只看top50
gc()
# scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top50")) 
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) 

gc()
# 回到单线程！否则容易报错！！！！！
scenicOptions@settings$nCores <- 1
library(foreach)

scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
saveRDS(scenicOptions, file="scenicOptions.Rds") # To save status

#============================================
scenicOptions <- readRDS("scenicOptions.Rds")
#下次可以直接读取他,我们已经越过最耗时间的步骤
# 如果经常报错，不要重新readRDS和load，应不间断运行

#scenicOptions <- readRDS("scenicOptions.Rds")
#load('scRNAsub.Rdata')
scenicOptions@settings$seed <- 123 # same seed for all of them

#Binarizing the network
runSCENIC_4_aucell_binarize(scenicOptions)

aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
aucell_regulonAUC
# AUC for 76 gene sets (rows) and 923 cells (columns).
aucell_regulonAUC.t <- t(aucell_regulonAUC@assays@data$AUC)
# colnames(aucell_regulonAUC.t) <- gsub(" ", "_", colnames(aucell_regulonAUC.t))
# rownames(aucell_regulonAUC.t) <- gsub("[.]", "-", rownames(aucell_regulonAUC.t))


load("../rdata/scenic_mac.rda")
set.seed(1314)
# setwd('./myelo_pseudotime/')
# scRNAsub <- subset(scRNA_N, downsample=1000)
table(scRNAsub@meta.data$condition)
# we need to reorder our group of condition
scRNAsub@meta.data$condition = factor(scRNAsub@meta.data$condition,levels = c("P","E","M","D"))
head(scRNAsub@meta.data)
mac.scenic <- scRNAsub
# rownames(fibro.scenic@meta.data) <-  gsub("[.]", "-", rownames(fibro.scenic@meta.data) )
# colnames(fibro.scenic@meta.data) <- gsub(" ", "_", fibro.scenic@meta.data)

# 结果导入endo.scenic
# please make sure the format of rownames
mac.scenic@meta.data <- cbind(mac.scenic@meta.data, 
                              aucell_regulonAUC.t[rownames(mac.scenic@meta.data),])
dev.off()

# set.seed(1314)
# mac.scenic <- FindVariableFeatures(mac.scenic, selection.method = "vst", nfeatures = 2000) 
# mac.scenic<- ScaleData(mac.scenic)
# mac.scenic<- RunPCA(mac.scenic) 
# mac.scenic <- FindNeighbors(mac.scenic, dims =1:30) 
# mac.scenic <- FindClusters(mac.scenic,resolution = 0.4)
# mac.scenic <- RunUMAP(mac.scenic)
# DimPlot(mac.scenic, reduction = "umap")
# 找6个
# FeaturePlot(mac.scenic, reduction = "umap", 
            # features = colnames(mac.scenic@meta.data)[c(6,9,13,14)], 
            # cols = c("yellow", "red"))
# FeaturePlot(scRNA_N,features = c("BRAC1","EZH2"))
#Regulon scores heatmap 热图
Idents(mac.scenic)=mac.scenic$condition
cells.ord.cluster <- mac.scenic@active.ident
cells.ord.cluster<- cells.ord.cluster[order(cells.ord.cluster)]
regulon.scores <- t(aucell_regulonAUC.t[names(cells.ord.cluster),])
regulon.scores.log <- log(regulon.scores +1)
regulon.scores.scaled <- scale(regulon.scores)

library(gplots)
library(pheatmap)

# 进一步标化
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
data_subset_norm <- t(apply(regulon.scores, 1, cal_z_score))
#pheatmap(data_subset_norm)

cluster.col <- data.frame(mac.scenic@active.ident, 
                          row.names = names(mac.scenic@active.ident))
colnames(cluster.col) <- "condition"


pheatmap::pheatmap(regulon.scores.scaled, cluster_rows = T,cluster_cols = F, annotation_col = cluster.col, show_colnames = F, fontsize_row=5)
#pheatmap::pheatmap(data_subset_norm, cluster_rows = T,cluster_cols = F, annotation_col = cluster.col, show_colnames = F)


#Average Regulon Activity 平均调控活性
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
cellInfo <-scRNAsub@meta.data
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$condition),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))

regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
library(ggplot2)
library(cowplot)
# "#377EB8","#F0F0F0","#E41A1C"
pheatmap::pheatmap(regulonActivity_byCellType_Scaled,
                   cluster_cols = F,
                   cluster_rows = F,
                   color = colorRampPalette(c(rep("#377EB8",1), "white", rep("#E41A1C",1)))(100))

dev.off()
pheatmap::pheatmap(regulonActivity_byCellType_Scaled,
                   cluster_cols = F,
                   cluster_rows = F,
                   color = colorRampPalette(c(rep("#377EB8",1), "white", rep("firebrick3",1)))(100))



auc <- aucell_regulonAUC@assays@data$AUC
auc[1:5,1:2]
load("../rdata/scenic_mac.rda")
scRNAsub@meta.data$condition = factor(scRNAsub@meta.data$condition,levels = c("P","E","M","D"))
head(scRNAsub@meta.data)
set.seed(1314)
# setwd('./myelo_pseudotime/')
# scRNAsub <- subset(scRNA_N, downsample=1000)
# table(scRNAsub@meta.data$celltype)
# head(scRNA_sub@meta.data)
# 计算Regulon Specificity Score (RSS)评估细胞特异性 # 计算RSS
rss <- calcRSS(AUC=auc, cellAnnotation=scRNAsub@meta.data$condition)
rss[1:5,]
rssPlot <- plotRSS(rss)
## RSS > 0.01被用于可视化
plotly::ggplotly(rssPlot$plot)
plotRSS_oneSet(rss,setName = 'E',n=10)




################################################################################
# GPNMB
rm(list = ls())
scenicOptions <- readRDS("scenicOptions.Rds")
#下次可以直接读取他,我们已经越过最耗时间的步骤
# 如果经常报错，不要重新readRDS和load，应不间断运行

#scenicOptions <- readRDS("scenicOptions.Rds")
#load('scRNAsub.Rdata')
scenicOptions@settings$seed <- 123 # same seed for all of them

#Binarizing the network
runSCENIC_4_aucell_binarize(scenicOptions)

aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
aucell_regulonAUC
aucell_regulonAUC.t <- t(aucell_regulonAUC@assays@data$AUC)

load("../rdata/scenic_mac.rda")
set.seed(1314)
# setwd('./myelo_pseudotime/')
# scRNAsub <- subset(scRNA_N, downsample=1000)
table(scRNAsub@meta.data$group)
head(scRNAsub@meta.data)
mac.scenic <- scRNAsub
# rownames(fibro.scenic@meta.data) <-  gsub("[.]", "-", rownames(fibro.scenic@meta.data) )
# colnames(fibro.scenic@meta.data) <- gsub(" ", "_", fibro.scenic@meta.data)

# 结果导入endo.scenic
mac.scenic@meta.data <- cbind(mac.scenic@meta.data, 
                              aucell_regulonAUC.t[rownames(mac.scenic@meta.data),])
dev.off()

# set.seed(1314)
# mac.scenic <- FindVariableFeatures(mac.scenic, selection.method = "vst", nfeatures = 2000) 
# mac.scenic<- ScaleData(mac.scenic)
# mac.scenic<- RunPCA(mac.scenic) 
# mac.scenic <- FindNeighbors(mac.scenic, dims =1:30) 
# mac.scenic <- FindClusters(mac.scenic,resolution = 0.4)
# mac.scenic <- RunUMAP(mac.scenic)
# DimPlot(mac.scenic, reduction = "umap")
# 找6个
# FeaturePlot(mac.scenic, reduction = "umap", 
# features = colnames(mac.scenic@meta.data)[c(6,9,13,14)], 
# cols = c("yellow", "red"))
# FeaturePlot(scRNA_N,features = c("BRAC1","EZH2"))
#Regulon scores heatmap 热图
Idents(mac.scenic)=mac.scenic$group
cells.ord.cluster <- mac.scenic@active.ident
cells.ord.cluster<- cells.ord.cluster[order(cells.ord.cluster)]
regulon.scores <- t(aucell_regulonAUC.t[names(cells.ord.cluster),])
regulon.scores.log <- log(regulon.scores +1)
regulon.scores.scaled <- scale(regulon.scores)

library(gplots)
library(pheatmap)

# 进一步标化
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
data_subset_norm <- t(apply(regulon.scores, 1, cal_z_score))
#pheatmap(data_subset_norm)

cluster.col <- data.frame(mac.scenic@active.ident, 
                          row.names = names(mac.scenic@active.ident))
colnames(cluster.col) <- "group"


pheatmap::pheatmap(regulon.scores.scaled, cluster_rows = T,cluster_cols = F, annotation_col = cluster.col, show_colnames = F, fontsize_row=5)
#pheatmap::pheatmap(data_subset_norm, cluster_rows = T,cluster_cols = F, annotation_col = cluster.col, show_colnames = F)


#Average Regulon Activity 平均调控活性
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
cellInfo <-scRNAsub@meta.data
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$group),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
library(ggplot2)
library(cowplot)
# "#377EB8","#F0F0F0","#E41A1C"
pheatmap::pheatmap(regulonActivity_byCellType_Scaled,
                   cluster_cols = F,
                   cluster_rows = F,
                   color = colorRampPalette(c(rep("#377EB8",1), "white", rep("firebrick3",1)))(100))

dev.off()


auc <- aucell_regulonAUC@assays@data$AUC
auc[1:5,1:2]
load("../rdata/scenic_mac.rda")
head(scRNAsub@meta.data)
set.seed(1314)
# setwd('./myelo_pseudotime/')
# scRNAsub <- subset(scRNA_N, downsample=1000)
# table(scRNAsub@meta.data$celltype)
# head(scRNA_sub@meta.data)
# 计算Regulon Specificity Score (RSS)评估细胞特异性 # 计算RSS
rss <- calcRSS(AUC=auc, cellAnnotation=scRNAsub@meta.data$group)
rss[1:5,]
rssPlot <- plotRSS(rss)
## RSS > 0.01被用于可视化
plotly::ggplotly(rssPlot$plot)
plotRSS_oneSet(rss,setName = 'Gpnmb+Mac',n=10)


