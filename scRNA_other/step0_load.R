library(Seurat)
rm(list = ls())
gc()
young1 <- Read10X(data.dir = "./GSE255690_RAW/GSM8077816_young_1/")
young1 <- CreateSeuratObject(counts = young1,project = "young1")

young2 <- Read10X(data.dir = "./GSE255690_RAW/GSM8077817_young_3/")
young2 <- CreateSeuratObject(counts = young2,project = "young2")

young3 <- Read10X(data.dir = "./GSE255690_RAW/GSM8077818_young_4/")
young3 <- CreateSeuratObject(counts = young3,project = "young3")

middle1 <- Read10X(data.dir = "./GSE255690_RAW/GSM8077819_middle_2/")
middle1 <- CreateSeuratObject(counts = middle1,project = "middle1")

middle2 <- Read10X(data.dir = "./GSE255690_RAW/GSM8077820_middle_3/")
middle2<- CreateSeuratObject(counts = middle2,project = "middle2")

middle3 <- Read10X(data.dir = "./GSE255690_RAW/GSM8077821_middle_4/")
middle3<- CreateSeuratObject(counts = middle3,project = "middle3")

old1 <- Read10X(data.dir = "./GSE255690_RAW/GSM8077822_old_1/")
old1<- CreateSeuratObject(counts = old1,project = "old1")

old2 <- Read10X(data.dir = "./GSE255690_RAW/GSM8077823_old_2/")
old2<- CreateSeuratObject(counts = old2,project = "old2")

old3 <- Read10X(data.dir = "./GSE255690_RAW/GSM8077824_old_3/")
old3<- CreateSeuratObject(counts = old3,project = "old3")

scRNA_ovary = merge(young1,list(young2,young3,middle1,middle2,middle3,old1,old2,old3))
save(scRNA_ovary,file = "rdata/ovary_raw.rda")



library(Seurat)
rm(list = ls())
gc()

set.seed(1314)
load("rdata/ovary_raw.rda")
################################################################################
# step2 Standard pre-processing workflow
################################################################################

library(ggplot2)
# 2.1. Screen
scRNA = scRNA_ovary
##mitochondrion, ribosome、erythrocyte
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
scRNA[["percent.rb"]] <- PercentageFeatureSet(scRNA, pattern = "^RP[SL]")
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(scRNA))
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes)
### QC violin picture
# possible theme
theme.set2 = theme(axis.title.x=element_blank())
# the elements of picture
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb", "percent.HB")
group = "orig.ident"
# vlnplot before qc
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
library(patchwork)
violin <- wrap_plots(plots = plots, nrow=2)    
ggsave("./preanalysis/vlnplot_before_qc.pdf", plot = violin, width = 12, height = 8)
### set QC strandards
minGene=300
maxGene=8000
maxUMI=50000
pctMT=50
### vlnplot after QC
scRNA <- subset(scRNA, subset =  nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT & nCount_RNA < maxUMI)
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2)     
ggsave("./preanalysis/vlnplot_after_qc.pdf", plot = violin, width = 12, height = 8) 
# scRNA@meta.data$
table(scRNA@meta.data$orig.ident)
# middle1 middle2 middle3    old1    old2    old3  young1  young2  young3 
# 11270   10080   10524   11577    9376    7555    8217   11342   10897 
save(scRNA,file = './rdata/ov_scRNA_qc.rda')


# 2.2 harmony for integrate
rm(list = ls())
gc()
library(harmony)
library(Seurat)
library(tidyverse)
library(patchwork)
load("rdata/ov_scRNA_qc.rda")
# scRNA = subset(scRNA,orig.ident != "T162")
set.seed(1314)
### standardize data with SCT v2!
# scRNA <- SCTransform(scRNA, vst.flavor = "v2")
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA<- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, features = scale.genes)
DefaultAssay(scRNA)

### PCA
scRNA <- RunPCA(scRNA, npcs=30, verbose=FALSE)

scRNA <- RunHarmony(scRNA, 
                    group.by.vars=c("orig.ident"),
                    assay.use="RNA", 
                    max_iter = 20) 
# group.by.vars: The parameter is to set which group to integrate by
# max.iter.harmony: Set the number of iterations, default is 10. When running RunHarmony, results will indicate how many iterations have elapsed before convergence is complete.


# dimensionality reduction --------------------------------------------------
# 一定要指定“harmony”！
scRNA <- FindNeighbors(scRNA, dims = 1:30, reduction = "harmony")
scRNA <- FindClusters(scRNA,resolution = 1)
# remove.packages("Matrix")
# install.packages("Matrix", type = "source"，)
# install.packages("irlba", type = "source")

scRNA <- RunUMAP(scRNA,reduction = 'harmony', dims = 1:30, verbose = FALSE)

# Plot final cluster resolution --------------------------------------------------------------------
mycolor = c("#A0C2E7","#6894B9","#386CB0","#8798A6","#E0D9E0",
            "#EDBAA7","#FADB7F","#F3B646","#EF9749","#B27466",
            "#646F3F","#899678","#C2BC9A","#868A63","#C4C3BE",
            "#DFA0A6","#98B3D9","#E4BE92","#CB6B7A","#D5CBDA",
            "#f1707d", "#f15536","#ef5767","#ae716e","#cb8e85",
            "#cf8878","#c86f67","#f1ccb8","#f2debd","#b8d38f",
            "#ddff95","#7FC97F","#ff9b6a","#f1b8f1","#d9b8f1",
            "#f1ccb8","#f1f1b8","#b8f1ed","#e7dbca","#e26538",
            "#f3d751","#fd803a","#fe997b","#c490a0"
)

final_embedding <- DimPlot(scRNA, group.by = "seurat_clusters",cols = mycolor,label = T,raster=FALSE)

final_embedding

ggsave(final_embedding,filename=paste0("preanalysis",'/umap_firstcluster.pdf'),width = 16,height = 12)
# \\d 表示数字
# + 多次
# $ 结尾
scRNA@meta.data$Type = gsub("\\d+$", "", scRNA@meta.data$orig.ident)
final_embedding <- DimPlot(scRNA, split.by = "Type",cols = mycolor[c(3:length(mycolor))],label = T,raster=FALSE)

final_embedding

ggsave(final_embedding,filename=paste0("preanalysis",'/umap_type.pdf'),width = 16,height = 6)

save(scRNA,file = "./rdata/ov_scRNA_clusters_29.rda")

################################################################################
# step3 Annotation
################################################################################
# Find HVG and cluster
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(stringr)
rm(list=ls())
options(stringsAsFactors = F)
load(file = "./rdata/ov_scRNA_clusters_29.rda")

library(SingleR)
refdata <- SingleR::HumanPrimaryCellAtlasData() # human
# refdata <- SingleR::MouseRNAseqData()
testdata <- GetAssayData(scRNA, slot="data")
clusters <- scRNA@meta.data$seurat_clusters # seurat_cluster
cellpred <- SingleR(test = testdata, ref = refdata,
                    labels =refdata$label.main,
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")

celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}

Idents(scRNA)=scRNA$celltype
scRNA@meta.data$Type = gsub("\\d+$", "", scRNA@meta.data$orig.ident)
table(scRNA@meta.data$Type)
# middle    old  young 
# 31874  28508  30456
p = DimPlot(scRNA, group.by="celltype", label=F, label.size=4.5, reduction='umap',cols = mycolor)
ggsave(p,filename=paste0("preanalysis",'/umap_singleR.pdf'),width = 16,height = 12)
save(scRNA,file = "./rdata/scRNA_singleR.rda")

FeaturePlot(scRNA,features = c("GPNMB"))
library(Seurat)
load("rdata/scRNA_singleR.rda")
levels(scRNA)


scRNA_N = subset(scRNA, idents=c('Neurons'))
# 一定要指定“harmony”！
scRNA_N <- FindNeighbors(scRNA_N, dims = 1:30, reduction = "harmony")
scRNA_N <- FindClusters(scRNA_N,resolution = 0.3)
scRNA_N <- RunUMAP(scRNA_N,reduction = 'harmony', dims = 1:30, verbose = FALSE)

################################################################################
# manually
rm(list=ls())
gc()
options(stringsAsFactors = F)
load(file = "./rdata/ov_scRNA_clusters_29.rda")
#################################################################################################
# cell annotation
# Specify genes  - Mac
library(ggplot2)
genes_to_check = c("ZP3","TUBB8", # Oocyte
                   "AMH","HSD17B1","SERPINE2","GSTA1", # Granulosa
                   "DCN","LUM","STAR", # theca & stroma
                   "TAGLN","RGS5","ACTA2","MUSTN1", # smooth muscle
                   "VWF","CLDN5","TM4SF1", # Endothelium
                   "PTPRC", # immune
                   "CD53","CXCR4","NKG7","IL7R",
                   "TYROBP","IFI30","CD14","CD68")
# featureplot
p <- FeaturePlot(scRNA, features = genes_to_check)
ggsave(p,filename=paste0("preanalysis",'/fearureplot_sepcify.png'),width = 16,height = 10)
# All on Dotplot 
p <- DotPlot(scRNA, features = genes_to_check) + coord_flip()
ggsave(p,filename=paste0("preanalysis",'/dotplot_sepcify_specific.pdf'),width = 16,height = 12)


#  Need to look at the picture, determine the cell subsets:
celltype=data.frame(ClusterID=0:29,
                    celltype='unkown')
celltype[celltype$ClusterID %in% c(27),2]='Oocyte' 
celltype[celltype$ClusterID %in% c(20,24),2]='Granulosa'
celltype[celltype$ClusterID %in% c(0,1,2,3,6,16,18,21),2]='Theca&stroma'
celltype[celltype$ClusterID %in% c(9,19,10,4,7,14,29),2]='SMC'

celltype[celltype$ClusterID %in% c(5,8,13,17),2]='Endo'
celltype[celltype$ClusterID %in% c(15,22),2]='Myeloid'
celltype[celltype$ClusterID %in% c(11,12,23,28),2]='Lymphocyte'

head(celltype)
celltype 
table(celltype$celltype)

new.cluster.ids <- celltype$celltype
names(new.cluster.ids) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, new.cluster.ids)
scRNA@meta.data$celltype = scRNA@active.ident
scRNA <-subset(scRNA,celltype != "unkown")
save(scRNA,file = 'rdata/scRNA_after_annotation.Rdata')

# nolegend
source("mycolor.R")
p<-DimPlot(scRNA, reduction = "umap", label = TRUE, pt.size = 0.5,cols = mycolor) + NoLegend()
ggsave(p,filename=paste0("preanalysis",'/umap_first_cluster.pdf'),width = 8,height = 6)

table(scRNA@meta.data$celltype)
# Theca&stroma          SMC         Endo   Lymphocyte      Myeloid    Granulosa       unkown       Oocyte 
# 51045        19727        11358         4817         1891         1227            0          185 
# legend
p<-DimPlot(ovary, reduction = "umap", label = TRUE, pt.size = 0.5)
ggsave(p,filename=paste0("preanalysis",'/umap_firstcluster_legend.pdf'),width = 16,height = 12)


FeaturePlot(scRNA,features = c("GPNMB"))
VlnPlot(scRNA,features = c("GPNMB"))
DotPlot(scRNA,features = c("GPNMB"))






