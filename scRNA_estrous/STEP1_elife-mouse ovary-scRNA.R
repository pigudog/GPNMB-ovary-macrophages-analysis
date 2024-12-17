####################STEP 1_数据预处理##############
##########加载包####################################
{library(openxlsx)
library(Seurat)
library(tidyverse)
library(scRNAstat) 
library(patchwork)
library(ggsci)
library(dplyr)
library(ggplot2)
  library(gplots)
library(clustree)
library(cowplot)
}
Sys.setenv(R_MAX_NUM_DLLS=999)
options(stringsAsFactors = F)

##########################STEP1.1——读取数据####################################
rm(list = ls())
setwd()
getwd()
dir()
dir.create("Monocle")

#因为给的RDS文件就是数据预处理之后的，所以现在我们直接读取提取cycling后的rds。
mus<- readRDS("../SCP1914_seurat_ovary_0.rds") 
#mus是包括了PDNL的，用subset提取cycling的
SC1 <- readRDS("SC1.RDS") #之前提取过保存在SC1.RDS
head(SC1)
#SC1是包含cycling的；SC2是取了immune的，macrophage_whole是只有巨噬细胞的。

#前面已经建立了seurat对象，现在对数据可视化下

{###Seurat的基本流程
VlnPlot(SC1, features = c("nFeature_RNA", "nCount_RNA", "fraction.mito"), group.by = "orig.ident", ncol = 3,pt.size = 0.1)
#这里有点小bug，因为identity画的是所有不同亚群的，想得到预处理之前总的,利用group.by。

#sce = basic_qc(SC1,org='mouse', dir = "./seurat_step1_QC/")
####这一步本来是数据质控了下，但是到后面发现过滤掉了卵母细胞，所以就还是用SC1。
#basic_qc和basic_workflow都不适用。




{mus <- ScaleData(object = mus,
                 vars.to.regress = c('nUMI'),
                 model.use = 'linear',
                 use.umi = FALSE)


mus <- FindVariableFeatures(object = mus,
                            mean.function = ExpMean,
                            dispersion.function = LogVMR,
                            x.low.cutoff = 0.0125,
                            x.high.cutoff = 4,
                            y.cutoff = 0.5)
{all.genes <- rownames(SC1)
SC2 <- ScaleData(SC1, features = all.genes) 
SC1 <- ScaleData(SC1, features = VariableFeatures(SC1))}


DefaultAssay(SC1) <- "RNA"
SC1[["percent.mt"]] <- PercentageFeatureSet(SC1, pattern = "^MT-")##[[]]相当于给meta.data中加了一列；$也是可以新添加一列

SC1 <- RunPCA(SC1, features = VariableFeatures(SC1))}




ElbowPlot(SC1)
ElbowPlot(SC1, ndims = 10, reduction = "pca")
VizDimLoadings(SC1, dims = 1:2, reduction = "pca")
DimPlot(SC1, reduction = "pca")
DimHeatmap(SC1, dims = 1, cells = 300, balanced = TRUE)
DimHeatmap(SC1, dims = 1:10, cells = 500, balanced = TRUE)

#聚类,主成分挑的太多也不好，选一个合适的
##试了很多次，只有dims为1:8，分辨率为0.15的时候，刚好细胞可以分为0-9共10群
sce <- FindNeighbors(SC1, dims = 1:10) 
sce <- FindClusters(sce, resolution = 0.1)
#resolution自己多试一下,clustertree这个包可以帮助确定
head(sce@meta.data)
levels(sce@meta.data$seurat_clusters)
##降维可视化
sce <- RunUMAP(sce, dims = 1:10)
sce <- RunTSNE(sce, dims = 1:10)
p1 <- DimPlot(sce, reduction = "umap")
p2 <- DimPlot(sce, reduction = "tsne")
VlnPlot(sce, features = c("Ddx4","Gdf9","Bmp15","Dazl","Zp3"))#都尝试了之后，发现seurat降维分群的都不行，还得是原始的level分群才对
ggsave(p1,filename = "plot_reduction_umap.pdf",width = 6,height = 8)
p1 + p2
saveRDS(sce, file = "ovary0_reduction.rds")}




###细胞注释####################################
{#####细胞注释
library(scMCA)
library(Seurat)
library(patchwork)
library(openxlsx)

{


# find all markers of cluster 2
cluster2.markers <- FindMarkers(sce, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 2)

write.xlsx(as.data.frame(cluster2.markers), 'cluster2.markers.xlsx', rowNames=TRUE)




# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(sce, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

write.xlsx(as.data.frame(cluster5.markers), 'cluster5.markers_DGE from cluster 0 and 3.xlsx', rowNames=TRUE)


### find markers for every cluster compared to all remaining cells, report only the positive ones
sce.markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
sce.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)

write.xlsx(as.data.frame(sce.markers), 'sce.markers.xlsx', rowNames=TRUE)


gene_list= list(
Immune =c("Lyz2","C1qb","C1qc","C1qa"),
GC=c("Hsd17b1","Grb14","Foxo1","Inhba","Inha"),
Endothelial =c("Pecam1","Cdh5","Cd34","Flt1"),
Oocytes=c("Gdf9","Bmp15","Dazl","Zp3"),
SMA =c("Des","Myh11","Acta2","Tagln"),
Stromal =c("Tcf21","Col1a1","Dcn","Pdgfra"),
Cumulus  =c("Kctd14","Pcsk6","Gatm","Slc18a2"),
Epithelial  =c("Upk1b","Upk3b","Krt7","Krt18"),
Theca =c("Cyp17a1","Fdx1","Star","Cyp11a1"))


p2 = DotPlot(sce, assay = "RNA", features = gene_list, 
             group.by = 'seurat_clusters') + theme(axis.text.x = element_text(angle = 45, 
                                                                              vjust = 0.5, hjust=0.5))
p2    
ggsave(p2,filename = 'markers.pdf',width = 10)}




####可视化
VlnPlot(sce, features = c("Des", "Star"))

# you can plot raw counts as well
VlnPlot(sce, features = c("Adgre1", "Enpp2"), slot = "counts", log = TRUE)

FeaturePlot(sce, features = c("Adgre1", "Enpp2", "Cd68", "Cd14", "Lpar1", "Ptgs2", "Lyz2", "Lpar6", 
                               "Fshr"))
table(sce@meta.data$seurat_clusters)


dir()
###统计不同发情周期
metadata <- SC1@meta.data
length(unique(metadata$condition)) 
table(metadata$condition)

library(dplyr)

Idents(SC1) <- "RNA" #修改活跃的idents，slot默认为标准化后的数据
DimPlot(SC1, reduction = "umap",split.by = 'condition')
#图像美学的调整
library(RColorBrewer)
colorplatter <- brewer.pal(9,"Set1")
colorplatter
display.brewer.pal(9,"Set1")
names(colorplatter) <- unique(SC1$Level0)
colorplatter

DimPlot(SC1, cols=colorplatter)

DimPlot(SC1, cols=colorplatter, pt.size = 0.8, label=TRUE, label.size = 4, label.box = TRUE)



##########细胞周期校正##########
##加载细胞周期marker基因，seurat自带了一个列表
s.genes <- mouse_cell_cycle_genes$s.genes
g2m.genes <- mouse_cell_cycle_genes$g2m.genes

##区别于scran包的另外重要的一点就是Seurat包仅提供了人类细胞有关的cell cycle related gene，没有小鼠的。
{convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}
m.s.genes <- convertHumanGeneList(cc.genes$s.genes)
m.g2m.genes <- convertHumanGeneList(cc.genes$g2m.genes)

}


library(SingleCellExperiment)
table(is.na(rowData(SC1)$SYMBOL))

rownames(SC1) <- uniquifyFeatureNames(rowData(SC1)$ENSEMBL, 
                                           rowData(SC1)$SYMBOL)
#若SYMBOL为NA值，则用对应的ENSEMBL替换
library(Seurat)
SC1 <- as.Seurat(SC1, counts = "counts", data = "logcounts")

#library(readr) #导入my_data.zip中的data1.csv
##df <- read_csv(unzip("mouse_cell_cycle_genes.zip", "mouse_cell_cycle_genes.rds"))


str(mouse_cell_cycle_genes)
#List of 2
# $ s.genes  : chr [1:42] "Mcm4" "Exo1" "Slbp" "Gmnn" ...
# $ g2m.genes: chr [1:52] "Nuf2" "Psrc1" "Ncapd2" "Ccnb2" ...






SC1 <- CellCycleScoring(SC1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(View(SC1))
RidgePlot(SC1, features = c("Slbp", "Top2a", "Mcm4", "Ccnb2"), ncol = 2)

SC1 <- RunPCA(SC1, features = c(s.genes, g2m.genes))
DimPlot(SC1)
FeaturePlot(SC1, features=c("Top2a", "Mcm4"))

##细胞周期影响校正
SC2 <- ScaleData(SC1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(SC1))
SC2 <- RunPCA(SC2, features = c(s.genes, g2m.genes))#需要的时间比较久
DimPlot(SC2)
FeaturePlot(SC2, features=c("Top2a", "Mcm4"))
}




##转录因子分析####################################
{
  suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(AUCell)
  library(Seurat)
  library(SCopeLoomR)
  library(SCENIC)
  library(GENIE3)
  library(RcisTarget)
  library(foreach)
  library(ComplexHeatmap)
  library(ggplot2)
  library(ggrepel)
  library(AnnotationHub)
  library(dplyr)
  })

SC1 <- readRDS("~/mus/mus_ovary/monocle/SC1.RDS")
class(SC1)

exprMat <- GetAssayData(SC1, slot='counts')
cellInfo <- SC1@meta.data

head(cellInfo)
table(cellInfo$Level0)
dim(exprMat)

# set the working directory to a new folder
dir.create("SCENIC_MouseOvary")
setwd("SCENIC_MouseOvary")

# Initialize SCENIC settings
org <- "mgi" # or hgnc, or dmel

myDatasetTitle <- "SCENIC example on Mouse ovary" # choose a name for your analysis
dbs <- c('500bp'='mm9-500bp-upstream-7species.mc9nr.feather',
         '10kb'='mm9-tss-centered-10kb-7species.mc9nr.feather')

###############
scenicOptions <- initializeScenic(org=org, dbDir="cisTarget_databases", datasetTitle=myDatasetTitle, nCores=1) 

# Co-expression network
## Gene filter/selection
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
## Correlation
runCorrelation(exprMat_filtered, scenicOptions)

## GENIE3
exprMat_filtered <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered, scenicOptions)
save.image(file = "my_scenic_workspace.RData")



# Build and score the GRN
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"]
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # ** Only for toy run!!

library(doParallel)

scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, log2(exprMat+1))
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings


saveRDS(scenicOptions, file="int/scenicOptions.Rds") # To save status

}

##################STEP2——免疫细胞亚群细分####################################
SC2= subset(x=SC1, idents =  "Immune") #共1323个免疫细胞
table(SC2$condition)  
#D   E   M   P 
#394 372 272 285 

Idents(SC2) <- 'Level1'
levels(SC2)
table(SC2$Level1)



###进一步提取免疫细胞里面的巨噬细胞
macrophage_whole= subset(x=SC2, idents =  "I_Macrophage") 
dim(macrophage_whole) #共941个巨噬细胞
table(macrophage_whole$condition) 
#D   E   M   P 
#282 266 172 221 



##这里调整下默认分群
Idents(macrophage_whole) <- 'condition'
levels(macrophage_whole)


##################绘制细胞比例图######################

#######堆叠柱状图
{
table(SC2$condition)
prop.table(table(Idents(SC2)))

table(SC2$condition,Idents(SC2))
Cellratio <- prop.table(table(Idents(SC2), SC2$condition), margin = 2)
Cellratio 

Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+
  scale_x_discrete(limits = c("D","M", "E","P"))
}

#######桑基图
{# 加载包
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggalluvial)
#构建数据

Ratio <-data.frame(condition=SC2$condition,seurat_clusters=SC2$Level1)
# 通过聚合计算细胞比例
Ratio <- Ratio  %>%
  group_by(seurat_clusters, condition) %>% # 分组
  summarise(n=n()) %>%
  mutate(relative_freq = n/sum(n))

# 将 condition列转换为因子类型，并指定因子级别的顺序
Ratio$condition <- factor(Ratio$condition , levels = c("D","M", "E", "P"))



#堆叠柱状图
mycolor = c('#efb306',
                     '#7db954',
                     '#852f88',
                     '#4e54ac',
                     '#0f8096',
                     'pink',
                     'green')
names(mycolor) <- c("P", "E", "M", "D")                   
                     
ggplot(Ratio, aes(x =reorder(seurat_clusters, -match(seurat_clusters, c("P", "E", "M", "D"))),
                  y= relative_freq, fill = condition,
                  stratum=condition, alluvium=condition)) +
  geom_col(width = 0.5, color='black')+
  geom_flow(width=0.5,alpha=0.4, knot.pos=0.5)+ # 参数knot.pos设置为0.5使连接为曲线面积，就像常见的桑基图
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  scale_fill_manual(values = mycolor, 
                    breaks = names(mycolor), 
                    labels = names(mycolor)) 
 
}


#######查看Enpp2+细胞比例变化

if ("Enpp2_exp" %in% names(SC1@meta.data)) {
  # 提取Enpp2+的细胞
  Enpp2_cells <- which(SC1@meta.data$Enpp2_exp > threshold) # 设定一个表达量阈值
  # 创建一个新的seurat对象，仅包含Enpp2+的细胞
  SC1_enpp2 <- SC1[Enpp2_cells, ]
} else {
  print("Enpp2_exp is not in the metadata of SC1")
}


# 计算Enpp2的表达量
Enpp2_exp <- AverageExpression(SC1, features = "Enpp2")
VlnPlot(SC1, features = "Enpp2")



# 提取Enpp2+的细胞
enpp2_cells <- which(unlist(Enpp2_exp) > 0.25)# 设定一个表达量阈值

# 创建一个新的seurat对象，仅包含Enpp2+的细胞
SC1_enpp2 <- SC1[enpp2_cells, ]

# 创建一个包含样本标识和Enpp2表达量的数据框
# 重新计算Enpp2表达量
Enpp2_exp <- AverageExpression(SC1_enpp2, features = "Enpp2")

# 提取样本标识和condition信息
id <- Idents(SC1_enpp2)
condition <- SC1_enpp2@meta.data$condition

# 创建包含样本标识、condition和Enpp2表达量的数据框
data <- data.frame(Sample = colnames(SC1_enpp2@meta.data), 
                   Condition = condition[id], 
                   Enpp2_exp = Enpp2_exp)
# 创建一个包含Enpp2+细胞的样本标识和表达量的数据框
enpp2_data <- data[data$Sample %in% enpp2_cells, ]

# 绘制Enpp2+细胞在不同样本中的表达变化趋势
library(ggplot2)
ggplot(enpp2_data, aes(x = Sample, y = Enpp2_exp)) + 
  geom_point(size = 3, alpha = 0.8, color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_classic() +
  labs(x = "Sample", y = "Enpp2 expression") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

################## STEP 3差异分析和功能分析############

{
 ####绘制Heatmap of the top 10 markers of each cluster by fold change.
macrophage_whole.markers <-FindAllMarkers(macrophage_whole,only.pos = T,min.pct = 0.1,logfc.threshold = 0.25) 
top10 <- macrophage_whole.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(macrophage_whole, features = top10$gene) + NoLegend()
  
#ggsave(p,filename = "TOP10 genes_heatmap.pdf",width = 6,height = 10)
  


#1、加载包
library(AnnotationDbi)
library(org.Mm.eg.db)#基因注释包
library(clusterProfiler)#富集包
library(dplyr)
library(ggplot2)#画图包

####2、基因id转换，kegg和go富集用的ID类型是ENTREZID）
  macrophage_whole.markers <- read_csv("macrophage_whole.markers.csv")
  dir()
  
  dim(clusterP.markers)
  library(org.Mm.eg.db)
  
  
  ids<-bitr(geneID=macrophage_whole.markers$gene,fromType='SYMBOL',toType='ENTREZID',OrgDb = 'org.Mm.eg.db')
  macrophage_whole.markers=merge(macrophage_whole.markers,ids,by.x='gene',by.y='SYMBOL')
  
  ## 函数split()可以按照分组因子，把向量，矩阵和数据框进行适当的分组。
  ## 它的返回值是一个列表，代表分组变量每个水平的观测。
  gcSample=split(macrophage_whole.markers$ENTREZID, macrophage_whole.markers$cluster) 
  ## KEGG
  xx <- compareCluster(gcSample,
                       fun = "enrichKEGG",
                       organism = "mouse", pvalueCutoff = 0.05
  )
  p <- dotplot(xx)
  p + theme(axis.text.x = element_text(
    angle = 45,
    vjust = 0.5, hjust = 0.5
  ))
xx@compareClusterResult$Description <- gsub("\\s*\\-\\s*Mus musculus.*$", "", xx@compareClusterResult$Description)
xx <- as.data.frame(xx)
write.table(xx, file = "compareClusterResult_go_CC.txt", sep = "\t", quote = FALSE, row.names = TRUE)



## GO
xx <- compareCluster(gcSample,
                     fun = "enrichGO",
                     OrgDb = "org.Mm.eg.db",
                     ont = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05
)
p <- dotplot(xx)
p + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))

##寻找两群之间的差异基因
diff.markers <- FindMarkers(macrophage_whole, ident.1 = "E", ident.2 = "P",test.use = 'MAST')
write.table(diff.markers, file = "E vs P diffGenes.txt", sep = "\t", quote = FALSE, row.names = TRUE)



## 获取上下调基因
gene_up=rownames(diff.markers[diff.markers$avg_log2FC > 0,])
gene_down=rownames(diff.markers[diff.markers$avg_log2FC < 0,])
## 把SYMBOL改为ENTREZID
library(org.Mm.eg.db)
gene_up=as.character(na.omit(AnnotationDbi::select(org.Mm.eg.db,
                                                   keys = gene_up,
                                                   columns = 'ENTREZID',
                                                   keytype = 'SYMBOL')[,2]))
gene_down=as.character(na.omit(AnnotationDbi::select(org.Mm.eg.db,
                                                     keys = gene_down,
                                                     columns = 'ENTREZID',
                                                     keytype = 'SYMBOL')[,2]))
library(clusterProfiler)
## 上调基因
## KEGG
gene_up <- unique(gene_up)
kk.up <- enrichKEGG(gene = gene_up,
                    organism = "mouse",
                    pvalueCutoff = 0.9,
                    qvalueCutoff = 0.9)

kk.up@result$Description <- gsub("\\s*\\-\\s*Mus musculus.*$", "", kk.up@result$Description)

dotplot(kk.up)
write.table(kk.up, file = "E vs P diffGenes_kk.up.txt", sep = "\t", quote = FALSE, row.names = TRUE)


## GO
go.up <- enrichGO(gene = gene_up,
                  OrgDb = org.Mm.eg.db,
                  ont = "CC" ,
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.99,
                  qvalueCutoff = 0.99,
                  readabl = TRUE)
dotplot(go.up)

write.table(go.up, file = "E vs P diffGenes_go_CC.up.txt", sep = "\t", quote = FALSE, row.names = TRUE)


## 下调基因
gene_down <- unique(gene_down)
kk.down <- enrichKEGG(gene = gene_down,
                    organism = "mouse",
                    pvalueCutoff = 0.9,
                    qvalueCutoff = 0.9)

kk.down@result$Description <- gsub("\\s*\\-\\s*Mus musculus.*$", "", kk.down@result$Description)

dotplot(kk.down)
write.table(kk.down, file = "E vs P diffGenes_kk.down.txt", sep = "\t", quote = FALSE, row.names = TRUE)

## GO
go.down <- enrichGO(gene = gene_down,
                  OrgDb = org.Mm.eg.db,
                  ont = "CC" ,
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.99,
                  qvalueCutoff = 0.99,
                  readabl = TRUE)
dotplot(go.down)

write.table(go.down, file = "E vs P diffGenes_go_CC_down.txt", sep = "\t", quote = FALSE, row.names = TRUE)

## 上一步差异分析得到差异基因列表deg后取出，p值和log2FC
nrDEG = diff.markers[,c('avg_log2FC', 'p_val')]
colnames(nrDEG)=c('log2FoldChange','pvalue') ##更改列名
library(org.Mm.eg.db)
library(clusterProfiler)
## 把SYMBOL转换为ENTREZID，可能有部分丢失
gene <- bitr(rownames(nrDEG),     
             fromType = "SYMBOL",     
             toType =  "ENTREZID",    
             OrgDb = org.Mm.eg.db)
## 基因名、ENTREZID、logFC一一对应起来
gene$logFC <- nrDEG$log2FoldChange[match(gene$SYMBOL,rownames(nrDEG))]
## 构建genelist
geneList=gene$logFC
names(geneList)=gene$ENTREZID 
geneList=sort(geneList,decreasing = T) # 降序，按照logFC的值来排序
## GSEA分析
kk_gse <- gseKEGG(geneList     = geneList,
                  organism     = 'mouse',
                  nPerm        = 1000,
                  minGSSize    = 10,
                  pvalueCutoff = 0.9,
                  verbose      = FALSE)
kk_gse=DOSE::setReadable(kk_gse, OrgDb='org.Mm.eg.db',keyType='ENTREZID')
sortkk<-kk_gse[order(kk_gse$enrichmentScore, decreasing = T),]
sortkk$Description <- gsub("\\s*\\-\\s*Mus musculus.*$", "", sortkk$Description)



library(enrichplot)
gseaplot2(kk_gse, 
          "mmu04350", 
          color = "firebrick",
          base_size = 11,
          title = "TGF-beta signaling pathway",
          rel_heights = c(1, .2, .6),
          pvalue_table = TRUE,
          ES_geom = "line")

kk_gse2 <- as.data.frame(kk_gse)
write.table(kk_gse2, file = "E vs P diffGenes_kk_gse.txt", sep = "\t", quote = FALSE, row.names = TRUE)

gseaplot2(kk_gse,row.names(sortkk)[1:5])  



#各种ID之间可以互相转换,toType可以是一个字符串，也可以是一个向量，看自己需求                     
            
                
result_go_clusterM <- enrichGO(clusterM.markers$ENTREZID,OrgDb = "org.Mm.eg.db",
                               keyType = "ENTREZID",ont = "ALL",readable = T,
                               pvalueCutoff = 0.05,pAdjustMethod = "BH",
                               qvalueCutoff = 0.05)




dotplot(result_go_clusterM,showCategory=15,font.size = 10,label_format = 30)
barplot(result_go_clusterM,showCategory=13)
result_go_clusterM<- as.data.frame(result_go_clusterM)
write.table(result_go_clusterM,"result_go_clusterM.xls",quote = F,sep = "\t")

write.table(macrophage_whole.markers,"macrophage_whole.markers.xls",quote = F,sep = "\t")

####KEGG富集分析
#1、KEGG富集
kk <- enrichKEGG(gene = clusterM.markers$ENTREZID,keyType = "kegg",organism= "mouse", qvalueCutoff = 0.05, pvalueCutoff=0.05)
kk@result$Description <- gsub("\\s*\\-\\s*Mus musculus.*$", "", kk@result$Description)

#2、可视化
###柱状图
hh <- as.data.frame(kk)#自己记得保存结果哈！
write.table(hh,"result_enrichKEGG_clusterD.xls",quote = F,sep = "\t")

rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
p= ggplot(hh,aes(y=order,x=Count,fill=p.adjust))+
  geom_bar(stat = "identity",width=0.7)+####柱子宽度
  #coord_flip()+##颠倒横纵轴
  scale_fill_gradient(low = "red",high ="blue" )+#颜色自己可以换
  labs(title = "KEGG Pathways Enrichment",
       x = "Gene numbers", 
       y = "Pathways")+
  theme(axis.title.x = element_text(face = "bold",size = 16),
        axis.title.y = element_text(face = "bold",size = 16),
        legend.title = element_text(face = "bold",size = 16))+
  theme_bw()
ggsave(p,filename = "barplot_enrichKEGG_clusterD.pdf",width =6,height = 9)



###气泡图
hh <- as.data.frame(kk)
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
p=ggplot(hh,aes(y=order,x=Count))+
  geom_point(aes(size=Count,color=-1*p.adjust))+# 修改点的大小
  scale_color_gradient(low="green",high = "red")+
  labs(color=expression(p.adjust,size="Count"), 
       x="Gene Number",y="Pathways",title="KEGG Pathway Enrichment")+
  theme_bw()
ggsave(p,filename = "Dotplot_enrichKEGG_clusterE.pdf",width =6,height = 9)


}






###############STEP 4   拟时序分析##################

{# 安装和加载monocle 2包

library(monocle)

# 读入单细胞测序数据和基因注释文件
countData <- read.table("path/to/countData.txt", header = TRUE, row.names = 1)
geneAnnotation <- read.table("path/to/geneAnnotation.txt", header = TRUE, row.names = 1)

# 创建CellDataSet对象
cds <- newCellDataSet(as.matrix(countData), phenoData = data.frame(condition = colnames(countData)), featureData = geneAnnotation)

# 进行基因筛选和归一化
cds <- preprocessCDS(HSMM, num_dim = 30, genes_use = rownames(geneAnnotation))

cds <- orderCells(cds, reverse = TRUE)
cds <- differentialGeneTest(cds, fullModelFormulaStr = "~condition")

# 可视化结果
plot_cells(cds, color_cells_by = "condition", label_cells_by = "condition")
plot_genes_branched_pseudotime(cds, "gene_of_interest")

}









###既然作者已经对巨噬细胞进行了recluster，那就直接用，然后做拟时序分析
VlnPlot(macrophage_whole, features = c("Top2a", "Gpnmb","Serpinb10"))

{
  #plot3<-DimPlot(macrophage_whole, reduction = "umap",group.by = 'ident',pt.size = 0.1)
##画出来的图奇丑无比，这篇文章也没有用tsne降维。另外把ENPP2划给了间质细胞，在巨噬细胞里根本找不到。


markers_df <- FindMarkers(object = macrophage_whole, ident.1 = "I_Macrophage", min.pct = 0.25)
print(x=head(markers_df))
markers_genes= rownames(head(markers_df,n=5))

DefaultAssay(object = macrophage_whole) <- "RNA"
DotPlot(macrophage_whole, features =markers_genes)
FeaturePlot(macrophage_whole, features =markers_genes)
}



###Monocle分析

{########## Create Monocle Object
  require(monocle) 
  require(openxlsx)
  require(Seurat)
  library(tidyverse)
  
  ## Seurat Data Analysis workflow
  # 检查变量是否已经添加到Seurat对象中
  "seurat_clusters" %in% colnames(macrophage_whole@meta.data)
  
  ###如果返回FALSE，则说明seurat_clusters变量没有被添加到SC1对象中。
  # 获取Seurat对象的聚类标记
  seurat_clusters <- Idents(macrophage_whole)
  table(seurat_clusters)
  
  # 变量添加到Seurat对象中
  macrophage_whole <- AddMetaData(macrophage_whole, metadata = seurat_clusters, col.name = "seurat_clusters")
  "seurat_clusters" %in% colnames(macrophage_whole@meta.data) #TRUE
  
  # 运行subset函数
  SObjectNew <- subset(macrophage_whole,seurat_clusters=="E"|seurat_clusters=="P"|seurat_clusters=="M")
  p1=DimPlot(SObjectNew, reduction = "umap",label=TRUE)
  #ggsave(p1,filename = "DimPlot.pdf",width = 6,height = 10)
  
  ########## Create Monocle Object
  Counts<-GetAssayData(object =SObjectNew, slot = 'counts',assay="RNA")
  gene_annotation <- data.frame(gene_short_name=rownames(Counts),biotype=rep("protein_coding",nrow(Counts)))
  rownames(gene_annotation)<-rownames(Counts)
  sample_sheet<-as.data.frame(SObjectNew@meta.data)
  sample_sheet<-cbind(rownames(sample_sheet),sample_sheet) #cbind函数
  colnames(sample_sheet)[1]<-"cell"
  sample_sheet<-as.data.frame(t(subset(t(sample_sheet),select=colnames(Counts))))##进行一下对齐，不对齐会报错。
  
  pd <- new("AnnotatedDataFrame", data = sample_sheet) ##pData可以查看注释信息，fDataK可以查看特征注释信息。
  fd <- new("AnnotatedDataFrame", data = gene_annotation)
  HSMM <- newCellDataSet(Counts,phenoData = pd, featureData = fd,expressionFamily = negbinomial.size())
  
  ########################
  
  HSMM <- estimateSizeFactors(HSMM)  ##尺度因子的估计，用来下游归一化；
  HSMM <- estimateDispersions(HSMM)  ##离散度的计算，用Monocle 不要装R 4.2以上的版本，有坑！！
  
  
  
  ######### filter Genes
  HSMM <- detectGenes(HSMM, min_expr = 1.0)  ##基因过滤
  expressed_genes <- row.names(subset(fData(HSMM),num_cells_expressed >= 50)) 
  length(expressed_genes) ##查看向量长度
  
  ######### Cluster DEGs
  
  diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],fullModelFormulaStr = "~seurat_clusters")
  
  head(diff_test_res)
  length(diff_test_res)
  ordering_genes <- row.names (subset(diff_test_res, qval < 0.001))
  HSMM <- setOrderingFilter(HSMM, ordering_genes)
  
  ################### assign color 
  library("scales")
  cluster_vec <- as.character(seq(from=0,to=10,by=1))
  hex_codes2 <- hue_pal()(length(cluster_vec))  #ggplot和pheatmap自带了绘图的包。这里是给细胞群赋值颜色，方便后面改颜色
  names(hex_codes2)<-cluster_vec
  hex_codes2<-hex_codes2[1:7]
  plotg<-plot_ordering_genes(HSMM)
  ggsave(filename = "plot_ordering_genes.png",width = 3,height = 4,plot = plotg)  #保存结果
  
  
  ##############################################
  #monocle降维，用的DDRTree，这里用的是Monocle 2版本。
  HSMM_myo <- reduceDimension(HSMM, max_components = 2, method = 'DDRTree', k = 10, n_iters = 100, init_method = 'pca')
  HSMM_myo <- orderCells(HSMM_myo)  #进行cell的指定
  p1<-plot_cell_trajectory(HSMM_myo, color_by = "seurat_clusters")+theme(legend.position="right")+ scale_color_manual(values = hex_codes2, name = "seurat_clusters")
  ggsave(filename = "plot_cell_trajectory.seurat_clusters.png",width = 6,height = 6,plot = p1)
  ##按时间
  p<-plot_cell_trajectory(HSMM_myo, color_by = "Pseudotime")+theme(legend.position="right")
  ggsave(filename = "plot_cell_trajectory.Pseudotime.png",width = 5,height = 4,plot = p)
  ##按state
  p<-plot_cell_trajectory(HSMM_myo, color_by = "State")# +facet_wrap(~State, nrow = 1)  ##facet_wrap 可以单独分面展示
  p<-p+theme(legend.position="right")
  ggsave(filename = "plot_cell_trajectory.State.png",width = 5,height = 4,plot = p)
  ##按分群
  p<-plot_cell_trajectory(HSMM_myo, color_by = "seurat_clusters") +facet_wrap(~seurat_clusters, nrow = 1)
  p<-p+theme(legend.position="right")+ scale_color_manual(values = hex_codes2, name = "seurat_clusters")
  ggsave(filename = "plot_cell_trajectory.seurat_clusters.facet_wrap.png",width = 13,height = 4,plot = p)
  
  ##计算出与假时间相关的基因
  ordering_genes <- row.names (subset(diff_test_res, qval < 0.0001))
  set.seed(123)
  head(pData(HSMM))
  
  marker_genes<-ordering_genes[sample(1:length(ordering_genes),20)]
  diff_test_res <- differentialGeneTest(HSMM_myo[marker_genes,],fullModelFormulaStr = "~sm.ns(Pseudotime)")
  sig_gene_names <- row.names(subset(diff_test_res, qval < 0.01))
  
  png(filename = "plot_pseudotime_heatmap.seurat_clusters.png",width = 1000,height = 1200,res = 300)
  plot_pseudotime_heatmap(HSMM_myo[sig_gene_names,],num_clusters = 3,cores = 1,show_rownames = T)
  dev.off()
  ##基因趋势变化的散点图
  cds_subset <- HSMM_myo[c("MKRN1", "TCEAL9","CD53" ,"CD3D"),]
  p<-plot_genes_in_pseudotime(cds_subset, color_by = "seurat_clusters",ncol=2)+ scale_color_manual(values = hex_codes2, name = "seurat_clusters")
  ggsave(filename = "plot_genes_in_pseudotime.seurat_clusters.png",width = 8,height = 4,plot = p)
  
  ######## BEAM找节点
  
  marker_genes<-ordering_genes[sample(1:length(ordering_genes),100)]  ##筛选出来与节点1相关的基因集，核心数用的1
  BEAM_res <- BEAM(HSMM_myo[marker_genes,], branch_point = 1, cores = 1)
  BEAM_res <- BEAM_res[order(BEAM_res$qval),]
  BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
  ##用Png函数进行保存。可以画热图和散点图。要用dev.off结束保存进程才能看到图
  png(filename = "plot_genes_branched_heatmap.seurat_clusters.png",width = 1000,height = 1200,res = 300)
  plot_genes_branched_heatmap(HSMM_myo[row.names(subset(BEAM_res,qval < 1e-2)),],branch_point = 1,num_clusters = 3,cores = 1,use_gene_short_name = T,show_rownames = T)
  dev.off()
  p<-plot_genes_branched_pseudotime(HSMM_myo[c("CD53", "TGM2","RGCC"),], color_by = "seurat_clusters",ncol=3)+ scale_color_manual(values = hex_codes2, name = "seurat_clusters")
  ggsave(filename = "plot_genes_branched_pseudotime.png",width = 8,height = 4,plot = p)
  
  ####### 
  
  ##个性化绘图，大数据集建议用monocle 3（几万个细胞的时候），或者用Python 的RNA 速率包。
  hex_codes2<-hex_codes2[1:2] ##重新指定想保留的2个cluster
  p<-plot_cell_trajectory(HSMM_myo, color_by = "seurat_clusters")+theme(legend.position="right")+ scale_color_manual(values = hex_codes2, name = "seurat_clusters")
  ggsave(filename = "plot_cell_trajectory1.png",width = 5,height = 4,plot = p)
  
  p<-plot_cell_trajectory(HSMM_myo, color_by = "seurat_clusters")+theme(legend.position="right")+ scale_color_manual(values = hex_codes2, name = "seurat_clusters",na.value = NA)##na.value = NA可以去掉你不想要的那个cluster。
  ggsave(filename = "plot_cell_trajectory2.png",width = 5,height = 4,plot = p)
  
  #########################
  ##画树状图
  pcct<-plot_complex_cell_trajectory(HSMM_myo, color_by = 'seurat_clusters', cell_size = 0.5, cell_link_size = 0.3)  + scale_size(range = c(0.2, 0.2)) +scale_color_manual(values = hex_codes2, name = "seurat_clusters",na.value = NA)+theme(axis.text.x = element_text(angle = 30, hjust = 1))+facet_wrap(~seurat_clusters, nrow = 1) #+ theme (legend.position="none", legend.title=element_blank())
  ggsave(filename = "plot_complex_cell_trajectory2.png",width = 12,height = 4,plot = pcct)
  
  ############################
  ##提取某一个单独的cluster画热图
  HSMM_myoSub<-HSMM_myo[,rownames(pData(HSMM_myo))[pData(HSMM_myo)[["seurat_clusters"]]=="1"]]
  plot_pseudotime_heatmap(HSMM_myoSub[sig_gene_names,],num_clusters = 3,cores = 1,show_rownames = T)
  
  ##展示某一个基因沿着发育轨迹的表达情况
  pData(HSMM_myo)$Genes<-HSMM_myo@assayData$exprs[rownames(HSMM_myo@assayData$exprs)=="CD53",]
  p<-plot_cell_trajectory(HSMM_myo, color_by = "Genes")+theme(legend.position="right")+scale_colour_gradient(low = "#46B8DA99",high = "#D43F3AFF")
  ggsave(filename = "plot_cell_trajectoryGene.png",plot = p,width = 5,height = 4)
  
}






{scRNAsub <- FindVariableFeatures(SC2, selection.method = "vst", nfeatures = 2000)
scale.genes <-  rownames(scRNAsub)

scRNAsub <- ScaleData(scRNAsub, features = scale.genes)

###非线性降维
scRNAsub <- RunPCA(scRNAsub, features = VariableFeatures(scRNAsub))

###制作滚石图，选择合适的PC
ElbowPlot(scRNAsub, ndims=20, reduction="pca")
pc.num=1:10


head(SC2@meta.data)


###标准后续分析流程
scRNAsub <- FindNeighbors(scRNAsub, dims = pc.num)
scRNAsub <- FindClusters(scRNAsub, resolution = 0.5)
scRNAsub = RunUMAP(scRNAsub,reduction="pca", dims = pc.num)
DimPlot(scRNAsub,reduction="umap", pt.size=0.1,label=TRUE,repel=TRUE)



table(scRNAsub@meta.data$seurat_clusters)


metadata <- scRNAsub@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,'cell_cluster.csv',row.names = F)


###非线性降维
scRNAsub = RunTSNE(scRNAsub, dims = pc.num)
embed_tsne <- Embeddings(scRNAsub, 'tsne')
write.csv(embed_tsne,'embed_tsne.csv')
plot1 = DimPlot(scRNAsub, reduction = "tsne")
ggsave("immune_recluster_tSNE.pdf", plot = plot1, width = 8, height = 7)


diff.wilcox = FindAllMarkers(scRNAsub)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>%dplyr:: top_n(n = 10, wt = avg_log2FC)
write.csv(all.markers, "diff_genes_wilcox.csv", row.names = F)
write.csv(top10, "top10_diff_genes_wilcox.csv", row.names = F)

save.image(file = "suRNAsub.RData")}

