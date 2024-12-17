################################################################################
# step1 subset
################################################################################
library(Seurat)
rm(list=ls())
options(stringsAsFactors = F)
load(file = "./rdata/scRNA_all_with_mac.rda")


scRNA_myelo=subset(scRNA,Level0 %in% c('Mac'))

dim(scRNA_myelo)
# [1] 23657   941

table(scRNA_myelo@meta.data$condition)
# D   E   M   P 
# 282 266 172 221 

## 重新聚类分群
set.seed(1314)
scRNA_myelo <- FindVariableFeatures(scRNA_myelo, selection.method = "vst", nfeatures = 2000) 
scRNA_myelo<- ScaleData(scRNA_myelo)
scRNA_myelo<- RunPCA(scRNA_myelo) 
plot1 <- DimPlot(scRNA_myelo, reduction = "pca", group.by="condition") 
plot2 <- ElbowPlot(scRNA_myelo, ndims=30, reduction="pca") 
plotc <- plot1+plot2
plotc


pc.num=1:10
set.seed(1314)
scRNA_myelo <- FindNeighbors(scRNA_myelo, dims = pc.num) 
scRNA_myelo <- FindClusters(scRNA_myelo,resolution = 0.4)
table(scRNA_myelo@meta.data$seurat_clusters)


#UMAP可视化
scRNA_myelo <- RunUMAP(scRNA_myelo, dims = pc.num)
DimPlot(scRNA_myelo, reduction = "umap",label=T,group.by = "condition") + DimPlot(scRNA_myelo, reduction = "umap",label=T)
DimPlot(scRNA_myelo, reduction = "umap",label=T,group.by = "condition")+ DimPlot(scRNA_myelo, reduction = "umap",label=T,group.by = "Level2")

FeaturePlot(scRNA_myelo,features = c("Gpnmb","Cd14","Top2a","S100a4"))
DimPlot(scRNA_myelo)

#11.1 识别每个类群的全部标记物【限速步骤】
cluster_markers <- FindAllMarkers(object = scRNA_myelo, #Seurat对象；
                                  test.use="wilcox",#检验方法
                                  only.pos = TRUE,#仅返回表达倍数大于0的基因（默认为 FALSE）
                                  logfc.threshold = 0.25,
) 
mt_genes <- row.names(cluster_markers)[grepl("^(MT|mt|Mt)-", row.names(cluster_markers))]
rps_genes <- row.names(cluster_markers)[grepl("^RPS", row.names(cluster_markers))]
mrp_genes <- row.names(cluster_markers)[grepl("^MRP", row.names(cluster_markers))]
rpl_genes <- row.names(cluster_markers)[grepl("^RPL", row.names(cluster_markers))]
filter_genes <- c(mt_genes,rps_genes, mrp_genes, rpl_genes)
a<-rownames(cluster_markers)%in%filter_genes
table(a)
cluster_markers<-cluster_markers[!a,]
#11.2 筛选每个cluster中表达前5的基因
# 筛选p_val<0.05的基因
library(tidyverse)
all.markers = cluster_markers %>% 
  dplyr::select(gene,everything()) %>%
  dplyr::filter(p_val<0.05)
write.csv(all.markers,"./files/dif_between_mac_cluster13.csv")



# 将avg_log2FC排名前5的基因筛选出来
top5 = all.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC)
write.table(top5,"./files/sig_dif_top5_mac_cluster13.xls",row.names = T,col.names = T,quote = F,sep = "\t")


FeaturePlot(scRNA_myelo,features = c("SPP1","TREM2"))

#  Need to look at the picture, determine the cell subsets:
celltype=data.frame(ClusterID=0:5,
                    celltype='unkown')
celltype[celltype$ClusterID %in% c(0,2),2]='Mono' 
celltype[celltype$ClusterID %in% c(5),2]='unkown' 
celltype[celltype$ClusterID %in% c(1),2]='Gpnmb+Mac'
celltype[celltype$ClusterID %in% c(3),2]='S100a4+Mac'
celltype[celltype$ClusterID %in% c(4),2]='Top2a+Mac'

table(celltype$celltype)
new.cluster.ids <- celltype$celltype
names(new.cluster.ids) <- levels(scRNA_myelo)
scRNA_myelo <- RenameIdents(scRNA_myelo, new.cluster.ids)
scRNA_myelo@meta.data$Myeloid_subtype <- Idents(scRNA_myelo)
table(scRNA_myelo@meta.data$Myeloid_subtype)


DimPlot(scRNA_myelo,label = T)

scRNA_myelo=subset(scRNA_myelo, Myeloid_subtype != c("unkown")) 
table(scRNA_myelo@meta.data$Myeloid_subtype)

DimPlot(scRNA_myelo, label = T,
        reduction = 'umap') 


scRNA_myelo@meta.data$Myeloid_subtype = droplevels(scRNA_myelo@meta.data$Myeloid_subtype, exclude = setdiff(levels(scRNA_myelo@meta.data$Myeloid_subtype),unique(scRNA_myelo@meta.data$Myeloid_subtype)))
table(scRNA_myelo@meta.data$Myeloid_subtype)


save(scRNA_myelo,file="./rdata/ov_scRNA_myelo.rda")

################################################################################
# step2 hvg and heatmap
################################################################################
rm(list = ls())
load("./rdata/ov_scRNA_myelo.rda")
scRNA_myelo@meta.data$condition = factor(scRNA_myelo@meta.data$condition,levels = c("P","E","M","D"))
# FeaturePlot(scRNA_myelo,features = c("VCAN","SPP1","TREM2","CXCL10","CXCL9","FOLR2","LYVE1","SELENOP","MKI67","TCOF1"))
# 
# FeaturePlot(scRNA_myelo,features = c("SPP1","IGF1","CTSD","CTSB")) 
# DotPlot(scRNA_myelo,features = c("SPP1","IGF1","CTSD","CTSB","CDKN2A"))
# 
# 
# 
# scRNA_myelo=subset(scRNA_myelo, Myeloid_subtype != c("DC")) 
# scRNA_myelo=subset(scRNA_myelo, Myeloid_subtype != c("Monocyte")) 
# scRNA_myelo@meta.data$Myeloid_subtype = droplevels(scRNA_myelo@meta.data$Myeloid_subtype, exclude = setdiff(levels(scRNA_myelo@meta.data$Myeloid_subtype),unique(scRNA_myelo@meta.data$Myeloid_subtype)))
# table(scRNA_myelo@meta.data$Myeloid_subtype)

DimPlot(scRNA_myelo, label = T,
        reduction = 'umap') 
source("utools.R")

p = DimPlot(scRNA_myelo, group.by="Myeloid_subtype", label=F, label.size=4.5,pt.size = 1.5, reduction='umap',cols = ov_palette[3:10] )
ggsave(p,filename=paste0("Mac",'/umap_anno.pdf'),width = 7,height = 6)
p
p = CellRatioPlot(object = scRNA_myelo,
                  sample.name = "condition",
                  celltype.name = "Myeloid_subtype",
                  fill.col = ov_palette[3:10] )
ggsave(p,filename=paste0("Mac",'/ratioplot.pdf'),width = 5,height = 6)
p
cluster_markers <- FindAllMarkers(object = scRNA_myelo, #Seurat对象；
                                  test.use="wilcox",#检验方法
                                  only.pos = TRUE,#仅返回表达倍数大于0的基因（默认为 FALSE）
                                  logfc.threshold = 0.25,
) 
mt_genes <- row.names(cluster_markers)[grepl("^(MT|mt|Mt)-", row.names(cluster_markers))]
rps_genes <- row.names(cluster_markers)[grepl("^RPS", row.names(cluster_markers))]
mrp_genes <- row.names(cluster_markers)[grepl("^MRP", row.names(cluster_markers))]
rpl_genes <- row.names(cluster_markers)[grepl("^RPL", row.names(cluster_markers))]
filter_genes <- c(mt_genes,rps_genes, mrp_genes, rpl_genes)
a<-rownames(cluster_markers)%in%filter_genes
table(a)
cluster_markers<-cluster_markers[!a,]
library(tidyverse)
all.markers = cluster_markers %>% 
  dplyr::select(gene,everything()) %>%
  dplyr::filter(p_val<0.05)
write.csv(all.markers,"./files/dif_between_mac.csv")



# 将avg_log2FC排名前5的基因筛选出来
top5 = all.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC)
write.table(top5,"./files/sig_dif_top5_mac.xls",row.names = T,col.names = T,quote = F,sep = "\t")

# 将avg_log2FC排名前5的基因筛选出来
top30 = all.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 30, wt = avg_log2FC)

source("utools.R")
genes = top30$gene
# genes[30] = "SPP1"
{
  ##美化热图
  library(RColorBrewer)
  mycolor = ov_palette[3:10]
  # mypalette <- brewer.pal(n = 9, name = "YlOrRd")
  p <- DoHeatmap(scRNA_myelo,
                 features = as.character(unique(genes)),
                 group.by = "Myeloid_subtype",
                 group.colors =mycolor)+
    scale_fill_gradientn(colors = c("white","grey90","firebrick3"))
  
  ggsave("./Mac/Top30DEG_seurat_cluster_heatmap.pdf", p, width = 8, height = 6)
}


############################
## PEMD
table(scRNA_myelo@active.ident)
Idents(scRNA_myelo) = "condition"

cluster_markers <- FindAllMarkers(object = scRNA_myelo, #Seurat对象；
                                  test.use="wilcox",#检验方法
                                  only.pos = TRUE,#仅返回表达倍数大于0的基因（默认为 FALSE）
                                  logfc.threshold = 0.25,
) 
mt_genes <- row.names(cluster_markers)[grepl("^(MT|mt|Mt)-", row.names(cluster_markers))]
rps_genes <- row.names(cluster_markers)[grepl("^Rps", row.names(cluster_markers))]
mrp_genes <- row.names(cluster_markers)[grepl("^Mrp", row.names(cluster_markers))]
rpl_genes <- row.names(cluster_markers)[grepl("^Rpl", row.names(cluster_markers))]
filter_genes <- c(mt_genes,rps_genes, mrp_genes, rpl_genes)
a<-rownames(cluster_markers)%in%filter_genes
table(a)
cluster_markers<-cluster_markers[!a,]
library(tidyverse)
all.markers = cluster_markers %>% 
  dplyr::select(gene,everything()) %>%
  dplyr::filter(p_val<0.05)
write.csv(all.markers,"./files/dif_between_mac_pemd.csv")



# 将avg_log2FC排名前5的基因筛选出来
top5 = all.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC)
write.table(top5,"./files/sig_dif_top5_mac_pemd.xls",row.names = T,col.names = T,quote = F,sep = "\t")



source("utools.R")
genes = top5$gene
# genes[30] = "SPP1"
{
  ##美化热图
  library(RColorBrewer)
  mycolor = ov_palette[3:10]
  # mypalette <- brewer.pal(n = 9, name = "YlOrRd")
  p <- DoHeatmap(scRNA_myelo,
                 features = as.character(unique(genes)),
                 group.by = "condition",
                 group.colors =mycolor)+
    scale_fill_gradientn(colors = c("white","grey90","firebrick3"))
  
  ggsave("./Mac/Top30DEG_seurat_cluster_heatmap_pemd.pdf", p, width = 8, height = 6)
}


################################################################################
# step3 scoring
################################################################################
scRNA = scRNA_myelo
rm(scRNA_myelo)
table(scRNA@meta.data$Myeloid_subtype)







################################################################################
# step3 scMetabolism
################################################################################
{
  library(Seurat)
  library(scMetabolism)
  library(ggplot2)
  library(rsvd)
}


DimPlot(scRNA, group.by = "Myeloid_subtype", 
        split.by = "condition", ncol = 4,pt.size = 0.2)


# devtools::install_github("YosefLab/VISION@v2.1.0")

countexp.Seurat <-scRNA
source("utools.R")
mouse_results <- Mouse.sc.metabolism(countexp.Seurat, metabolism.type = 'KEGG')
# countexp.Seurat@meta.data$cell_type <- countexp.Seurat@active.ident

# countexp.Seurat<-sc.metabolism.Seurat(obj = countexp.Seurat, method = "AUCell",rna_slot = "RNA",
#                                       imputation = F, ncores = 2, metabolism.type = "KEGG")
save(mouse_results,file="rdata/ov_mac_scmetabolism.rda")
countexp.Seurat = mouse_results
metabolism.matrix <- countexp.Seurat@assays$METABOLISM$score
metabolism.matrix



library(ggplot2)
p.dimplot<-DimPlot.metabolism(obj = countexp.Seurat, 
                              pathway = "Arachidonic acid metabolism", 
                              dimention.reduction.type = "umap", 
                              dimention.reduction.run = F, size = 1)
p.dimplot

input.pathway <- rownames(countexp.Seurat@assays[["METABOLISM"]][["score"]])[1:85] # 通路序号更改
p_dotplot <- DotPlot.metabolism(obj = countexp.Seurat, 
                                pathway = input.pathway, 
                                phenotype = "Myeloid_subtype", 
                                norm = "y")
p_dotplot
# ggsave(p_dotplot,filename=paste0("scMetabolism",'/dotplot_celltype_metabolism.png'),width = 12,height = 9)


p_dotplot <- DotPlot.metabolism(obj = countexp.Seurat, 
                                pathway = input.pathway, 
                                phenotype = "condition", 
                                norm = "y")
p_dotplot
# ggsave(p_dotplot,filename=paste0("scMetabolism",'/dotplot_disease_metabolism.png'),width = 12,height = 9)


sce_Metal_exp = countexp.Seurat
mscore_data = data.frame(t(sce_Metal_exp@assays[["METABOLISM"]][["score"]]),sce_Metal_exp$Myeloid_subtype)
avg_sM=aggregate(mscore_data[,1:ncol(mscore_data)-1],list(mscore_data$sce_Metal_exp.Myeloid_subtype),mean)
rownames(avg_sM) = avg_sM$Group.1
avg_sM=data.frame(t(avg_sM[,-1]))
avg_sM$KEGG = rownames(sce_Metal_exp@assays[["METABOLISM"]][["score"]])
rownames(avg_sM)=avg_sM$KEGG

c_k_l = c()
for(c in c(1:ncol(avg_sM))){
  c_k=avg_sM[order(avg_sM[,c]),]$KEGG[1:5]
  c_k_l=c(c_k_l,c_k)
}
c_k_l= unique(c_k_l)
c_k_d = avg_sM[avg_sM$KEGG %in%c_k_l,]

rownames(c_k_d) = c_k_d$KEGG
p_heatmap<-pheatmap::pheatmap(c_k_d[,-ncol(c_k_d)],show_colnames = T,scale='row')
p_heatmap
# ggsave(p_heatmap,filename=paste0("scMetabolism",'/heatmap_mac_metabolism.png'),width = 12,height = 9)




################################################################################
# step4 ## 拟时序分析 slingshot and monocle2
################################################################################
## 拟时序分析
library(Seurat)
#没有monocle要先安装 BiocManager::install
library(monocle)
library(tidyverse)
library(patchwork)
library(Seurat)
rm(list = ls())
gc()
load("./rdata/ov_scRNA_myelo.rda")
##开始monocle
library(monocle)
# setwd('./myelo_pseudotime/')
scRNAsub <- scRNA_myelo

#scRNAsub <- readRDS("subcluster/scRNAsub.rds")  #scRNAsub是上一节保存的Mastcell细胞子集seurat对象
data <- as(as.matrix(scRNAsub@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNAsub@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
## 以下代码一律不得修改
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())


mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)


##使用monocle选择的高变基因，不修改
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)
p3 <- plot_ordering_genes(mycds)
p3

#降维 我们将减少空间到一个二维，我们将能够轻松地可视化和解释
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
#排序
# trace("project2MST", edit = T, where = asNamespace("monocle"))
mycds <- orderCells(mycds)
mycds <- orderCells(mycds,root_state = 3)


## 查看关键差异基因在拟时序中所处的时间
load("./rdata/monocle2_mac.rda")
# library(readr)
# scenscence_marker <- read_csv("files/scenscence_marker.csv")
dev.off()
# my_pseudotime_cluster <- plot_pseudotime_heatmap(mycds[c(scenscence_marker$...1),],
#                                                  # num_clusters = 2, # add_annotation_col = ac,
#                                                  show_rownames = TRUE,
#                                                  return_heatmap = TRUE)
# 
# my_pseudotime_cluster 





source("utools.R")
#State轨迹分布图
# plot_cell_trajectory(monocle_data, color_by = "angio_score1") +
#   scale_color_gradient2 (low="#fcfbfd", mid="#9e9ac8", high="#3f007d")
# plot_cell_trajectory(
#   monocle_data, color_by = "angio_score1") + scale_color_manual(breaks = waiver(),values=coul)
plot1 <- plot_cell_trajectory(mycds, color_by = "State",)
plot1
# ggsave("./State.pdf", plot = plot1, width = 6, height = 5)
# ggsave("./State.png", plot = plot1, width = 6, height = 5)


##Cluster轨迹分布图
mycol =c(
  "Mono"=     '#E069A6',
  "Gpnmb+Mac"=    "#f1707d",
  "S100a4+Mac"=   "#AFC2D9",
  "Top2a+Mac"=  "#6894B9"
)
plot2 <- plot_cell_trajectory(mycds, color_by = "Myeloid_subtype")+scale_color_manual(values = mycol)
plot2
# ggsave("./Cluster.pdf", plot = plot2, width = 10, height = 5)
# ggsave("./Cluster.png", plot = plot2, width = 10, height = 5)


##Pseudotime轨迹图
library(viridis)
plot3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")+scale_colour_gradientn(colours = magma(10)[3:10])
plot3
# ggsave("./Pseudotime.pdf", plot = plot3, width = 6, height = 5)
# ggsave("./Pseudotime.png", plot = plot3, width = 6, height = 5)

plot4 <-  plot_cell_trajectory(mycds, color_by = "condition")+scale_color_manual(values = ov_palette)
plot4
# ggsave("./group.pdf", plot = plot4, width = 6, height = 5)
# ggsave("./group.png", plot = plot4, width = 6, height = 5)


source("utools.R")
mat<- as.matrix(scRNA_myelo@assays$RNA@counts)
results <- CytoTRACE(mat = mat)
phe <- scRNA_myelo$Myeloid_subtype
phe = as.character(phe)
names(phe) <- rownames(scRNA_myelo@meta.data)
# plotCytoGenes(results, numOfGenes = 10)
# plotCytoTRACE(results, phenotype = phe,gene = "CCR7",outputDir = "DC_trace")#单个基因的表达量展示
scRNA_myelo@meta.data$cytotrace=results[["CytoTRACE"]]
FeaturePlot(scRNA_myelo,features = c("cytotrace"))

mycds$cytotrace=results[["CytoTRACE"]]#添加cytotrace结果到cds对象中
cols <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
# pdf("cytotrace.pdf",6,6)
plot5 = plot_cell_trajectory(mycds,color_by="cytotrace", size=1,show_backbone=TRUE)+
  scale_colour_gradientn(colours = cols)
plot5
# dev.off()
save(mycds,results,file = "rdata/monocle2_mac.rda")
load("./rdata/monocle2_mac.rda")
##合并出图
#fig3.BC
plotc <- (plot3|plot1)/(plot4|plot2)
plotc

ggsave("./Mac/monocle2.pdf", plot = plotc, width =8, height = 8)

# write.csv(pData(mycds), "./lesson/result/pseudotime_mac.csv")


# heatmap_with_mac_genes

genes_table = data.table::fread("./files/dif_between_mac.csv")
diff.genes <- genes_table
sig_diff.genes <- subset(diff.genes,p_val_adj<0.0001&abs(avg_log2FC)>1)$gene
sig_diff.genes <- unique(as.character(sig_diff.genes))
diff_test <- differentialGeneTest(mycds[sig_diff.genes,], cores = 1, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test, qval < 0.01))

p1 = plot_pseudotime_heatmap(mycds[sig_gene_names,], num_clusters=6,
                             show_rownames=T, return_heatmap=T)

# dev.off()
ggsave("./Mac/monocle2_heatmap.pdf", plot = p1, width = 5, height = 7)


#  two fate
pseudotime_de <- diff_test[order(diff_test$qval), ]
states_de <- differentialGeneTest(mycds[sig_diff.genes,],
                                  fullModelFormulaStr = "~State")
states_de <- states_de[order(states_de$qval), ]

source("./BEAM.R")
BEAM_res=BEAM(mycds,branch_point = 2,cores = 1)
#会返回每个基因的显著性，显著的基因就是那些随不同branch变化的基因
#这一步很慢
BEAM_res=BEAM_res[,c("gene_short_name","pval","qval")]
saveRDS(BEAM_res, file = "BEAM_res.rds")
画图

tmp1=plot_genes_branched_heatmap(test[row.names(subset(BEAM_res,qval<1e-4)),],
                                 branch_point = 1,
                                 num_clusters = 4, #这些基因被分成几个group
                                 cores = 1,
                                 branch_labels = c("Cell fate 1", "Cell fate 2"),
                                 #hmcols = NULL, #默认值
                                 hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
                                 branch_colors = c("#979797", "#F05662", "#7990C8"), #pre-branch, Cell fate 1, Cell fate 2分别用什么颜色
                                 use_gene_short_name = T,
                                 show_rownames = F,
                                 return_heatmap = T #是否返回一些重要信息
)

pdf("branched_heatmap.pdf",width = 5,height = 6)
tmp1$ph_res
dev.off()


# 获取state
state=pData(mycds)
identical(rownames(scRNAsub@meta.data),rownames(state))
# 把state加到metadata里面以供差异分析
scRNAsub@meta.data$monocle_state=state$State
Idents(scRNAsub)=Idents(scRNAsub)=scRNAsub@meta.data$monocle_state



library(Seurat)
library(slingshot)
rm(list = ls())
gc()
load("./rdata/ov_scRNA_myelo.rda")
# 我们把 counts 矩阵提取，用于构建 sce 的对象。
# 我们提取高边基因
# counts <- scRNA_myelo@assays$RNA@counts
counts <- counts[rownames(scRNA_myelo@assays$RNA@scale.data),]

# 将表达矩阵转换为SingleCellExperiment对象
# 输入需要counts矩阵，否则影响下游分析
sim <- SingleCellExperiment(assays = List(counts = counts)) 
# 2. 提取其他有用的信息
# 为了与前面画图一致，我们把降维坐标和细胞标签也提取并加入 sce 对象。

# umap reduction
umap = scRNA_myelo@reductions$umap@cell.embeddings
colnames(umap) = c('UMAP-1', 'UMAP-2')
# 将降维的结果添加到SingleCellExperiment对象中
reducedDims(sim) = SimpleList(UMAP = umap)

# metadata
meta = scRNA_myelo@meta.data
# colData(sim)相当于meta.data，但他不是data.frame格式
# 所以需要分步赋予
colData(sim)$sampleId = meta$orig.ident
colData(sim)$celltype = meta$Myeloid_subtype
colData(sim)$seurat_clusters = meta$seurat_clusters
# 我们可以直接用 plot() 函数画一下看看：

rd = umap
plot(rd, col = rgb(0,0,0,.5), pch=16, asp = 1)
DimPlot(scRNA_myelo,group.by = "seurat_clusters")

#   3. 细胞亚群划分

sim <- slingshot(sim, 
                 clusterLabels = 'celltype',  # 选择colData中细胞注释的列名
                 reducedDim = 'UMAP',  
                 start.clus= c("Mono"),  # 选择起点
                 end.clus = NULL     # 这里我们不指定终点
)     
colnames(colData(sim))

# 首先我们看看我们各个 lineage 是什么样的，我们使用官方提供的代码
library(RColorBrewer)
plot(reducedDims(sim)$UMAP, pch=16, asp = 1)
lines(SlingshotDataSet(sim), lwd=2, col=brewer.pal(9,"Set1"))
legend("right",
       legend = paste0("lineage",1:3),
       col = unique(brewer.pal(6,"Set1")),
       inset=0.8,
       pch = 16)


summary(sim$slingPseudotime_2)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  0.000   2.207   4.050   6.339  10.041  18.065    9446 

# 我们做一个绚丽的渐变色彩虹色
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100) # 我们把这些颜色变成100个梯度，模拟渐变色
plotcol <- colors[cut(sim$slingPseudotime_2, breaks=100)] # 这里我们用cut函数把 lineage2 分割成100个区间，同一区间的细胞给与同一个颜色
plotcol[is.na(plotcol)] <- "lightgrey" # 不属于 lineage3 的点为NA，我们把他们变成灰色
plotcol

# pdf(file = "figures/04_lineage2.pdf",width = 7,height = 5)
plot(reducedDims(sim)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sim), lwd=2, col=brewer.pal(9,"Set1"))
legend("right",
       legend = paste0("lineage",1:6),
       col = unique(brewer.pal(6,"Set1")),
       inset=0.8,
       pch = 16)
# dev.off()

coldata <- colData(sim)
coldata <- data.frame(celltype = coldata@listData$celltype, 
                      sampleId = coldata@listData$sampleId,
                      plotcol = plotcol)
rownames(coldata) = sim@colData@rownames

# 我们把 lineage3 信息筛选出来
# 细胞 barcodes
filter_cell <- dplyr::filter(coldata, plotcol != "lightgrey")
filter_cell <- rownames(filter_cell)
head(filter_cell)
# 矩阵 counts
counts <- sim@assays@data@listData$counts
filter_counts <- counts[,filter_cell]
filter_counts[1:5,1:5]

set.seed(111)
scell <- sample(colnames(filter_counts), size = 500)

# 最后的到抽样后的矩阵
filter_counts = filter_counts[, scell]
dim(filter_counts)

# 重新将表达矩阵转换为SingleCellExperiment对象
filter_sim <- SingleCellExperiment(assays = List(counts = filter_counts))

# colData
filter_coldata = colData(sim)[scell, 1:3]
filter_sim@colData = filter_coldata

# reduction
rd = reducedDim(sim)
filter_rd <- rd[scell,]
reducedDims(filter_sim) <- SimpleList(UMAP = filter_rd)

# K-Means
set.seed(111)
cl <- kmeans(filter_rd, centers = 4)$cluster # 可以指定分群数量
head(cl)
colData(filter_sim)$kmeans <- cl

## 可视化聚类分群的结果
library(RColorBrewer)
mycolors = brewer.pal(4,"Set1") # colors
library(ggplot2)
plt_k = data.frame(filter_rd,
                   kmeans_clusters = factor(cl, levels = sort(unique(cl))))
ggplot(plt_k, aes(UMAP.1, UMAP.2))+
  geom_point(aes(color = kmeans_clusters), size = .5)+
  scale_color_manual(values = mycolors)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  xlab('UMAP_1')+
  ylab('UMAP_2')+
  guides(color = guide_legend(override.aes = list(size = 5)) # 调整图例展示的点大小
  )  


# lineage3_kmeans
# 使用 K-Means 划分的亚群进行轨迹推断
filter_sim <- slingshot(filter_sim,
                        clusterLabels = 'kmeans',
                        reducedDim = 'UMAP',
                        start.clus= "3", # 指定起点
                        end.clus = '1' # 指定终点
)
head(colnames(colData(filter_sim)))

plot(reducedDims(filter_sim)$UMAP, pch=16, asp = 1)
lines(SlingshotDataSet(filter_sim), lwd=2, col=mycolors)


devtools::install_github("statOmics/tradeSeq")
# BiocManager::install("tradeSeq")

library(tradeSeq)

# Fit negative binomial model
counts <- filter_sim@assays@data$counts
crv <- SlingshotDataSet(filter_sim)



set.seed(111)
icMat <- evaluateK(counts = counts, 
                   sds = crv, 
                   k = 3:10,    # no more than 12
                   nGenes = 500, # 每个细胞纳入分析的基因数量，默认是500，这里为了节省示例计算时间
                   verbose = T)



# evaluateK
# we pick nknots = 6.
set.seed(111)
pseudotime <- slingPseudotime(crv, na = FALSE)
cellWeights <- slingCurveWeights(crv)
# fit negative binomial GAM
# 2k cells ~13 min
# system.time()这个函数可以计算运行时间
system.time({ 
  sce <- fitGAM(counts = counts, 
                pseudotime = pseudotime, 
                cellWeights = cellWeights,
                nknots = 6, 
                verbose = FALSE)
})
table(rowData(sce)$tradeSeq$converged)

assoRes <- associationTest(sce)
head(assoRes)
# waldStat df pvalue  meanLogFC
# X1110008L16Rik       NA NA     NA 19.0735745
# X1110008P14Rik       NA NA     NA  1.4417585
# X1110019D14Rik       NA NA     NA  1.1954504
# X1110051M20Rik       NA NA     NA  0.9184634
# X1190005I06Rik       NA NA     NA 41.0395194
# X1500009L16Rik       NA NA     NA 36.0545822

startRes <- startVsEndTest(sce)
head(startRes)
#             waldStat df       pvalue logFClineage1
# HES4    1.011937e+00  1 3.144393e-01   -10.1403413
# ISG15   2.027631e+01  1 6.702522e-06     1.6526103
# MXRA8   7.526157e-15  1 9.999999e-01    -0.3659606
# SMIM1   1.267845e+02  1 0.000000e+00     3.3647451
# RBP7    3.080566e-01  1 5.788752e-01   332.3079668
# PLA2G2A 7.526157e-15  1 9.999999e-01    -0.3659606
# startVsEndTest 使用Wald检验来评估零假设，即平滑器起始点（祖细胞群）的平均表达是否等于平滑器终点（分化细胞群）的平均表达。该检验基本上涉及每个亚群的两个平滑器系数的比较。默认情况下，startVsEndTest 函数在所有亚群中执行全局检验（即同时比较起始点和终点），但你也可以通过设置lineages=TRUE来单独评估每个 lineage。

# 在这个例子中，我们已经把 lineage3 单独挑选出来运算，并根据 K-Means 复现了轨迹，因此实际上我们只有一个 lineage 。所以（5）和（6）的运行结果是相似的。如果你有多个 lineage ，应该在这一步加上 lineages=TRUE。但不建议同时有多个轨迹。

# 按相关性进行排序
oStart <- order(startRes$waldStat, decreasing = TRUE)

# 挑选相关性最强的基因，并可视化
sigGeneStart <- names(sce)[oStart[1]]
plotSmoothers(sce, counts, gene = sigGeneStart)

# top gene
# 我们也可以用UMAP图展示

plotGeneCount(crv, counts, gene = sigGeneStart)

# umap top
# （7）美化基因表达量变化图
# 上述两个基因表达量变化图是通过 ggplot2 系统绘制的，我们可以通过添加一些自定义参数来美化他。
# 
# 我们提取需要的数据，并构建成清洁数据框。

# 取celltpye和配色信息
coldata <- data.frame(celltype = sim@colData$celltype,
                      plotcol = plotcol)
rownames(coldata) = colnames(sim)

# 把sce中的3000个细胞对应信息取出
filter_coldata <- coldata[colnames(sce),]

# 添加拟时序信息
filter_coldata$Pseudotime = sce$crv$pseudotime.Lineage1

# top20 genes
top6 <- names(sce)[oStart[1:30]]
top6_exp = sce@assays@data$counts[top6,] 
library(dplyr)
top6_exp = log2(top6_exp + 1) %>% t()

# 获得最终的清洁数据
plt_data = cbind(filter_coldata, top6_exp)
colnames(plt_data)
# [1] "celltype"   "plotcol"    "Pseudotime" "CENPF"      "HBB"       
# [6] "HMGB2"      "CD74"       "HBA1"       "MKI67"  
# 画图

library(ggplot2)
library(ggsci)
library(ggpubr)
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
mycolors = getPalette(length(unique(plt_data$celltype)))

# 为了拼图美观，我们把legend隐藏掉
plt_list = list()
for (gene in top6) {
  print(gene)
  p = ggscatter(data = plt_data,
                x = 'Pseudotime',
                y = gene,
                color = 'celltype',
                size = 0.6)+
    geom_smooth(se = F, color = 'orange')+
    theme_bw()+
    scale_color_manual(values = mycolors)+
    theme(legend.position = 'none')
  plt_list[[gene]] = p
}

library(patchwork)
wrap_plots(plt_list)
ggsave(filename = 'Mac/05_top6_genes.pdf', width = 24, height = 15)

# top genes patchworks
# 这里单独save一张有legend的图
gene = 'Spp1'
p_test = ggscatter(data = plt_data,
                   x = 'Pseudotime',
                   y = gene,
                   color = 'celltype',
                   size = 0.6)+
  geom_smooth(se = F, color = 'orange')+
  theme_bw()+
  scale_color_manual(values = mycolors)+
  theme(legend.position = 'right')+
  guides(color = guide_legend(ncol = 1,  # 将图例分为两列显示
                              override.aes = list(size = 3)) )  # 调整图例展示的点大小
p_test
ggsave(plot = p_test, filename = 'Mac/05_Spp1_legend.pdf', width = 5, height = 5)
# （8）美化拼图









################################################################################
## cellchat



library(Seurat)
rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
load("./rdata/scRNA_estrous.rda")
scRNA@meta.data$condition = factor(scRNA@meta.data$condition,levels = c("P","E","M","D"))
load("./rdata/scRNA_all_with_mac.rda")
table(scRNA@meta.data$Level0)
# Endothelium   Epithelium    Granulosa          Mac   Mesenchyme       Oocyte Other_Immune 
# 2668          845        14231          941         8934           13          382 
scRNA_other = subset(scRNA,Level0 != c("Mac"))
table(scRNA_other@meta.data$Level0)
# Endothelium   Epithelium    Granulosa   Mesenchyme       Oocyte Other_Immune 
# 2668          845        14231         8934           13          382
scRNA_myelo@meta.data$Level0 = scRNA_myelo@meta.data$Myeloid_subtype
scRNA_chat = merge(scRNA_other,scRNA_myelo)
table(scRNA_chat@meta.data$celltype)

# scRNA_chat = subset(scRNA_chat,Type %in% "Cancer")

library(CellChat)
meta =scRNA_chat@meta.data # a dataframe with rownames containing cell mata data
data_input <- as.matrix(scRNA_chat@assays$RNA@data)
#data_input=data_input[,rownames(meta)]
identical(colnames(data_input),rownames(meta))

library(CellChat)
cellchat <- createCellChat(object = data_input, meta = meta, group.by = "Level0")

CellChatDB <- CellChatDB.mouse
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use 

dplyr::glimpse(CellChatDB$interaction)##配体-受体分析
# 提取数据库支持的数据子集
cellchat <- subsetData(cellchat)
# 识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
# 识别配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
# 将配体、受体投射到PPI网络
cellchat <- projectData(cellchat, PPI.human)
unique(cellchat@idents)
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
df.net<- subsetCommunication(cellchat)

write.csv(df.net,file ='./files/cellchat_mac_all.csv',quote=F)

cellchat <- aggregateNet(cellchat)
save(cellchat,file = "./rdata/cellchat_mac.rda")
load("./rdata/cellchat_mac.rda")
#cellchat_OV <-cellchat
groupSize <- as.numeric(table(cellchat@idents))
source("utools.R")
par(mfrow = c(1,2), xpd=TRUE)
#fig3.DE
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, color.use = ov_palette,
                 label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,color.use = ov_palette, label.edge= F, title.name = "Interaction weights/strength")

dev.off()

table(cellchat@idents)
#指定受体-配体细胞类型 "CX3CR1+Mac"  , "CXCL10+Mac", "MK167+Mac", "SELENOP+Mac" ,"SPP1+Mac" ,"TCOF1+Mac" ,   "TIMP1+Mac"
netVisual_bubble(cellchat, sources.use = c( "Gpnmb+Mac","Top2a+Mac","S100a4+Mac"), 
                 targets.use = c("Granulosa","Mesenchyme"),
                 remove.isolate = TRUE)

netVisual_bubble(cellchat, 
                 sources.use = c("Granulosa","Mesenchyme"), 
                 targets.use = c( "Gpnmb+Mac","Top2a+Mac","S100a4+Mac","Mono"),
                 remove.isolate = T)

netVisual_bubble(cellchat, 
                 # sources.use = c( "Epi"), 
                 targets.use = c("Epi"),
                 remove.isolate = FALSE)
