##加载会用到的包
{
  library(formatR)
  library(stringr)
  library(Seurat)
  library(dplyr)
  library(patchwork)
  library(ggplot2)
}


library(future)
plan("multisession", workers = 8)

# 安装Ringo 版本
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("kableExtra", "webshot"))


dir()

markers <- read.csv("Markers.csv",header = TRUE,stringsAsFactors = FALSE)

##====================================================================================##
##===========================Step01：数据预处理=========================================##
##====================================================================================##
###加载包
{
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(ggplot2)
  # install.packages('tidydr')
  library(tidydr)
}

#rm(list = ls())

options(future.globals.maxSize = 20 * 1024^3)#设置了全局变量的最大大小限制为20GB
###### 创建Seurat空转对象
dir()
SectionID = c('CTR','LIP')
# DataDir = c('D:/R Notebook/10X-mouse ovary-my data ananlysis/project/CTR/', 'D:/R Notebook/10X-mouse ovary-my data ananlysis/project/LIP/')
DataDir = c('~/spatial/raw_data/project/CTR/', '~/spatial/raw_data/project/LIP/')

stolist = list()
for(i in seq_along(SectionID)){
  # 读取数据
  sto <- Load10X_Spatial(data.dir = DataDir[i], slice = SectionID[i])
  # spot重命名
  sto <- RenameCells(sto, add.cell.id = SectionID[i])
  # 添加切片ID
  sto$SectionID <- SectionID[i]
  # 计算细胞中线粒体基因比例
  sto[["percent.mt"]] <- PercentageFeatureSet(sto, pattern = "^(MT|mt|Mt)-") #该格式为人鼠通用
  # 计算细胞中核糖体基因比例
  sto[["percent.rb"]] <- PercentageFeatureSet(sto, pattern = "^Rp[sl]")
  stolist[[SectionID[i]]] <- sto
}

### 给列表命名并保存数据
dir.create("./Integrate")
# setwd("../Integrate")


##保存读取的数据
# saveRDS(stolist, "D:/R Notebook/10X-mouse ovary-my data ananlysis/QC/2022-10-18-New-QC/stolist18.qc.rds")  
save(stolist,file = "./Integrate/stolist_qc.rda")
rm(list = ls())
gc()


###### 数据质控
load("Integrate/stolist_qc.rda")
table(stolist$CTR$orig.ident)
# SeuratProject 
# 388 
table(stolist$LIP$orig.ident)
# SeuratProject 
# 403 
### 绘制质控小提琴图
lapply(c("QC","Data","Cluster_Bayes","Cluster_Louvain","Deconvolution","Function","stLearn","SVGs"), 
       function(x){if(!file.exists(x)) dir.create(file.path( "Integrate/",x ))})


setwd("./Integrate/QC")


if(T){
  ### 质控前指标
  # 分子数
  p.list <- lapply(stolist, function(x){
    p1 <- VlnPlot(x, features = "nCount_Spatial") + ggtitle(unique(x$SectionID)) +
      theme(legend.position = "none", axis.text.x = element_blank())
    p2 <- SpatialFeaturePlot(x, features = "nCount_Spatial") + theme(legend.position = "right")
    p <- p1|p2
    p
  })
  p <- wrap_plots(p.list, ncol = 1)
  ggsave("nCount_Spatial_before.pdf", p, width = 10, height = 9)
  # 基因数
  p.list <- lapply(stolist, function(x){
    p1 <- VlnPlot(x, features = "nFeature_Spatial") + ggtitle(unique(x$SectionID)) +
      theme(legend.position = "none", axis.text.x = element_blank())
    p2 <- SpatialFeaturePlot(x, features = "nFeature_Spatial") + theme(legend.position = "right")
    p <- p1|p2
    p
  })
  p <- wrap_plots(p.list, ncol = 1)
  ggsave("nFeature_Spatial_before.pdf", p, width = 10, height = 9)
  # 线粒体
  p.list <- lapply(stolist, function(x){
    p1 <- VlnPlot(x, features = "percent.mt") + ggtitle(unique(x$SectionID)) +
      theme(legend.position = "none", axis.text.x = element_blank())
    p2 <- SpatialFeaturePlot(x, features = "percent.mt") + theme(legend.position = "right")
    p <- p1|p2
    p
  })
  p <- wrap_plots(p.list, ncol = 1)
  ggsave("percent_mt_before.pdf", p, width = 10, height = 9)
  # 核糖体
  p.list <- lapply(stolist, function(x){
    p1 <- VlnPlot(x, features = "percent.rb") + ggtitle(unique(x$SectionID)) +
      theme(legend.position = "none", axis.text.x = element_blank())
    p2 <- SpatialFeaturePlot(x, features = "percent.rb") + theme(legend.position = "right")
    p <- p1|p2
    p
  })
  p <- wrap_plots(p.list, ncol = 1)
  ggsave("percent_rb_before.pdf", p, width = 10, height = 9)  
  
  ### 质控
  minCount = 1500
  minFeature = 500
  maxmt = 15  ##gene symbol中人线粒体以MT-开头，小鼠为mt-
  stolist <- lapply(stolist, function(x){
    subset(x, nCount_Spatial>minCount&nFeature_Spatial>minFeature&percent.mt<maxmt)
  })
  
  ### 质控后指标
  p.list <- lapply(stolist, function(x){
    p1 <- VlnPlot(x, features = "nCount_Spatial") + ggtitle(unique(x$SectionID)) +
      theme(legend.position = "none", axis.text.x = element_blank())
    p2 <- SpatialFeaturePlot(x, features = "nCount_Spatial") + theme(legend.position = "right")
    p <- p1|p2
    p
  })
  p <- wrap_plots(p.list, ncol = 1)
  ggsave("nCount_Spatial_after.pdf", p, width = 10, height = 9)
  # 基因数
  p.list <- lapply(stolist, function(x){
    p1 <- VlnPlot(x, features = "nFeature_Spatial") + ggtitle(unique(x$SectionID)) +
      theme(legend.position = "none", axis.text.x = element_blank())
    p2 <- SpatialFeaturePlot(x, features = "nFeature_Spatial") + theme(legend.position = "right")
    p <- p1|p2
    p
  })
  p <- wrap_plots(p.list, ncol = 1)
  ggsave("nFeature_Spatial_after.pdf", p, width = 10, height = 9)
  # 线粒体
  p.list <- lapply(stolist, function(x){
    p1 <- VlnPlot(x, features = "percent.mt") + ggtitle(unique(x$SectionID)) +
      theme(legend.position = "none", axis.text.x = element_blank())
    p2 <- SpatialFeaturePlot(x, features = "percent.mt") + theme(legend.position = "right")
    p <- p1|p2
    p
  })
  p <- wrap_plots(p.list, ncol = 1)
  ggsave("percent_mt_after.pdf", p, width = 10, height = 9)
  # 核糖体
  p.list <- lapply(stolist, function(x){
    p1 <- VlnPlot(x, features = "percent.rb") + ggtitle(unique(x$SectionID)) +
      theme(legend.position = "none", axis.text.x = element_blank())
    p2 <- SpatialFeaturePlot(x, features = "percent.rb") + theme(legend.position = "right")
    p <- p1|p2
    p
  })
  p <- wrap_plots(p.list, ncol = 1)
  ggsave("percent_rb_after.pdf", p, width = 10, height = 9)  
  
  
}  

##====================================================================================##
##===========================Step02：聚类分析=======
# setwd("~/spatial/Script/")
### seurat V4 降维聚类,使用sctransform进行数据标准化，并执行标准的 scRNA-seq 降维和聚类工作流
# 
# stolist <- readRDS("D:/R Notebook/10X-mouse ovary-my data ananlysis/QC/2022-10-18-New-QC/stolist18.qc.rds")

####修改orig.ident为CTR和LIP
stolist[[1]]$orig.ident <- "CTR"
stolist[[2]]$orig.ident <- "Lip"



table(stolist$CTR$orig.ident)
# CTR 
# 383 
table(stolist$LIP$orig.ident)
# Lip 
# 384 
stRNA <- merge(stolist[[1]], stolist[[2]])
table(stRNA$orig.ident)
# CTR Lip 
# 383 384 


set.seed(1314)
stRNA <- SCTransform(stRNA, assay = "Spatial", verbose = FALSE)
DefaultAssay(stRNA) <- "SCT"
SpatialFeaturePlot(stRNA, features = c("Cd68", "Inha"))

SpatialFeaturePlot(stRNA, features = c("Enpp2"))

stRNA <- RunPCA(stRNA, assay = "SCT", verbose = FALSE)
ElbowPlot(stRNA, ndims = 50)

stRNA <- FindNeighbors(stRNA, reduction = "pca",dims = 1:30)
stRNA <- FindClusters(stRNA, verbose = FALSE)
stRNA <- RunUMAP(stRNA,  dims = 1:30)
#tsne非线性降维
stRNA <- RunTSNE(stRNA, dims =  1:30) 
DimPlot(stRNA, group.by = "orig.ident", reduction = "umap")
p1 <- DimPlot(stRNA, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(stRNA, label = TRUE, label.size = 3)
p <- p1 + p2
p
# umap图美化

{
  stRNA <- AddMetaData(stRNA,stRNA@reductions$umap@cell.embeddings,col.name = colnames(stRNA@reductions$umap@cell.embeddings))
  
  mytheme <- theme_void()+theme(plot.margin = margin(5.5,15,5.5,5.5))
  source("utools.R")
  mycolor <- c('#66C5CC','#F6CF71','#F89C74','#DCB0F2','#87C55F','#9EB9F3','#FE88B1','#C9DB74','#8BE0A4','#A5AA99')
  mycolor = ov_palette
  p <- ggplot(stRNA@meta.data, aes(x=UMAP_1,y=UMAP_2))+
    geom_point(aes(color=seurat_clusters), size=2)+
    theme_bw(base_size = 16)+mytheme+
    scale_color_manual(values = mycolor)+
    scale_fill_manual(values = mycolor)+
    theme_dr(xlength=0.2,
             ylength=0.2,
             arrow=grid::arrow(length = unit(0.1,"inches"),
                               ends = 'last',
                               type = "closed"))+
    #stat_ellipse(aes(color=seurat_clusters),level = 0.95,
    #linetype=2,show.legend = F)+
    guides(color=guide_legend(override.aes = list(size=4)))
  p <- p+ theme(panel.grid = element_blank())
  p
  ggsave("./Integrate/Cluster_Louvain/seurat_cluster_UMAP.pdf", p, width = 5, height = 5)
}

### 结果展示
p1 <- DimPlot(stRNA, reduction = "umap", label = TRUE, split.by = "SectionID") + plot_layout(guides = "collect")
p1 <-  ggplot(stRNA@meta.data, aes(x=UMAP_1,y=UMAP_2))+
  geom_point(aes(color=orig.ident), size=2)+
  theme_bw(base_size = 16)+mytheme+
  scale_color_manual(values = mycolor)+
  scale_fill_manual(values = mycolor)+
  guides(color=guide_legend(override.aes = list(size=4)))+
  theme_dr(xlength=0.2,
           ylength=0.2,
           arrow=grid::arrow(length = unit(0.1,"inches"),
                             ends = 'last',
                             type = "closed"))
p1 <- p1+ theme(panel.grid = element_blank())
p1
ggsave("./Integrate/Cluster_Louvain/seurat_cluster_orig.ident.pdf", p1, width = 5, height = 5)
p2 <- SpatialDimPlot(stRNA, label = F, repel=F, label.size = 2) + plot_layout(guides = "collect")
p <- p1/p2
p
ggsave("seurat_cluster_SpatialDimPlot_withoutlabel.pdf", p2, width = 9, height = 9)


save(stRNA,file = './Integrate/stRNA_orig.Rdata')


#################################################################################
# begin
rm(list = ls())
gc()
load("./Integrate/stRNA_orig.Rdata")
####高变基因的选择
stRNA <- FindSpatiallyVariableFeatures(stRNA, assay = "SCT", features = VariableFeatures(stRNA)[1:1000], 
                                       selection.method = "markvariogram")


## 保存临时结果
save(stRNA, file = "stRNA.tmp.rda")


SpatialDimPlot(stRNA)





###可视化
{
  #用TSNEPlot函数绘制tsne图
  # tsneplot<-TSNEPlot(stRNA,label = TRUE, pt.size = 1.5)
  # tsneplot
  # 
  # #绘制 Marker 基因的 tsne 图； 
  # FeaturePlot(stRNA, reduction = "tsne", features = c("Enpp2", "Gdf9","Des","Lyz2","Inha","Inhbb")) 
  # 
  # if (!require("BiocManager", quietly = TRUE))
  #   install.packages("BiocManager")
  # 
  # BiocManager::install("BiocFileCache",force = TRUE)
  
  
  library("Nebulosa")
  library("Seurat")
  library("BiocFileCache")
  
  
  # FeaturePlot(stRNA, reduction = "tsne",features=c("Gdf9", "Amhr2","Tcf21","Des","Cdh5","Krt19","Ptprc","Ptgfr"),ncol = 4 )
  
  #绘制 Marker 基因的星云 图；
  
  p <-plot_density(stRNA, c("Gdf9", "Amhr2","Tcf21","Des","Cdh5","Krt19","Ptprc","Ptgfr"),
                   pal = "magma") 
  
  p1 <- p+plot_layout(ncol = 4)
  p1
  ggsave("./Integrate/marker/markergenes-plot_density_tsne.pdf", p1, width = 17, height =6)
  
  
  p <-plot_density(stRNA, c("Enpp2", "Cd68"), joint = TRUE)
  p <- p+plot_layout(ncol = 3)
  ggsave("./Integrate/marker/Enpp2_Cd68-plot_density_tsne.pdf", p, width =15.5, height =4)
  
  # Return joint density plot
  p <-plot_density(stRNA, c("Gpnmb", "Cd68"),pal = "magma", joint = TRUE)
  # +scale_colour_gradientn(colours =  c("grey90","#f1707d"))
  p <- p+plot_layout(ncol = 3)
  p
  ggsave("./Integrate/marker/Gpnmb_Cd68-plot_density_tsne.pdf", p, width =12, height =4)
  
  #绘制 Marker 基因的 SpatialPlot 图；
  p <-SpatialPlot(stRNA,features="Gpnmb",
                  image.alpha = 1,
                  crop = TRUE,
                  # cols.highlight = c("#DE2D26", "grey50"),
                  label.size = 5,
                  pt.size.factor = 1.6,
                  alpha = c(1, 1),
                  stroke = 0.25,
                  label.box = TRUE)+scale_colour_gradientn(colours =  c("grey90","#f1707d"))
  
  ggsave("./Integrate/marker/Gpnmb-SpatialPlot_merged.pdf", p, width = 7, height = 4)
  
  #绘制 Marker 基因的 SpatialPlot 图；
  p <-SpatialPlot(stRNA,features="Cd68",
                  image.alpha = 1,
                  crop = TRUE,
                  # cols.highlight = c("#DE2D26", "grey50"),
                  label.size = 5,
                  pt.size.factor = 1.6,
                  alpha = c(1, 1),
                  stroke = 0.25,
                  label.box = TRUE)+scale_colour_gradientn(colours =  c("grey90","#f1707d"))
  
  ggsave("./Integrate/marker/Cd68-SpatialPlot_merged.pdf", p, width = 7, height = 4)
  
  
  features = c("Enpp2", "Gdf9","Amhr2","Des","Lyz2","Inha","Inhbb","Cyp19a1","Cyp11a1","Cd68","Wt1","Slc38a3","Wnt6","Lpar6","Prss23",'Hsd17b1','Ptgfr','Sfrp4','Lum','Cyp17a1','Enpep','Aldh1a2','Tcf21','Star','Ptgs2','Plin2')
  #绘制 Marker 基因的 SpatialFeaturePlot 图；        
  p <- SpatialFeaturePlot(stRNA,features="Enpp2", alpha = c(0.1, 1)) 
  ggsave("Enpp2-SpatialFeaturePlot_merged.pdf", p, width = 7, height = 4.5)
  
  
  
  #降维可视化
  p1 <- DimPlot(stRNA, reduction = "tsne", group.by = "orig.ident", pt.size = 1.5)
  p2 <- DimPlot(stRNA, reduction = "tsne", label = TRUE, repel = TRUE,pt.size = 1.5)
  p2 
  #(1) 亚群间差异基因分析,同时分析所有亚群的marker基因
  dif<-FindAllMarkers(stRNA,only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25,test.use = "wilcox")
  write.csv(dif,"dif between ST0-9.csv")
  
  #logfc.threshold定义上调倍数阈值，min.pct定义基因至少在细胞亚群中多少细胞中表达，only.pos确定只筛选上调基因
  #通过dplyr包来完成对top基因的筛选
  sig.dif<-dif%>%group_by(cluster)%>%top_n(n = 6,wt = avg_log2FC)
  #保存差异分析结果
  write.table(sig.dif,"sig.dif.xls",row.names = T,col.names = T,quote = F,sep = "\t")
  
  ##查看每个亚群top5差异基因
  topgene<-dif%>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  topgene
  
  ##单独两个或多个亚群比较分析marker基因分析
  markers <- FindMarkers(stRNA, ident.1 = 3, ident.2 = 9)
  
  ##绘制每个亚群top5差异基因热图
  {
    p6 <- DoHeatmap(stRNA, features = topgene$gene,size = 2) + NoLegend()
    
    ##美化热图
    library(RColorBrewer)
    mypalette <- brewer.pal(n = 9, name = "YlOrRd")
    p6 <- DoHeatmap(stRNA,
                    features = as.character(unique(topgene$gene)),
                    group.by = "seurat_clusters",
                    group.colors =mycolor)+
      scale_fill_gradientn(colors = c("white","grey90","firebrick3"))
    
    ggsave("Top5DEG_seurat_cluster_heatmap.pdf", p6, width = 12, height = 9)
  }
  
  ##做不同ST cluster的火山图
  {
    
    dif$label <- ifelse(dif$p_val<0.05,"adjust P-val<0.05","adjust P-val>=0.05") 
    filtered_dif <- dif[dif$p_val < 0.05&abs(dif$avg_log2FC) >1,]
    dim(filtered_dif)
    #[1] 3082    8
    
    
    p <- ggplot()+
      geom_jitter(data = dif,
                  aes(x = cluster, y = avg_log2FC,color=label),
                  size = 0.5,
                  width =0.4)+
      geom_jitter(data = topgene,
                  aes(x = cluster, y = avg_log2FC),
                  size = 1,
                  width =0.4)
    
    #根据图p中log2FC区间确定背景柱长度：
    dfbar<-data.frame(x=c("0","1","2","3","4","5","6","7","8","9"),
                      y=c(1.8,5.3,5,4.2,2.5,2.1,2.1,5,5,4))
    dfbar1<-data.frame(x=c("0","1","2","3","4","5","6","7","8","9"),
                       y=c(-3,-2.8,-3.1,-3.2,-3,-2,-2.5,-3.2,-3,-1.5))
    #绘制背景柱：
    p1 <- ggplot()+
      geom_col(data = dfbar,
               mapping = aes(x = x,y = y),
               fill = "#dcdcdc",alpha = 0.5)+
      geom_col(data = dfbar1,
               mapping = aes(x = x,y = y),
               fill = "#dcdcdc",alpha = 0.5)
    
    #添加X轴的cluster色块标签：
    dfcol<-data.frame(x=c("0","1","2","3","4","5","6","7","8","9"),
                      y=0,
                      label=c("0","1","2","3","4","5","6","7","8","9"))
    
    p <- ggplot()+
      geom_col(data = dfbar,
               mapping = aes(x = x,y = y),
               fill = "#dcdcdc",alpha = 0.5)+
      geom_col(data = dfbar1,
               mapping = aes(x = x,y = y),
               fill = "#dcdcdc",alpha = 0.5)+
      geom_jitter(data = dif,
                  aes(x = cluster, y = avg_log2FC,color=label),
                  size = 0.4,
                  width =0.4)+
      geom_jitter(data = topgene,
                  aes(x = cluster, y = avg_log2FC),
                  size = 0.8,
                  width =0.4)+
      geom_tile(data = dfcol,
                aes(x=x,y=y),
                height=0.8,
                color = "black",
                fill = mycolor,
                alpha = 0.6,
                show.legend = F)
    
    library(ggrepel)
    #给每个Cluster差异表达前Top5基因加上标签：
    p2 <- p + geom_tile(data = dfcol,
                        aes(x=x,y=y),
                        height=0.8,
                        color = "black",
                        fill = mycolor,
                        alpha = 0.6,
                        show.legend = F)+
      geom_text_repel(
        data=topgene,
        aes(x=cluster,y=avg_log2FC,label=gene),
        size =3,
        arrow = arrow(length = unit(0.008, "npc"),
                      type = "open", ends = "last")
      )
    
    
    p2 <- p2 +
      scale_color_manual(name=NULL,
                         values = c("#f8766d","grey20"))
    
    
    
    
    
    p2 <- p2+
      labs(x="Cluster",y="avg_log2FC")+
      geom_text(data=dfcol,
                aes(x=x,y=y,label=label),
                size =3,
                color ="white")
    
    
    
    p2 <- p2 +
      geom_text(data=dfcol,
                aes(x=x,y=y,label=label),
                size =4,
                color ="white")+
      theme_minimal()+
      theme(
        axis.title = element_text(size = 12,
                                  color = "black",
                                  face = "bold"),
        axis.line.y = element_line(color = "black",
                                   size = 0.8),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.direction = "vertical",
        legend.justification = c(1,0),
        legend.text = element_text(size = 12)
      )
    
    ggsave("Volcano_DEG_Cluster0_9.pdf",p2,width = 12,height = 8)
  }
  
  ##dotplot
  genelist=c('Tcf21','Cd34','Des','Inhba','Inhbb','Hsd17b1','Ptgfr','Sfrp4','Aldh1a1','Mgarp','Adh1','Grb14','Prss23','Ddx4','Gdf9','Lyz2','Cd68','Adgre1','Epcam','Krt19')
  p <- DotPlot(stRNA,features = genelist)
  ggsave("0-9-markergenes-DotPlot.pdf", p, width = 14, height = 9)
  
  p <-  FeaturePlot(stRNA, reduction = "tsne", features ='Col1a1') 
  p <- p+mytheme
  p
  
  rm(p)
  
  
  
  
  
  
  ###### 使用marker基因确定细胞分布,已完成，保存在8个主成分那里
  {
    DefaultAssay(stRNA) <- "SCT"
    markerlist <- list(
      Stromal_cells = c("Tcf21","Dcn"),
      GC = c("Inha", "Inhba"),
      SMA = c("Des","Col1a1"),
      Luteal_cells= c("Ptgfr","Sfrp4"),
      Theca_cells = c("Aldh1a1","Star"),
      Cumulus_cells = c("Grb14", "Prss23"),
      Epithelial_cells = c("Epcam","Krt19"),
      Endothelial_cells = c("Cd34","Vwf"),
      Oocytes = c("Ddx4","Gdf9","Bmp15","Dazl"),
      Immune_cells = c("Lyz2", "Cd52","Cd68")
    )
    
    for(i in names(markerlist)){
      markers <- markerlist[[i]]
      p <- SpatialFeaturePlot(stRNA, features = markers) + plot_layout()&theme(legend.position = "right")
      ggsave(paste0(i, ".pdf"), p, width = 10, height = 8)
    }
    {
      markers <- do.call("c",markerlist)
      p <- SpatialFeaturePlot(stRNA, features = markers,images = "CTR", ncol = 4) + plot_layout()&theme(legend.position = "right")
      ggsave("Markers_CTR1.pdf", p, width = 18, height = 30)
      
      p <- SpatialFeaturePlot(stRNA, features = markers, images = "LIP", ncol = 4) + plot_layout()&theme(legend.position = "right")
      ggsave("Markers_LIP2.pdf", p, width = 18, height = 30)
    }
  }
}


###细胞堆叠图
{
  table(stRNA$orig.ident)
  # CTR Lip 
  # 379 382
  
  
  prop.table(table(Idents(stRNA)))
  table(Idents(stRNA), stRNA$orig.ident)
  
  
  Cellratio <- prop.table(table(Idents(stRNA), stRNA$SectionID), margin = 2)#计算各组样本不同细胞群比例
  Cellratio <- as.data.frame(Cellratio)
  write.csv(Cellratio,"Cellration_stRNA.csv")
  
  allcolour=c("#00468B","#925E9F","#759EDD","#0099B4","#76D1B1","#42B540","#B8D24D",
              "#EDE447","#FAB158","#FF7777","#FD0000","#AD002A","#AE8691","#CE9573",
              "#756455")
  library(ggplot2)
  p6 <- ggplot(Cellratio) + 
    geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.6,size = 0.5,colour = '#222222')+ 
    theme_classic() +
    labs(x='Sample',y = 'Ratio')+
    scale_fill_manual(values = mycolor)+
    theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
  p6
  ggsave("细胞堆叠图_CTR1_LIP2.pdf", p6, width = 8, height = 6)
  
  
  p <- VlnPlot(stRNA,features=c("nFeature_Spatial","nCount_Spatial","percent.mt"),ncol=3)
  p
  ggsave("总的Lip组—nFeature—nCount-percent.mt.pdf", p, width = 14, height = 7)
  plot1 <- FeatureScatter(stRNA,feature1 = "nCount_Spatial",feature2 = "percent.mt")
  plot2 <- FeatureScatter(stRNA,feature1 = "nCount_Spatial",feature2 = "nFeature_Spatial")
  p3 <- plot1/plot2
  ggsave("Lip组—nFeature—nCount-percent.mt-散点图.pdf", p3, width = 7, height = 14)
  
  
  
}




##====================================================================================##
##===========================Step03：去卷积分析=======================================##
##====================================================================================##
library(Seurat)
library(tidyverse)
library(patchwork)

rm(list = ls())
# stRNA <- readRDS("D:/R Notebook/10X-mouse ovary-my data ananlysis/QC/2022-10-18-New-QC/sto_cluster18.rds")
load("~/spatial/processed_data/Integrate/QC/stRNA_orig.Rdata")


##### SPOTlight去卷积##### 
#跨平台能力比较差，比较适合配对的10X单细胞和空转

if(F) {
  library(SPOTlight)
  library(NMF)
  library(ggplot2)
  library(scater)
  library(Seurat)
  library(SingleCellExperiment)
  library(scater)
  library(scran)
  library(Matrix)
  setwd("/home/data/t070429/spatial/processed data/Integrate/Deconvolution/SPOTlight")
  
  
  #SPOTLight安装方式
  {if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
    BiocManager::install("SingleCellExperiment")
    
    
    install.packages("BiocManager")
    BiocManager::install("SPOTlight")
    # Or the devel version
    BiocManager::install("SPOTlight", version = "devel")
    
    
    
    
  }
  ### 去卷积
  Idents(SC1) <- "Level0"
  ct.marker <- FindAllMarkers(SC1, logfc.threshold = 0.5)
  save(ct.marker, file = "celltype_marker.rda")
  set.seed(123)
  
  
  
  
  scRNA <- as.SingleCellExperiment(SC1)
  
  
  ######下面是一些官网的核心步骤
  
  scRNA <- logNormCounts(scRNA)
  
  
  dec <- modelGeneVar(scRNA)
  plot(dec$mean, dec$total, xlab = "Mean log-expression", ylab = "Variance")
  curve(metadata(dec)$trend(x), col = "blue", add = TRUE)
  #####计算高变基因，官方文档建议3000个
  hvg <- getTopHVGs(dec,n=3000)
  
  saveRDS(hvg, file = "spotlight_hvg3000.rds")
  options(stringsAsFactors = FALSE)
  
  #加上细胞注释信息
  colLabels(scRNA) <- SC1$Level0
  #去掉核糖体和线粒体基因
  genes <- !grepl(pattern = "^Rp[l|s]|mt",x=rownames(scRNA))
  #计算并保留最显著的marker基因
  mgs <- scoreMarkers(scRNA,subset.row=genes)
  #保留最相关的marker基因
  mgs_fil <- lapply(names(mgs),function(i){
    x <- mgs[[i]]
    # Filter and keep relevant marker genes, those with AUC > 0.8
    x <- x[x$mean.AUC > 0.7, ]
    # Sort the genes from highest to lowest weight
    x <- x[order(x$mean.AUC, decreasing = TRUE), ]
    # Add gene and cluster id to the dataframe
    x$gene <- rownames(x)
    x$cluster <- i
    data.frame(x)
  })
  mgs_df <- do.call(rbind, mgs_fil)
  
  ######每种细胞类型提取100个，加快计算时间，一般来说细胞数母越多越好
  
  ####cell downsampling
  # split cell indices by identity
  idx <- split(seq(ncol(scRNA)), scRNA$Level0)
  # downsample to at most 20 per identity & subset
  # We are using 5 here to speed up the process but set to 75-100 for your real
  # life analysis
  n_cells <- 100
  cs_keep <- lapply(idx, function(i) {
    n <- length(i)
    if (n < n_cells)
      n_cells <- n
    sample(i, n_cells)
  })
  sce <- scRNA[, unlist(cs_keep)]
  
  ###去卷积分析 
  res <- SPOTlight(
    x = as.matrix(sce@assays@data@listData$counts),
    y = as.matrix(stRNA@assays$Spatial@counts),
    groups = as.character(sce$Level0),
    mgs = mgs_df,
    hvg = hvg,
    weight_id = "mean.AUC",
    group_id = "cluster",
    gene_id = "gene")
  
  
  saveRDS(res, file = "spotlight_res.rds")
  
  ### Extract deconvolution matrix
  head(mat <- res$mat)[,seq_len(3)]
  
  # Extract NMF model fit
  mod <- res$NMF
  
  ###结果可视化
  ###Topic profiles
  plotTopicProfiles(
    x = mod,
    y = sce$Level0,
    facet = FALSE,
    min_prop = 0.01,
    ncol = 1) +
    theme(aspect.ratio = 1)
  
  library(NMF)
  
  ## 检查主题对应的基因
  sign <- basis(mod)
  colnames(sign) <- paste0("Topic", seq_len(ncol(sign)))
  write.csv(sign, file = "spotlight_Topic_signautre.csv")
  head(sign)
  
  ###结果可视化
  stRNA[["SPOTlight"]] <- CreateAssayObject(t(res$mat))
  DefaultAssay(stRNA) <- "SPOTlight"
  celltypes= rownames(stRNA)
  p <- SpatialFeaturePlot(stRNA,features = celltypes, pt.size.factor=1.6,ncol=4, crop=TRUE)
  ggsave("SPOTlight_Deconvolution_test1.pdf",p,width=18,height=32, limitsize=F)
  
  #细胞相关性
  library(ggcorrplot)
  
  p <- plotCorrelationMatrix(res$mat)
  ggsave("SPOTlight_plotCorrelationMatrix_test1.pdf",p,width=8,height=6)
  
  #细胞共定位
  p <-plotInteractions(res$mat,which = "heatmap", metric = "prop")
  ggsave("SPOTlight_plotInteractions_prop_test1.pdf",p,width=8,height=6)
  p <- plotInteractions(res$mat,which = "heatmap", metric = "jaccard")
  ggsave("SPOTlight_plotInteractions_jaccard_test1.pdf",p,width=8,height=6)
  p <-plotInteractions(res$mat,which = "network")
  ggsave("SPOTlight_plotInteractions_test1.pdf",p,width=8,height=6)
  
  #spot成分饼图
  library(ggsci)
  library(scatterpie)
  
  library(SpatialExperiment)
  # Get spatial coordinates
  
  coords2 <- stRNA@images[["LIP2"]]@coordinates
  coords2<- coords2[,c(4:5)]
  coords1 <- stRNA@images[["CTR1"]]@coordinates
  coords1<- coords1[,c(4:5)]
  
  
  
  
  ##获取切片图像坐标和key，这一步花了我很多时间，还是基础太薄弱了
  {Key(object = stRNA@images$CTR1)
    img1 <- GetTissueCoordinates(stRNA)
    Key(object = stRNA@images$LIP2)
    img2 <- GetTissueCoordinates(stRNA)
    
    head(img)
  }
  
  
  mat <- res$mat
  mat[mat< 0.1] <- 0
  mat_matrix <- as.matrix(mat)
  #出现了Error in plotSpatialScatterpie(x = stRNA, y = mat_matrix, cell_types = colnames(mat_matrix),  : 
  #nrow(x) == nrow(y) is not TRUE报错时，先检查是否对齐
  nrow(merged_df) == nrow(df)
  #[1] FALSE
  ncol(merged_df) == ncol(df)
  
  df <- as.data.frame(mat)
  df1 <- df[1:379, ]
  df2 <- df[380:761, ]
  p <- plotSpatialScatterpie(x=coords2, y=df2, cell_types= colnames(df2),img=FALSE,scatter_alpha=1, pie_scale=0.4)
  p
  
  
  ggsave("SPOTlight_composition_LIP2.pdf",p,width=9,height=6, limitsize=F)
  
  
  
  
  
  #安装下包
  {if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
    BiocManager::install("SpatialExperiment")
    library(SpatialExperiment)
  }
  
  ## 保存结果
  
  save(ct.marker, df1, df2, coords1,coords2,sce,res,mgs_fil, hvg, mod,scRNA, stRNA, file = "spotlight_data.rda")
}

##### RCTD去卷积######  
if(F){
  
  # install.packages("devtools")
  options(timeout = 600000000) ### set this to avoid timeout error
  devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
  library(ggplot2)
  library(spacexr)
  library(Matrix)
  library(doParallel) 
  #设置工作路径
  setwd("/home/data/t070429/spatial/processed data/Integrate/Deconvolution/RCTD")
  
  ### 创建单细胞分析对象
  
  
  #用elife的单细胞数据提取counts时一直报错，要不是counts矩阵不是integer，要不稀疏矩阵转换为dataframe的时候，提示内存耗尽
  
  {
    rm(sce, i,idx,df1,df2,cs_keep,merged_df1,mgs_fil)
    counts = SC1@assays$RNA@counts
    cell_types = factor(structure(SC1$Level0, names=colnames(SC1)))
    nUMI = structure(SC1$nCount_RNA, names=colnames(SC1))
    reference <- Reference(counts, cell_types, nUMI)
    #Error in check_counts(counts, "Reference", require_2d = T, require_int = require_int) : 
    # Reference: counts does not contain integers
  }
  rm(counts,cell_types,nUMI) #移除之前的一些变量
  
  #所以后面用MCA的那个reference尝试了一下。 
  reference <- readRDS("/home/data/t070429/MCA-ovary/RCTD_MCA_reference.rds")
  
  sto.list <- SplitObject(stRNA, split.by = "SectionID")
  sto.CTR1 <- sto.list[["CTR1"]]
  sto.LIP2 <- sto.list[["LIP2"]]
  
  
  query.CTR1 <- SpatialRNA(coords = coords1, counts = sto.CTR1@assays$Spatial@counts, 
                           nUMI = structure(sto.CTR1$nCount_Spatial, names=colnames(sto.CTR1)))
  RCTD.CTR1 <- create.RCTD(spatialRNA = query.CTR1, reference = reference, max_cores = 10)
  
  
  query.LIP2 <- SpatialRNA(coords = coords2, counts = sto.LIP2@assays$Spatial@counts, 
                           nUMI = structure(sto.LIP2$nCount_Spatial, names=colnames(sto.LIP2)))
  RCTD.LIP2  <- create.RCTD(spatialRNA = query.LIP2, reference = reference, max_cores = 10)
  
  ### 去卷积
  # doublet对应每个spot 1-2个细胞，multi mode对应3-4个细胞，full mode没有限制
  RCTD.CTR1 <- run.RCTD(RCTD.CTR1, doublet_mode = 'full')
  RCTD.LIP2  <- run.RCTD(RCTD.LIP2 , doublet_mode = 'full')
  save(RCTD.CTR1, RCTD.LIP2 , file = "ovaryspatial_RCTD.rda")
  
  ### 可视化
  ## 细胞类型预测结果
  pred.CTR1 <- as.matrix(t((normalize_weights(RCTD.CTR1@results$weights))))
  pred.LIP2 <- as.matrix(t((normalize_weights(RCTD.LIP2 @results$weights))))
  pred.mat <- cbind(pred.CTR1, pred.LIP2)
  stRNA[["RCTD"]] <- CreateAssayObject(pred.mat[,colnames(stRNA)])
  DefaultAssay(stRNA) <- "RCTD"
  library(ggplot2)
  library(patchwork)
  celltypes = rownames(stRNA)
  p <- SpatialFeaturePlot(stRNA, features = celltypes, images = "CTR1", pt.size.factor = 1.6, combine = F, crop = TRUE)
  p <- lapply(p, function(x){x+labs(title = x$labels$fill)+theme(legend.position="right", legend.title=element_blank())})
  p <- wrap_plots(p, ncol = 3)
  ggsave("RCTD_celltype_predicted_CTR1.pdf", p, width = 12, height = 24)
  
  p <- SpatialFeaturePlot(stRNA, features = celltypes, images = "LIP2", pt.size.factor = 1.6, combine = F, crop = TRUE)
  p <- lapply(p, function(x){x+labs(title = x$labels$fill)+theme(legend.position="right", legend.title=element_blank())})
  p <- wrap_plots(p, ncol = 3)
  ggsave("RCTD_celltype_predicted_LIP2 .pdf", p, width = 12, height = 24)
  
  ## spot成分饼图
  pred.mat[pred.mat<0.1] <- 0
  pred.df <- t(pred.mat)
  pred.df <- data.frame(pred.df)
  colnames(stRNA@meta.data)
  stRNA@meta.data <- stRNA@meta.data[,-c(18:34)]
  stRNA <- AddMetaData(stRNA, metadata = pred.df)
  ct <- colnames(pred.df)
  library(STdeconvolve)
  library(ggsci)
  packageVersion("STdeconvolve")
  
  weights <- RCTD.LIP2@results$weights
  norm_weights <- normalize_weights(weights)
  m <- as.matrix(norm_weights)
  colnames(coords2) <- c("x", "y")
  my_palette <- colorRampPalette(c("#F0E68C", "#90EE90", "#00BFFF"))
  plt <- vizAllTopics(theta=m,
                      pos=coords2,
                      topicOrder=seq(ncol(m)),
                      topicCols=my_palette(ncol(m)),
                      groups = NA,
                      group_cols = NA,
                      r = 39,
                      lwd = 0.2,
                      showLegend = TRUE,
                      overlay = NA
                      #plotTitle = "scatterpies"
  )  #CTR的r设置为3比较合适
  
  plt
  plt <- plt+ggplot2::guides(fill=ggplot2::guide_legend(ncol = 2))
  ggsave("RCTD_SpatialScatterpie_LIP2.pdf", width = 12,height = 8, plot = plt,bg="white")
  
  ##### 输出去卷积矩阵
  decon.mat <- as.matrix(stRNA@assays$RCTD@data)
  decon.mat <- t(decon.mat)
  decon.mat <- decon.mat[grepl("^LIP2", rownames(decon.mat)),]
  rownames(decon.mat) <- gsub("LIP2_", "", rownames(decon.mat))
  decon.mat <- as.data.frame(decon.mat)
  decon.mat$`Macrophage` <- 0.4*(decon.mat$`Macrophage`)
  decon.mat <- data.frame(decon.mat, celltype=colnames(decon.mat)[apply(decon.mat,1,which.max)],check.names = F)
  table(decon.mat$celltype)
  write.csv(decon.mat, "RCTD_IP2.csv")
  
  rm(decon.mat,img1,img2,mat,mgs,mod,pred.CTR1,pred.LIP2,pred.mat,Q_mat,query.CTR1,query.LIP2,RCTD.CTR1,RCTD.LIP2,S_mat,SQ_mat,sto.CTR1,sto.LIP2,stolist,K_val,labers,N_X,X_vals)
  
  
}



##===========================Step04：空间变异基因===================================
setwd("/home/data/t070429/spatial/processed data/Integrate/Deconvolution/SVGs")
lapply(c("DomainDEG","Seurat","SPARK","SpatialDE"), FUN = dir.create)

##### DomainDEG 区域差异分析检测SVGs 

{
  if(F){
    library(Seurat)
    library(tidyverse)
    library(patchwork)
  }
  
  colnames(stRNA@meta.data)
  DefaultAssay(stRNA) <- "Spatial"
  p <- SpatialDimPlot(stRNA, group.by = "seurat_clusters", label = T, label.size = 5) + NoLegend()
  ggsave("SpatialSeurat_clusters.pdf", p, width = 15, height = 8)
  
  ### 差异分析
  Idents(stRNA) <- "seurat_clusters"
  domainDEGs <- FindAllMarkers(stRNA, only.pos = T)
  saveRDS(domainDEGs, "domainDEGs.rds")
  
  ss_markers <- group_by(domainDEGs, cluster) %>% top_n(1, avg_log2FC) %>% pull(gene)
  p <- SpatialFeaturePlot(stRNA, features = ss_markers, images = "CTR1", alpha = c(0.1, 1), ncol = 4)
  p <- p + plot_layout()&theme(legend.position = "right")
  ggsave("svg_seuratcluster_deg_top1_CTR1.pdf", p, width = 18, height = 16)
  
  ss_markers <- group_by(domainDEGs, cluster) %>% top_n(1, avg_log2FC) %>% pull(gene)
  p <- SpatialFeaturePlot(stRNA, features = ss_markers, images = "LIP2", alpha = c(0.1, 1), ncol = 4)
  p <- p + plot_layout()&theme(legend.position = "right")
  ggsave("svg_seuratcluster_deg_top1_LIP2.pdf", p, width = 18, height = 16)
}
###### Seurat方法检测SVGs
{
  setwd("~/spatial/processed data/Integrate/Deconvolution/SVGs/Seurat")
  set.seed(235)
  DefaultAssay(stRNA) <- "SCT"
  features <- sample(VariableFeatures(stRNA), 500)
  
  # markvariogram(标记变异函数)受Trendsceek的启发
  T0 <- Sys.time()
  stRNA <- FindSpatiallyVariableFeatures(stRNA, assay = "SCT", 
                                         features = features,
                                         #features = VariableFeatures(stRNA),
                                         selection.method = "markvariogram")
  T1 <- Sys.time() - T0  #Time difference of 2.600504 mins
  
  # moransi(莫兰指数)
  T0 <- Sys.time()
  stRNA <- FindSpatiallyVariableFeatures(stRNA, assay = "SCT", 
                                         features = features,
                                         #features = VariableFeatures(stRNA),
                                         selection.method = "moransi")
  T2 <- Sys.time() - T0  # Time difference of 41.07003 secs
  
  svgs <- na.omit(stRNA@assays$SCT@meta.features)
  saveRDS(svgs, "svgs_markvario_morans.rds")
  
  ### 筛选SVGs可视化
  svgs.sub <- svgs[order(svgs$markvariogram.spatially.variable.rank),] %>% rownames() %>% head(9)
  p <- SpatialFeaturePlot(object = stRNA, features = svgs.sub, images = "CTR1", alpha = c(0.1, 1), ncol = 3)
  p <- p + plot_layout()&theme(legend.position = "right")
  ggsave("svg_markvariogram_top9_CTR1.pdf", p, width = 12, height = 9)
  p <- SpatialFeaturePlot(object = stRNA, features = svgs.sub, images = "LIP2", alpha = c(0.1, 1), ncol = 3)
  p <- p + plot_layout()&theme(legend.position = "right")
  ggsave("svg_markvariogram_top9_LIP2.pdf", p, width = 12, height = 9)
  
  
  #莫兰指数筛选的TOP9可以筛选到ENPP2
  svgs.sub <- svgs[order(svgs$moransi.spatially.variable.rank),] %>% rownames() %>% head(9)
  p <- SpatialFeaturePlot(object = stRNA, features = svgs.sub, images = "CTR1", alpha = c(0.1, 1), ncol = 3)
  p <- p + plot_layout()&theme(legend.position = "right")
  ggsave("svg_moransi_top9_CTR1.pdf", p, width = 12, height = 9)
  p <- SpatialFeaturePlot(object = stRNA, features = svgs.sub, images = "LIP2", alpha = c(0.1, 1), ncol = 3)
  p <- p + plot_layout()&theme(legend.position = "right")
  ggsave("svg_moransi_top9_LIP2.pdf", p, width = 12, height = 9)
  
}

##===========================Step05：功能富集分析===================================
#准备工作，加载包，设置路径和文件夹
{
  library(future)
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("GSVA")
  
  
  library(GSVA)
  library(AUCell)
  
  library(msigdbr)
  
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("clusterProfiler")
  library(clusterProfiler)
  library(tidyverse)
  library(patchwork)
  
  setwd("~/spatial/processed data/Integrate/Function")
  
  lapply(c("AUCell","GSVA","GSEA_KEGG","GSEA_HALLMARK","ORA_GO","ORA_KEGG"), FUN=dir.create)
  
}
### 组间差异基因
stRNA <- readRDS("~/project/Data/sto_cluster.rds")
colnames(stRNA@meta.data)
stRNA.list <- SplitObject(stRNA, split.by = "seurat_clusters")
#这里assay用的spatial的话需要再进行一步标准化，因为用的原始数据，但是SCT转化后的已经是标准化后的，scale.data可以看出来
#所以后面我用的SCT的，没有进行标准化。回头可以看看用spatial再进行normalization富集出来的结果会不会不一样
stRNA.list <- lapply(stRNA.list, function(sco){
  DefaultAssay(sco) <- "Spatial"
  sco <- NormalizeData(sco)
  return(sco)
})
# 检测细胞数量
table(stRNA$seurat_clusters, stRNA$SectionID)
# 剔除不分析的细胞
#stRNA.list <- stRNA.list[-which(names(stRNA.list)=="15")]  
## 差异分析
T0 <- Sys.time()
deg.list <- list()
for(i in names(stRNA.list)){
  plan(multisession, workers = 15) 
  deg.list[[i]] <- FindMarkers(stRNA.list[[i]], ident.1 = "CTR1", ident.2 = "LIP2", #test.use = "MAST",
                               assay = "SCT", group.by = "SectionID", logfc.threshold = 0, min.pct = 0)
  plan(sequential) 
}
Sys.time()-T0  #Time difference of 2.815803 hours
deg.list <- deg.list[as.character(0:9)]
saveRDS(deg.list, file = "deg.list.group.comparision.rds") 


### cluster特异性基因
T0 <- Sys.time()
deg.list <- list()
for(i in sort(unique(stRNA$seurat_clusters))){
  plan(multisession, workers = 15) 
  deg.list[[i]] <- FindMarkers(stRNA, ident.1 = i, ident.2 = NULL, #test.use = "MAST",
                               assay = "SCT", group.by = "seurat_clusters", logfc.threshold = 0, min.pct = 0)
  plan(sequential) 
}
#ident.2是对照组，ident.1与ident.2相比的差异基因
Sys.time()-T0  #Time difference of 28.24903 mins
names(deg.list) <- as.character(0:9)

saveRDS(deg.list, file = "deg.list.cluster.specific.rds")

###### GO与KEGG分析
# 调整名称
names(deg.list) <- paste0("cluster", 0:9)
# 筛选基因
gene.list <- lapply(deg.list, function(x){
  subset(x, p_val_adj < 0.05&avg_log2FC>0.5) %>% row.names()
})
options(stringsAsFactors = F)
Sys.setenv(R_MAX_NUM_DLLS=999)

### GO分析
{
  library(AnnotationDbi)
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("org.Mm.eg.db")
  library(org.Mm.eg.db)#基因注释包
  
  library(clusterProfiler)#富集包
  library(dplyr)
  library(ggplot2)#画图包
  source("/home/data/t070429/spatial/processed data/Integrate/Function/st_function .txt")
  
  for(i in names(gene.list)){
    genelist <- gene.list[[i]]
    if(length(genelist)>=10) {
      #runGO(genelist, output="ORA_GO", name = paste0(i ,"_TopMarkers"))
      runGO(genelist, output="ORA_GO", orgdb='org.Mm.eg.db', name = paste0(i ,"_TopMarkers"))
    }
  }
}

### KEGG分析

# 网络有问题时尝试调整下载设置
# R.utils::setOption('clusterProfiler.download.method', 'auto')
for(i in names(gene.list)){
  genelist <- gene.list[[i]]
  if(length(genelist)>=10) {
    runKEGG(genelist = genelist, output = "ORA_KEGG", name = paste0(i ,"_TopMarkers"))
  }
} 

### hallmark数据集
###### GSEA

#设置工作路径

deg.list <- readRDS("Function/deg.list.cluster.specific.rds")
deg.list <- lapply(deg.list, function(x){subset(x, pct.1>0|pct.2>0)})
# 调整名称,filenames顺序一定要和deg.list一致，否则会报错
names(deg.list) <- paste0("cluster", 0:9)
filename = names(deg.list)

# GSEA分析
{
  #msigdbr提供多个物种的基因集数据
  # 设置基因集
  view(msigdbr_collections()) #查看msigdbr包中所有的基因集，species能查看物种
  view(msigdbr_species())
  
  genesets = msigdbr(species = "Mus musculus", category = "C2")
  unique(genesets$gs_subcat)  # 有多个数据库来源的基因集可选
  genesets <- subset(genesets, gs_subcat=="CP:KEGG", select = c("gs_name", "gene_symbol"))
  unique(genesets$gs_name) #查看有多少条通路（）
  
  setwd("~/spatial/processed data/Integrate/Function")
  #按差异倍数降序并返回GSEA结果
  res.list <- lapply(deg.list, FUN = function(x){
    # x = deg.list[[1]]
    x = x[order(x$avg_log2FC, decreasing = T),]
    genelist <- structure(x$avg_log2FC, names = rownames(x))
    res <- GSEA(genelist, TERM2GENE = genesets, eps = 0)
    return(res)
  })
  # 导出表格
  for(i in seq_along(res.list)){
    res <- data.frame(res.list[[i]])
    write.csv(res, paste0("GSEA_KEGG", "/",names(res.list)[i], ".csv"), row.names = F)
  }
  # 导出图形
  
  for(i in seq_along(res.list)){
    res <- res.list[[i]]
    for(j in seq_along(res@result$ID)){
      #p <- gseaplot(res, geneSetID = j, title = res@result$ID[j], by = "runningScore")
      p <- enrichplot::gseaplot2(res, geneSetID = j, title = res@result$ID[j])
      filename <- paste0("GSEA_KEGG/", names(res.list)[i], "_", res@result$ID[j], '.pdf')
      
      # 绘制图形并保存
      pdf(file = filename, width = 8, height = 6)
      print(p)
      dev.off()
    }
  }
  
  #这里本来用ggsave保存的，但是会出现报错，所以后面就用pdf来保存。Error in UseMethod("grid.draw") :
  #no applicable method for 'grid.draw' applied to an object of class "c('gglist', 'list')"
  
}


### hallmark数据集
{
  # 设置基因集
  genesets = msigdbr(species = "Mus musculus", category = "H")
  genesets <- subset(genesets, select = c("gs_name", "gene_symbol"))
  # GSEA分析
  res.list <- lapply(deg.list, FUN = function(x){
    x = x[order(x$avg_log2FC, decreasing = T),]
    genelist <- structure(x$avg_log2FC, names = rownames(x))
    res <- GSEA(genelist, TERM2GENE = genesets, eps = 0)
    return(res)
  })
  # 导出表格
  for(i in seq_along(res.list)){
    res <- data.frame(res.list[[i]])
    write.csv(res, paste0("GSEA_HALLMARK/", names(res.list)[i], ".csv"), row.names = F)
  }
  # 导出图形
  for(i in seq_along(res.list)){
    res <- res.list[[i]]
    for(j in seq_along(res@result$ID)){
      #p <- gseaplot(res, geneSetID = j, title = res@result$ID[j], by = "runningScore")
      p <- enrichplot::gseaplot2(res, geneSetID = j, title = res@result$ID[j])
      filename <- paste0("GSEA_HALLMARK/", names(res.list)[i], "_", res@result$ID[j], '.pdf')
      # 绘制图形并保存
      pdf(file = filename, width = 8, height = 6)
      print(p)
      dev.off()
    }
  }
  
}

###### 分组平均水平GSVA
{
  setwd("~/project/Function")
  ### 平均表达矩阵
  stRNA <- readRDS("~/project/Data/sto_cluster.rds")
  DefaultAssay(stRNA) <- "SCT"
  colnames(stRNA@meta.data)
  expr <- AverageExpression(stRNA, assays = "SCT", slot = "data", group.by = "seurat_clusters")[[1]]
  expr <- expr[rowSums(expr)>0,]  #选取非零基因
  expr <- as.matrix(expr)
  
  ### Hallmark基因集的变异分析
  # 选择基因集
  genesets <- msigdbr(species = "Mus musculus", category = "H") 
  genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
  genesets <- split(genesets$gene_symbol, genesets$gs_name)
  # gsva默认开启全部线程计算
  #芯片数据属于连续型数据，故用正态分布（即kcdf = Gaussian）
  #RNA高通量测序数据属于离散型数据，故用泊松分布（即kcdf = Poisson）
  #如果RNA高通量测序数据经过log-CPMs, log-RPKMs 或者 log-TPMs处理，那kcdf需设置为Gaussian
  gsva.res <- gsva(expr, genesets, method="ssgsea", kcdf="Gaussian", parallel.sz=1) 
  saveRDS(gsva.res, "GSVA/ssGSEA_hallmark.rds")
  gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
  write.csv(gsva.df, "GSVA/ssGSEA_hallmark.csv", row.names = F)
  # 热图展示
  pheatmap::pheatmap(gsva.res, scale = "row", filename = "GSVA/ssGSEA_hallmark.pdf", width = 12, height = 10)
  
  ### KEGG基因集的变异分析
  # 选取基因集
  genesets <- msigdbr(species = "Mus musculus", category = "C2") 
  genesets <- subset(genesets, gs_subcat=="CP:KEGG", select = c("gs_name", "gene_symbol")) %>% as.data.frame()
  genesets <- split(genesets$gene_symbol, genesets$gs_name)
  # gsva默认开启全部线程计算
  gsva.res <- gsva(expr, genesets, method="ssgsea", kcdf="Gaussian", parallel.sz=1) 
  saveRDS(gsva.res, "GSVA/ssGSEA_kegg.rds")
  gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
  write.csv(gsva.df, "GSVA/ssGSEA_kegg.csv", row.names = F)
  # 热图展示
  pheatmap::pheatmap(gsva.res, scale = "row", filename = "GSVA/ssGSEA_kegg.pdf", width = 12, height = 30)
  
  ###### 单细胞水平基因集AUCell评分
  setwd("/home/data/t070429/spatial/processed data/Integrate/Function/AUCell")
  
  ### 提取表达矩阵
  #stRNA <- readRDS("~/project/Data/sto_cluster.rds")
  expr <- as.matrix(stRNA@assays$SCT@data)
  expr <- expr[rowSums(expr)>0,]  #选取非零基因
  
  ### 选择Hallmark基因集
  genesets <- msigdbr(species = "Mus musculus", category = "H") 
  genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
  genesets <- split(genesets$gene_symbol, genesets$gs_name)
  
  ### ssGSEA评分
  #gsva.res <- gsva(expr, genesets, method="ssgsea", kcdf="Gaussian", parallel.sz=1) 
}


### AUCell评分
{
  # 排序
  T0 = Sys.time()
  cells_rankings <- AUCell_buildRankings(expr, nCores = 2)
  Sys.time() - T0  #Time difference of 53.65928 secs
  # 评分
  install.packages("doMC")
  library("doMC")
  T0 = Sys.time()
  cells_AUC <- AUCell_calcAUC(genesets, cells_rankings, nCores = 4)
  Sys.time() - T0  #Time difference of 24.97522 secs
  
  ### 结果可视化
  aucMat <- getAUC(cells_AUC)
  write.csv(aucMat, file = "AUCell_hallmark.csv")
  annodf <- data.frame(cluster = paste0("cluster", stRNA$spatial.cluster), row.names = colnames(stRNA))
  pheatmap::pheatmap(aucMat, annotation_col = annodf, show_colnames = F, cluster_cols = T,
                     color = colorRampPalette(c("#344CB7", "white", "#CD1818"))(50),
                     scale = "row", breaks = unique(seq(-3, 3, length=50)), fontsize_col = 5,
                     filename = "AUCell_hallmark.pdf", width = 12, height = 8)
  
  
}





##===========================Step06：区域分析=======
##new cluster ids
{
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
}

## Assign new cluster IDs
{
  new.cluster.ids <- c("Follicle", "Oviduct smooth muscle", "Antral", "CL", "Atretic", "Stroma",
                       "Antral", "Oocyte", "CL",
                       "Epithelium")
  
  names(new.cluster.ids) <- levels(stRNA)
  
  stRNA <- RenameIdents(stRNA, new.cluster.ids)
  stRNA[["new.cluster.ids"]] <- Idents(object = stRNA)
  
  
  table(stRNA$new.cluster.ids)
  #Follicle Oviduct smooth muscle                Antral 
  #164                   107                   145 
  #CL               Atretic                Stroma 
  #124                    79                    63 
  #Oocyte            Epithelium 
  #50                    29 
  
  
  ## Save new ID analysis
  saveRDS(stRNA, file = "stRNA.newIDs.rds")
  setwd("~/spatial/processed data/Integrate/Area")
  
}

##TSNE with labels
{
  DimPlot(stRNA, reduction = "tsne", label = TRUE, label.size = 4.75) +
    plot_annotation(title = '') +
    labs(x = "TSNE1", y = "TSNE2", title = "") +
    NoLegend() +
    expand_limits(x = 11) +
    theme(text = element_text(size=20),
          axis.text = element_text(size = 20))
  
}


## 细胞比例图summarizing mouse proportions per cluster
{
  clusters <- as.data.frame(stRNA@meta.data) %>%
    group_by(orig.ident, new.cluster.ids, .drop = FALSE) %>%
    summarize(count = n()) %>%
    group_by(orig.ident) %>%
    mutate(Proportion = count/sum(count))
  
  
  ggplot(clusters, aes(x = orig.ident, y = Proportion)) +
    geom_col(aes(fill = new.cluster.ids)) +
    labs(fill = "seurat_clusters", x = "Sample") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_classic() +
    theme(text = element_text(size=20))+
    scale_fill_manual(values = mycolor)
  
}

## spot比例图
## group proportions per cluster
{
  aclusters <- filter(clusters, orig.ident == "CTR")
  yclusters <- filter(clusters, orig.ident == "Lip")
  aclusters$remaining <- sum(aclusters$count)-aclusters$count
  yclusters$remaining <- sum(yclusters$count)-yclusters$count
  clusters <- left_join(yclusters, aclusters, by = "new.cluster.ids")
  
  clusters$total <- clusters$count.x + clusters$count.y
  clusters$total_prop <- clusters$total/(sum(clusters$count.x + clusters$count.y))
  
  
  ggplot(clusters, aes(x = new.cluster.ids, y = total_prop, fill = new.cluster.ids)) +
    geom_col() +
    labs(x = "", y = "Proportion of Total Spots") +
    theme_classic() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    theme(text = element_text(size=20)) +
    theme(plot.margin = margin(10, 10, 10, 50))+scale_fill_manual(values = mycolor)
  
}

###more abundant in CTR and LIP,有统计学的
{
  ## fisher.exact on all rows
  fisher <- apply(clusters, 1, 
                  function(x) {
                    tbl <- matrix(as.numeric(x[c(3,5,7,9)]), ncol=2, byrow=T)
                    fisher.test(tbl, alternative="two.sided")
                  })
  clusters$p <- unlist(lapply(fisher, function(i) i$p.value))
  clusters$OR <- unlist(lapply(fisher, function(i) i$estimate))
  clusters$FDR <- p.adjust(clusters$p, method = "BH")
  clusters$logOR <- log10(clusters$OR) 
  clusters$logFDR <- -log10(clusters$FDR)
  
  ## replace infinite ORs with +-2 for plotting purposes
  
  clusters$logOR <- ifelse(clusters$logOR == "Inf", 2,
                           ifelse(clusters$logOR == "-Inf", -2, clusters$logOR))
  
  
  ggplot(clusters, aes(x = logOR, y = logFDR, color = new.cluster.ids, label = new.cluster.ids)) +
    geom_point(size = 5) +
    geom_text_repel(min.segment.length = 0.5,
                    box.padding = 0.5, point.padding = 0.5, max.overlaps = 20,
                    size = 4, fontface = "bold", segment.size = 1, nudge_y = .4) +
    geom_hline(yintercept = 1.3) +
    geom_vline(xintercept = c((-log10(0.67)), (-log10(1.5))), linetype = 2) +
    labs(y = "-log10 FDR", x = "log10 Odds Ratio") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none") +
    theme(text = element_text(size=20))+
    ggtitle("more abundant in CTR")+
    scale_fill_manual(values = mycolor)
  
}

#spatialDimplot

SpatialDimPlot(stRNA, images = "CTR1", crop = FALSE, pt.size.factor = 1.5) +
  labs(fill = "")

SpatialDimPlot(stRNA, images = "LIP2", crop = FALSE, pt.size.factor = 1.5) +
  labs(fill = "")

#area DEG

#如果省略了ident.2参数或将其设置为NULL，FindMarkers函数将对指定的ident.1组(对照组)与其他所有组之间进行差异表达分析。



# 比较Antral和follicle的差异基因
gene_list_antral <- FindMarkers(stRNA, assay = "Spatial", ident.1 = "Follicle", ident.2 = "Antral", 
                                logfc.threshold = 0, latent.vars = "SectionID", test.use = "LR")
write.csv(gene_list_antral,"DEG Antral VS follicle.csv")
# 比较Antral和CL的差异基因 (对照是窦卵泡)
gene_list_CL <- FindMarkers(stRNA, assay = "Spatial", ident.1 = "Antral", ident.2 = "CL", 
                            logfc.threshold = 0, latent.vars = "SectionID", test.use = "LR")


write.csv(gene_list_CL,"CL VS Antral.csv")

#韦恩图（VennDiagram 包，适用样本数 2-5）

{library(VennDiagram)
  
  #读入作图文件，all.txt即上述提到的记录group1-4的元素名称的文件
  dat <- read.table('all.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
  
  #以2个分组为例
  #指定统计的分组列，并设置作图颜色、字体样式等
  venn_list <- list(group1 = gene.list$cluster0, group2 = gene.list$cluster2,group3=gene.list$cluster6)
  
  venn.diagram(venn_list, 
               filename = 'venn_Follicle_and_antral.png', imagetype = 'png', 
               resolution = 300,
               fill =mycolor[1:3], alpha = 0.50, cat.col = rep('black', 3), 
               col = 'black', cex = 1.5, fontfamily = 'serif', 
               main = "Venn C0-C2-C6 ",main.cex = 2,
               cat.cex = 1.5, cat.fontfamily = 'serif',
               category.names = c("cluster 0", "cluster 2","cluster 6"
               ),
               output=TRUE)
  
  #继续以上述3个分组为例，组间交集元素获得
  inter <- get.venn.partitions(venn_list)
  for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
  write.table(inter[-c(5, 6)], 'venn_Follicle_and_antral.txt', row.names = FALSE, sep = '\t', quote = FALSE)
  
  #TRUE代表该组中出现的元素，FALSE则代表未出现的元素
  #count为交集元素数量，values为交集元素名称
  #画窦卵泡和黄体的差异基因图
  venn_list <- list(group1 = gene.list$cluster2, group2 = gene.list$cluster6,group3=gene.list$cluster3,group4=gene.list$cluster8)
  
  venn.diagram(venn_list, 
               filename = 'venn_Antral_and_CL.png', imagetype = 'png', 
               resolution = 300,
               fill =mycolor[1:4], alpha = 0.50, cat.col = rep('black', 4), 
               col = 'black', cex = 1.5, fontfamily = 'serif', 
               main = "Venn C2_6-C3-C8 ",main.cex = 2,
               cat.cex = 1.5, cat.fontfamily = 'serif',
               category.names = c("cluster 2","cluster 6", "cluster 3","cluster 8"),
               output=TRUE)
  
  
  inter <- get.venn.partitions(venn_list)
  for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
  write.table(inter[-c(5, 6)], 'venn_CL and_antral.txt', row.names = FALSE, sep = '\t', quote = FALSE)
  
  rm(aclusters,clusters,inter,yclusters,venn_list,i,SectionID,venn.plot)
  rm(new.cluster.ids)
  
  ##画5群的
  
  venn_list <- list(group1 = gene.list$cluster2, group2 = gene.list$cluster6,group3=gene.list$cluster3,group4=gene.list$cluster8,group5 = gene.list$cluster0 )
  
  venn.diagram(venn_list, 
               filename = 'venn_Follicle Antral_and_CL.png', imagetype = 'png', 
               resolution = 300,
               fill =mycolor[1:5], alpha = 0.50, cat.col = rep('black', 5), 
               col = 'black', cex = 1.5, fontfamily = 'serif', 
               main = "VennC0-C2_6-C3-C8 ",main.cex = 2,
               cat.cex = 1.5, cat.fontfamily = 'serif',
               category.names = c("cluster 2","cluster 6", "cluster 3","cluster 8","cluster 0"),
               output=TRUE)
  
  # 使用cells.highlight参数高亮感兴趣的一些细胞 
  Idents(stRNA)
  SpatialDimPlot(stRNA, cells.highlight = CellsByIdentities(object = stRNA, idents = c("Follicle", "Antral","CL")), facet.highlight = TRUE, ncol = 3) 
  
}

#堆叠小提琴图
{
  VlnPlot(stRNA, 
          features = c("Cd68","Enpp2","Cd52"),
          pt.size = 0,
          ncol = 3)
  
  movlnplot<- function(obj,
                       feature,
                       pt.size = 0,
                       plot.margin = unit(c(-1, 1, -0.5, 1), "cm"),
                       ...
  ) {
    p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  +
      xlab("") + ylab(feature) + ggtitle("")+ 
      scale_fill_manual(values = paletteer::paletteer_d('ggsci::category20c_d3'))+
      theme(legend.position = "none",
            panel.spacing = unit(x = 0, units = 'lines'),
            axis.line = element_blank(),
            panel.background = element_rect(fill = NA, color = 'black'),
            strip.background = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(),
            axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
            axis.ticks.length = unit(x = 0, units = 'cm'))
    return(p)
  }
  
  FVlnPlot<- function(obj, features,
                      pt.size = 0,
                      plot.margin = unit(c(0, 1, 0,1), "cm"),
                      ...) {
    
    plot_list<- purrr::map(features, function(x) movlnplot(obj = obj,feature = x, ...))
    plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
      theme(panel.spacing = unit(x = 0, units = 'lines'),
            axis.text.x=element_text(), 
            axis.ticks.x = element_line(),
            plot.margin =plot.margin,
            strip.background = element_blank())
    
    p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
    return(p)
  }
  
  FVlnPlot(stRNA, 
           features = c("Cd68","Enpp2","Cd52"), pt.size=0)
  
  # 从stRNA对象中选择包含follicle、antral和CL的
  stRNA_subset <- subset(stRNA, idents = c("Follicle", "Antral", "CL"))
  stRNA_subset1 <- subset(stRNA, idents = c("Follicle", "Antral"))
  
  
  # 绘制小提琴图
  FVlnPlot(stRNA_subset1, 
           features = c("Cd68","Enpp2","Cd52"),
           pt.size = 0
  )
  
  
  #两组之间比较
  VlnPlot(stRNA_subset, features = c("Cd68","Enpp2","Cd52"),
          stack=T,pt.size=0,
          split.by = 'orig.ident',
          flip = T,
          add.noise = T)+#横纵轴不标记任何东西
    theme(axis.text.y = element_blank(), #不显示坐标刻度
          axis.ticks.y = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(colour = 'black',size = 10,angle = 90),
          legend.position = 'top',
          legend.title=element_blank(),
          legend.box.background = element_blank(),
          legend.text = element_text(color="black",size=10),
          legend.spacing.x=unit(0.2,'cm'),
          legend.key.width=unit(0.4,'cm'),
          legend.key.height=unit(0.4,'cm'),
          legend.background=element_blank())
}

#区域热图
{
  allMarkersSpatial <- FindAllMarkers(stRNA_subset, assay = "Spatial", only.pos = TRUE)
  topgene <- allMarkersSpatial %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
  
  write.csv(allMarkersSpatial,"allMarkersFollicle_Antral_CL.csv")
  DoHeatmap(stRNA_subset,
            features = as.character(unique(topgene$gene)),
            group.colors =mycolor)+
    scale_fill_gradientn(colors = c("white","grey90","firebrick3"))
  
}



library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)



##===========================Step07：区域cellchat====
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Seurat)

##每组之间单独比较
{
  # Prepare input data for CelChat analysis
  stRNA_CTR <- sto.list[[1]]  
  stRNA_LIP <- sto.list[[2]]
  new.cluster.ids <- c("Follicle", "Oviduct smooth muscle", "Antral", "CL", "Atretic", "Stroma",
                       "Antral", "Oocyte", "CL",
                       "Epithelium")
  
  names(new.cluster.ids) <- levels(stRNA_LIP)
  stRNA_LIP <- RenameIdents(stRNA_LIP, new.cluster.ids)
  stRNA_LIP[["new.cluster.ids"]] <- Idents(object = stRNA_LIP)
  
  
  data.input = GetAssayData(stRNA_LIP, slot = "data", assay = "SCT") # normalized data matrix
  meta = data.frame(labels = Idents(stRNA_LIP), row.names = names(Idents(stRNA_LIP))) # manually create a dataframe consisting of the cell labels
  unique(meta$labels) # check the cell labels
  
  # load spatial imaging information
  # Spatial locations of spots from full (NOT high/low) resolution images are required
  
  coords2 <- stRNA@images[["LIP2"]]@coordinates
  coords2<- coords2[,c(4:5)]
  spatial.locs <- coords2
  spatial.locs = GetTissueCoordinates(stRNA_LIP, scale = NULL, cols = c("imagerow", "imagecol")) 
  # Scale factors and spot diameters of the full resolution images 
  scale.factors = jsonlite::fromJSON(txt = file.path("/home/data/t070429/spatial/raw data/project/LIP/spatial/", 'scalefactors_json.json'))
  scale.factors = list(spot.diameter = 55, spot = scale.factors$spot_diameter_fullres, # these two information are required
                       fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef, lowres = scale.factors$tissue_lowres_scalef # these three information are not required
  )
  #Create a CellChat object
  
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                             datatype = "spatial", coordinates = spatial.locs, scale.factors = scale.factors)
  
  
  CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data
  
  # use a subset of CellChatDB for cell-cell communication analysis
  CellChatDB.use <- CellChatDB # simply use the default CellChatDB
  #CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
  
  
  cellchat@DB <- CellChatDB.use
  # subset the expression data of signaling genes for saving computation cost
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  future::plan("multiprocess", workers = 4) # do parallel
  
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  # 检查cellchat数据框中是否包含缺失值、无穷大或未定义的值
  summary(cellchat)
  
  
  
  cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
                                distance.use = TRUE, interaction.length = 200, scale.distance = 0.01)
  
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  cellchat <- computeCommunProbPathway(cellchat)
  
  cellchat <- aggregateNet(cellchat)
  groupSize <- as.numeric(table(cellchat@idents))
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  
  groupLIP.net <- subsetCommunication(cellchat) 
  write.csv(groupLIP.net, file = "groupLIP_net_inter.csv")
  
  save(cellchat, stRNA_LIP, groupLIP.net, file = "groupLIP_cellchat.RData")
  
  cellchat@netP$pathways
  
  pathways.show <- c("CXCL") 
  # Circle plot
  par(mfrow=c(1,1))
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
  
  # Spatial plot
  par(mfrow=c(1,1))
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5)
  
  # Compute the network centrality scores
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
  # Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
  par(mfrow=c(1,1))
  netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
  
  
  # USER can visualize this information on the spatial imaging, e.g., bigger circle indicates larger incoming signaling
  par(mfrow=c(1,1))
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", edge.width.max = 2, alpha.image = 0.2, vertex.weight = "incoming", vertex.size.max = 3, vertex.label.cex = 3.5)
  
  netAnalysis_contribution(cellchat, signaling = pathways.show)
  
  plotGeneExpression(cellchat, signaling = "CXCL")
  
  
  gg1 <- netAnalysis_signalingRole_scatter(cellchat)
  gg1 <- netAnalysis_signalingRole_scatter(cellchat)
  plot=gg1 + gg2
  ggsave(filename = "netAnalysis_signalingRole_scatter_LIP.pdf",width = 9,height = 4,plot =gg1 + gg2 )
  
  
  library(NMF)
  library(ggalluvial)
  selectK(cellchat, pattern = "outgoing")
  
  Cairo::CairoPNG(filename = "identifyCommunicationPatterns_outgoing_LIP.png",width = 1800,height = 900,dpi = 300)
  nPatterns = 3
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
  dev.off()
  
  Cairo::CairoPNG(filename = "netAnalysis_river_outgoing_LIP.png",width = 1800,height = 900,dpi = 300)
  netAnalysis_river(cellchat, pattern = "outgoing")
  dev.off()
  
  Cairo::CairoPNG(filename = "netAnalysis_dot_outgoing_LIP.png",width = 1800,height = 900,dpi = 300)
  netAnalysis_dot(cellchat, pattern = "outgoing")
  dev.off()
  
  
  selectK(cellchat, pattern = "incoming")
  nPatterns = 4
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
  Cairo::CairoPNG(filename = "netAnalysis_river_incoming_LIP.png",width = 1800,height = 900,dpi = 300)
  netAnalysis_river(cellchat, pattern = "incoming")
  dev.off()
  
  Cairo::CairoPNG(filename = "netAnalysis_dot_incoming_LIP.png",width = 1800,height = 900,dpi = 300)
  netAnalysis_dot(cellchat, pattern = "incoming")
  dev.off()
  
  
  ##Identify signaling groups based on their functional similarity
  
  cellchat <- computeNetSimilarity(cellchat, type = "functional")
  cellchat <- netEmbedding(cellchat, type = "functional",umap.method = "uwot")
  cellchat <- netClustering(cellchat, type = "functional",k = 5,do.parallel = FALSE)
  netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
  
  
  ####Identify signaling groups based on structure similarity
  
  
  cellchat <- computeNetSimilarity(cellchat, type = "structural")
  cellchat <- netEmbedding(cellchat, type = "structural",umap.method = "uwot")
  #> Manifold learning of the signaling networks for a single dataset
  cellchat <- netClustering(cellchat, type = "structural",do.parallel = FALSE)
  #> Classification learning of the signaling networks for a single dataset
  # Visualization in 2D-space
  netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
  LIP_cellchat <- cellchat
  
}

##两组之间比较
{
  object.list <- list(CTR = CTR_cellchat, LIP = LIP_cellchat)
  cellchat <- mergeCellChat(object.list, add.names = names(object.list))
  #Merge the following slots: 'data.signaling','images','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
  
  #比较交互总数和交互强度
  gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
  gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
  p<-gg1 + gg2
  
  p
  
  gg1 <- netVisual_heatmap(cellchat)
  #> Do heatmap based on a merged object
  gg2 <- netVisual_heatmap(cellchat, measure = "weight")
  #> Do heatmap based on a merged object
  gg1 + gg2
  
  #可视化两组之间信号接受与发出的主要细胞群
  for (i in 1:length(object.list)) {
    object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]])
  }
  num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
  weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
  gg <- list()
  # Plot the signaling role scores for each cell type in each network
  gg <- list()
  for (i in 1:length(object.list)) {
    gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
  }
  patchwork::wrap_plots(plots = gg)
  
  #比较每个信号通路的整体信息流
  gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
  gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
  ggsave(filename = "rankNetstatCTR_LIP.pdf",width = 6,height = 8,plot = gg1 + gg2)
  
  
  gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = FALSE)
  gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = FALSE)
  ggsave(filename = "rankNetnonstat.pdf",width = 6,height = 9,plot = gg1 + gg2)
}






##===========================Step08：SPATA空间向量分析=======
library(SPATA2)
library(metR)
library(magrittr)
library(ggplot2)
library(dplyr)
library(patchwork)
library(monocle3)

library(Seurat)
###空间向量轨迹分析
stRNA_CTR@images[[2]]<-NULL

{
  
  
  runVectorFields <- function(object,
                              features,
                              cut_off=NULL,
                              normalize=T,
                              smooth=T,
                              smooth_span=NULL,
                              dist.spot=10,
                              run.mcor=T,
                              workers=8,
                              ram=50){
    #get the data:
    df <- SPATA2::hlpr_join_with_aes(object,
                                     df=SPATA2::getCoordsDf(object),
                                     color_by = features,
                                     normalize=normalize,
                                     smooth=smooth,
                                     smooth_span=smooth_span)
    if(!is.null(cut_off)){df[df[,features]<cut_off,features]=0}
    
    #prepare data
    NN.file <- SPATAwrappers::getSurroundedSpots(object)
    
    if(run.mcor==T){
      base::options(future.fork.enable=TRUE)
      future::plan("multisession",workers=workers)
      future::supportsMulticore()
      base::options(future.globals.maxSize=ram*100*1024^2)
      message("...Run multicore...")
    }
    
    if(dist.spot<min(NN.file%>%filter(distance!=0)%>%pull(distance))){
      dist.spot<min(NN.file%>%filter(distance!=0)%>%pull(distance))
      message(paste0("The distance was adopted for the minimal distance:", dist.spot,"px"))
    }
    
    VF <- furrr::future_map(.x=1:nrow(df),.f=function(i){
      #spot Def.
      bc <- df[i,c("barcodes")]
      cc <- df[i,c("x","y")]
      
      #neighbour spots
      NN <- 
        NN.file%>%
        dplyr::filter(xo<cc$x+dist.spot&xo>cc$x-dist.spot)%>%
        dplyr::filter(yo<cc$y+dist.spot&yo>cc$y-dist.spot)%>%
        dplyr::pull(bc_destination)
      
      
      #filter input DF
      NN.df <- df%>%dplyr::filter(barcodes %in% NN)%>%as.data.frame()
      
      parameter <- features
      #Create Vector
      
      V <- -c(as.numeric(cc)-c(NN.df$x[which.max(NN.df[,parameter])],NN.df$y[which.max(NN.df[,parameter])]))
      
      if(length(V)==0){out <- data.frame(barcodes=bc,t.x=0,t.y=0)}else {out <- data.frame(barcodes=bc,t.x=V[1],t.y=V[2])}
      return(out)
      
    }, .progress=T)%>%
      do.call(rbind,.)%>%
      as.data.frame()
    
    out <- cbind(df,VF)
    out [is.na(out)] <- 0
    
    return(out)
  }
  
  plotVectorFields <- function(VF,parameter,pt.size=6,pt.alpha=0.8,
                               color.extern=NULL,skip=1){
    VF <- 
      VF %>%
      dplyr::select(x,y,{{parameter}},t.x,t.y)%>%
      dplyr::rename("parameter":=!!sym(parameter))
    
    color.points <- VF$parameter
    if(!is.null(color.extern)){color.points <- color.extern}
    
    if(color.points%>% class()=="factor"){
      p <- 
        ggplot2::ggplot(data = VF,aes(x,y))+
        ggplot2::geom_point(data = VF,mapping = aes(x,y,color=color.points),size=pt.size,alpha=pt.alpha)+
        metR::geom_vector(aes(dx=t.x,dy=t.y),skip=skip)+
        metR:scale_mag()+
        ggplot2::theme_void()+
        Seurat::NoLegend()
    }else{
      p <- 
        ggplot2::ggplot(data = VF,aes(x,y))+
        ggplot2::geom_point(data = VF,mapping = aes(x,y,color=color.points),size=pt.size,alpha=pt.alpha)+
        ggplot2::scale_color_viridis_c(guide = "none")+
        metR::geom_vector(aes(dx=t.x,dy=t.y),skip=skip)+
        ggplot2::theme_void()+
        Seurat::NoLegend()+
        metR:scale_mag()
      
    }
    
    return(p)
  }
  x <- stRNA_CTR
  
  y <- transformSeuratToSpata(seurat_object = stRNA_CTR,
                              sample_name = "CTR1",
                              #method = "spatial",
                              #coords_from = "pca",
                              assay_name = "SCT")
  #可以是某个基因，也可以是某个细胞或者通路（目前后面两种方式还没研究到）
  c('Fst','Cd68','Enpp2','Gdf9','Esr2','Gja1','Lyz2','Sfrp4','Ptgfr','Star','Cd52','Spp1','Ctss','Wnt4','Ccl6','Ccl2','Gpnmb','Cd36','Foxo1','Runx2','Cd74')
  
  data=runVectorFields(y,'Cd68')
  head(data)
  
  
  
  #方法一：一个一个手动保存
  p=plotVectorFields(data,"Cd68")
  p
  
  
  
  #方法二：写一个函数，然后自动输出保存
  
  pdf('Cd68-RNA trajectory.pdf')
  print(p)
  dev.off()
  
}





save(stRNA_subset, stRNA_subset1, sto.list, SC1, mycolor, mypalette, stRNA, FVlnPlot, file = "OvarySpatial.Rdata")

