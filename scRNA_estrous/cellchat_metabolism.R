################################################################################
# step1 scoring
################################################################################
library(Seurat)
rm(list=ls())
options(stringsAsFactors = F)
load(file = "./rdata/scRNA_all_with_mac.rda")


################################################################################
# step2 cellchat
################################################################################

library(Seurat)
library(CellChat)
rm(list=ls())
options(stringsAsFactors = F)
# load(file = "./rdata/ov_scRNA_anno.rda")

# 我们先选取ovarian cancer观察一下细胞通讯
scRNA_chat <- scRNA
{
  meta =scRNA_chat@meta.data # a dataframe with rownames containing cell mata data
  data_input <- as.matrix(scRNA_chat@assays$RNA@data)
  identical(colnames(data_input),rownames(meta))
  # "Create a CellChat object from a data matrix"
  cellchat <- createCellChat(object = data_input, # support normalized expression matrix
                             meta = meta, # meta.data
                             group.by = "Level0") # you need to change the celltype name for your object
  
  
  CellChatDB <- CellChatDB.mouse # (CellChatDB.human)&(CellChatDB.mouse)
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # you can choose Secreted Signaling、ECM-Receptor or Cell-Cell Contact
  cellchat@DB <- CellChatDB.use 
  cellchat <- subsetData(cellchat) #choose the subset from CellChatDB.use->cellchat@data.Signaling
  ## calculate overexpressed genes
  cellchat <- identifyOverExpressedGenes(cellchat,thresh.pc = 0, #细胞比例阈值
                                         thresh.fc = 0, #差异倍数
                                         thresh.p = 0.05) #P-Value
  # head(cellchat@var.features$features.info) #差异计算结果表
  
  # Identify ligand-receptor interactions of overexpressed genes
  cellchat <- identifyOverExpressedInteractions(cellchat)
  # head(cellchat@LR$LRsig) #计算结果赋值位置
  
  # Project ligands and receptors into the PPI network
  cellchat <- projectData(cellchat, PPI.human)
  
  # Cell communication prediction
  ## Calculation of cell communication probability
  cellchat <- computeCommunProb(cellchat)  
  
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10) 
  
  # 提取信号通路水平的细胞通讯表
  cellchat <- computeCommunProbPathway(cellchat)
  
  # 提取配受体对细胞通讯结果表
  df.net<- subsetCommunication(cellchat)
}

write.csv(df.net,file ='./files/cellchat_OV.csv',quote=F)

cellchat_OV = cellchat
save(cellchat_OV,file = "./rdata/cellchat_OV.rda")



#计算整合的细胞类型之间通信结果
cellchat <- aggregateNet(cellchat_OV)


# 画图
source("utools.R")
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,color.use = ov_palette,
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,color.use = ov_palette,
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()



#指定受体-配体细胞类型
netVisual_bubble(cellchat, sources.use = c("Mac","Other_Immune"), 
                 targets.use = c("Granulosa","Mesenchyme"), remove.isolate = FALSE)

netVisual_bubble(cellchat, sources.use = c("Granulosa","Mesenchyme"), 
                 targets.use = c("Mac","Other_Immune"), remove.isolate = FALSE)

########################################################
rm(list = ls())
gc()
load(file = "./rdata/scRNA_all_with_mac.rda")
library(Seurat)
library(CellChat)

# load(file = "./rdata/ov_scRNA_anno.rda")

# 我们先选取ovarian cancer观察一下细胞通讯
# scRNA_chat <- scRNA

scRNA <-subset(scRNA,Level0 != "Oocyte")
# scRNA_chat@meta.data$Level0 = droplevels(scRNA_chat@meta.data$Level0 , exclude = setdiff(levels(scRNA_chat@meta.data$Level0 ),unique(scRNA_chat@meta.data$Level0)))
table(scRNA@meta.data$Level0 )
scRNA_chat <- subset(scRNA, condition =='P')
{
  meta =scRNA_chat@meta.data # a dataframe with rownames containing cell mata data
  data_input <- as.matrix(scRNA_chat@assays$RNA@data)
  identical(colnames(data_input),rownames(meta))
  # "Create a CellChat object from a data matrix"
  cellchat <- createCellChat(object = data_input, # support normalized expression matrix
                             meta = meta, # meta.data
                             group.by = "Level0")
  
  
  CellChatDB <- CellChatDB.mouse # (CellChatDB.human)&(CellChatDB.mouse)
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # you can choose Secreted Signaling、ECM-Receptor or Cell-Cell Contact
  cellchat@DB <- CellChatDB.use 
  cellchat <- subsetData(cellchat) #choose the subset from CellChatDB.use->cellchat@data.Signaling
  ## calculate overexpressed genes
  cellchat <- identifyOverExpressedGenes(cellchat,thresh.pc = 0, #细胞比例阈值
                                         thresh.fc = 0, #差异倍数
                                         thresh.p = 0.05) #P-Value
  # head(cellchat@var.features$features.info) #差异计算结果表
  
  # Identify ligand-receptor interactions of overexpressed genes
  cellchat <- identifyOverExpressedInteractions(cellchat)
  # head(cellchat@LR$LRsig) #计算结果赋值位置
  
  # Project ligands and receptors into the PPI network
  cellchat <- projectData(cellchat, PPI.human)
  
  # Cell communication prediction
  ## Calculation of cell communication probability
  cellchat <- computeCommunProb(cellchat)  
  
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10) 
  
  # 提取信号通路水平的细胞通讯表
  cellchat <- computeCommunProbPathway(cellchat)
  
  # 提取配受体对细胞通讯结果表
  df.net<- subsetCommunication(cellchat)
}

write.csv(df.net,file ='./files/cellchat_P.csv',quote=F)

cellchat_P = cellchat
save(cellchat_P,file = "./rdata/cellchat_P.rda")

# PEMD -E
# 我们先选取ovarian cancer观察一下细胞通讯
scRNA_chat <- subset(scRNA, condition =='E')
{
  meta =scRNA_chat@meta.data # a dataframe with rownames containing cell mata data
  data_input <- as.matrix(scRNA_chat@assays$RNA@data)
  identical(colnames(data_input),rownames(meta))
  # "Create a CellChat object from a data matrix"
  cellchat <- createCellChat(object = data_input, # support normalized expression matrix
                             meta = meta, # meta.data
                             group.by = "Level0")
  
  
  CellChatDB <- CellChatDB.mouse # (CellChatDB.human)&(CellChatDB.mouse)
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # you can choose Secreted Signaling、ECM-Receptor or Cell-Cell Contact
  cellchat@DB <- CellChatDB.use 
  cellchat <- subsetData(cellchat) #choose the subset from CellChatDB.use->cellchat@data.Signaling
  ## calculate overexpressed genes
  cellchat <- identifyOverExpressedGenes(cellchat,thresh.pc = 0, #细胞比例阈值
                                         thresh.fc = 0, #差异倍数
                                         thresh.p = 0.05) #P-Value
  # head(cellchat@var.features$features.info) #差异计算结果表
  
  # Identify ligand-receptor interactions of overexpressed genes
  cellchat <- identifyOverExpressedInteractions(cellchat)
  # head(cellchat@LR$LRsig) #计算结果赋值位置
  
  # Project ligands and receptors into the PPI network
  cellchat <- projectData(cellchat, PPI.human)
  
  # Cell communication prediction
  ## Calculation of cell communication probability
  cellchat <- computeCommunProb(cellchat)  
  
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10) 
  
  # 提取信号通路水平的细胞通讯表
  cellchat <- computeCommunProbPathway(cellchat)
  
  # 提取配受体对细胞通讯结果表
  df.net<- subsetCommunication(cellchat)
}

write.csv(df.net,file ='./files/cellchat_E.csv',quote=F)

cellchat_E = cellchat
save(cellchat_E,file = "./rdata/cellchat_E.rda")

# PEMD -M
# 我们先选取ovarian cancer观察一下细胞通讯
scRNA_chat <- subset(scRNA, condition =='M')
{
  meta =scRNA_chat@meta.data # a dataframe with rownames containing cell mata data
  data_input <- as.matrix(scRNA_chat@assays$RNA@data)
  identical(colnames(data_input),rownames(meta))
  # "Create a CellChat object from a data matrix"
  cellchat <- createCellChat(object = data_input, # support normalized expression matrix
                             meta = meta, # meta.data
                             group.by = "Level0")
  
  
  CellChatDB <- CellChatDB.mouse # (CellChatDB.human)&(CellChatDB.mouse)
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # you can choose Secreted Signaling、ECM-Receptor or Cell-Cell Contact
  cellchat@DB <- CellChatDB.use 
  cellchat <- subsetData(cellchat) #choose the subset from CellChatDB.use->cellchat@data.Signaling
  ## calculate overexpressed genes
  cellchat <- identifyOverExpressedGenes(cellchat,thresh.pc = 0, #细胞比例阈值
                                         thresh.fc = 0, #差异倍数
                                         thresh.p = 0.05) #P-Value
  # head(cellchat@var.features$features.info) #差异计算结果表
  
  # Identify ligand-receptor interactions of overexpressed genes
  cellchat <- identifyOverExpressedInteractions(cellchat)
  # head(cellchat@LR$LRsig) #计算结果赋值位置
  
  # Project ligands and receptors into the PPI network
  cellchat <- projectData(cellchat, PPI.human)
  
  # Cell communication prediction
  ## Calculation of cell communication probability
  cellchat <- computeCommunProb(cellchat)  
  
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10) 
  
  # 提取信号通路水平的细胞通讯表
  cellchat <- computeCommunProbPathway(cellchat)
  
  # 提取配受体对细胞通讯结果表
  df.net<- subsetCommunication(cellchat)
}

write.csv(df.net,file ='./files/cellchat_M.csv',quote=F)

cellchat_M = cellchat
save(cellchat_M,file = "./rdata/cellchat_M.rda")

# PEMD-D
scRNA_chat <- subset(scRNA, condition =='D')
{
  meta =scRNA_chat@meta.data # a dataframe with rownames containing cell mata data
  data_input <- as.matrix(scRNA_chat@assays$RNA@data)
  identical(colnames(data_input),rownames(meta))
  # "Create a CellChat object from a data matrix"
  cellchat <- createCellChat(object = data_input, # support normalized expression matrix
                             meta = meta, # meta.data
                             group.by = "Level0")
  
  
  CellChatDB <- CellChatDB.mouse # (CellChatDB.human)&(CellChatDB.mouse)
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # you can choose Secreted Signaling、ECM-Receptor or Cell-Cell Contact
  cellchat@DB <- CellChatDB.use 
  cellchat <- subsetData(cellchat) #choose the subset from CellChatDB.use->cellchat@data.Signaling
  ## calculate overexpressed genes
  cellchat <- identifyOverExpressedGenes(cellchat,thresh.pc = 0, #细胞比例阈值
                                         thresh.fc = 0, #差异倍数
                                         thresh.p = 0.05) #P-Value
  # head(cellchat@var.features$features.info) #差异计算结果表
  
  # Identify ligand-receptor interactions of overexpressed genes
  cellchat <- identifyOverExpressedInteractions(cellchat)
  # head(cellchat@LR$LRsig) #计算结果赋值位置
  
  # Project ligands and receptors into the PPI network
  cellchat <- projectData(cellchat, PPI.human)
  
  # Cell communication prediction
  ## Calculation of cell communication probability
  cellchat <- computeCommunProb(cellchat)  
  
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10) 
  
  # 提取信号通路水平的细胞通讯表
  cellchat <- computeCommunProbPathway(cellchat)
  
  # 提取配受体对细胞通讯结果表
  df.net<- subsetCommunication(cellchat)
}

write.csv(df.net,file ='./files/cellchat_D.csv',quote=F)

cellchat_D = cellchat
save(cellchat_D,file = "./rdata/cellchat_D.rda")


#############################################################
# PEMD
source("utools.R")
# calculate numbers
load("./rdata/cellchat_P.rda")
load("./rdata/cellchat_E.rda")
load("./rdata/cellchat_M.rda")
load("./rdata/cellchat_D.rda")

cellchat_P <- aggregateNet(cellchat_P)
cellchat_E <- aggregateNet(cellchat_E)
cellchat_M<- aggregateNet(cellchat_M)
cellchat_D <- aggregateNet(cellchat_D)
# calculate strength

# Compute the network centrality scores
cellchat_P <- netAnalysis_computeCentrality(cellchat_P, slot.name = "netP")
cellchat_E <- netAnalysis_computeCentrality(cellchat_E, slot.name = "netP")
cellchat_M <- netAnalysis_computeCentrality(cellchat_M, slot.name = "netP")
cellchat_D <- netAnalysis_computeCentrality(cellchat_D, slot.name = "netP")

dev.off()
source("utools.R")
par(mfrow = c(1,4), xpd=TRUE)

groupSize <- as.numeric(table(cellchat_P@idents))

netVisual_circle(cellchat_P@net$count, vertex.weight = groupSize,color.use = ov_palette, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")

groupSize <- as.numeric(table(cellchat_E@idents))

netVisual_circle(cellchat_E@net$count, vertex.weight = groupSize, color.use = ov_palette,
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")

groupSize <- as.numeric(table(cellchat_M@idents))

netVisual_circle(cellchat_M@net$count, vertex.weight = groupSize, color.use = ov_palette,
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")

groupSize <- as.numeric(table(cellchat_D@idents))

netVisual_circle(cellchat_D@net$count, vertex.weight = groupSize, color.use = ov_palette,
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")


# netVisual_circle(cellchat_P@net$weight, vertex.weight = groupSize, 
#                  weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")




object.list <- list(P = cellchat_P, E = cellchat_E,M=cellchat_M,D=cellchat_D)
cellchat <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)

# 如果数据有三组及以上， group = c(1,2)改为 group = c(1,2,3)，以此类推。
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3,4),color.use =ov_palette)
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3,4), measure = "weight",color.use =ov_palette)
gg1 + gg2

#差异性heatmap
gg1 <- netVisual_heatmap(cellchat,comparison = c(1,2,3,4))
gg2 <- netVisual_heatmap(cellchat, measure = "weight",comparison = c(1,2,3,4))
gg1 + gg2



library(ComplexHeatmap)
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all",# 还可以设置成outgoing和ingoing
                                        signaling = pathway.union, title = names(object.list)[i], 
                                        width = 5, height = 12, color.heatmap = "RdPu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all",
                                        signaling = pathway.union, title = names(object.list)[i+1], 
                                        width = 5, height = 12, color.heatmap = "RdPu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,comparison = c(1,2,3,4),color.use =ov_palette[1:4])
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,comparison = c(1,2,3,4),color.use =ov_palette[1:4])
gg1 + gg2


# output
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2),,color.use = c("#6894B9",  "#f1707d"))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,comparison = c(1,2),color.use = c("#6894B9",  "#f1707d"))
gg1 + gg2


#fig2.D
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}

## Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
## Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}

## Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
## Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways

#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)









################################################################################
# cellcall


load("./rdata/ov_scRNA_anno.rda")
library(Seurat)
library(cellcall)

scRNA_cellcall=scRNA

# 整理cellcall所需要的形式！
counts=as.data.frame(scRNA_cellcall@assays$RNA@counts)
colnames(counts)
colnames(counts)=stringr::str_replace(colnames(counts),pattern = '_',replacement = '-')
colnames(counts)
colnames(counts)=paste0(colnames(counts),'_',scRNA_cellcall$celltype)
colnames(counts)

gc()

# names.field标识-所分开的第几个元素
# names.delin最好为-

mt <- CreateNichConObject(data=counts, min.feature = 3,
                          names.field = 2,
                          names.delim = "_",   # 双下划线
                          source = "UMI",
                          scale.factor = 10^6,
                          Org = "Homo sapiens",
                          project = "Microenvironment")



mt <- TransCommuProfile(object = mt,
                        pValueCor = 0.05,
                        CorValue = 0.1,
                        topTargetCor=1,
                        p.adjust =0.05,  ## 若无通路，可放大到1
                        use.type="median",
                        probs = 0.9,
                        method="weighted",
                        IS_core = TRUE,
                        Org = 'Homo sapiens')

# 运行时间长，记得保存
save(mt,file ='./rdata/cellcall_result_all.Rdata')
load('./rdata/cellcall_result_all.Rdata')

n <- mt@data$expr_l_r_log2_scale

colnames(n)

names = colnames(n)[c(17:18,20:24)]
names
## 只关心GT，所以需要筛选
## 挑选colnames!!!!
pathway.hyper.list <- lapply(names, function(i){
  print(i)
  tmp <- getHyperPathway(data = n, object = mt, cella_cellb = i)
  return(tmp)
})
myPub.df <- getForBubble(pathway.hyper.list, cella_cellb=names)


# 可视化！-----------------
library(ggplot2)
plotBubble(myPub.df)+scale_colour_gradientn(colours =  c("#6894B9","grey90","#f1707d"))                       


names = colnames(n)[9:16]
names
## 只关心GT，所以需要筛选
## 挑选colnames!!!!
pathway.hyper.list <- lapply(names, function(i){
  print(i)
  tmp <- getHyperPathway(data = n, object = mt, cella_cellb = i)
  return(tmp)
})
myPub.df <- getForBubble(pathway.hyper.list, cella_cellb=names)
library(ggplot2)
plotBubble(myPub.df)+scale_colour_gradientn(colours =  c("#6894B9","grey90","#f1707d"))                       



table(mt@meta.data$celltype)
source("utools.R")
col =c("Fibro"='#A499CC',
       "T cell"='#E0A7C8',
       "Epi"=     '#E069A6',
       "NK"=    "#f1707d",
       "Mac/Mono"=   "#AFC2D9",
       "B cell"=  "#6894B9",
       "Endo"= "#79B99D",
       "SMC/MyoFibro"= "#F5D2A8"
)
col = data.frame(col)
# cellColor <- data.frame(color=ov_palette[1:8], stringsAsFactors = FALSE)
# rownames(cellColor) <- c(names(table(mt@meta.data$celltype)))
pdf("./function/cellcall_circle.pdf")
ViewInterCircos(object = mt, font = 2,cellColor = col,
                lrColor = c("#f1707d", "#6894B9"),
                arr.type = "big.arrow",arr.length = 0.04,
                trackhight1 = 0.05, slot="expr_l_r_log2_scale",
                linkcolor.from.sender = TRUE,
                linkcolor = NULL, gap.degree = 0.5, #细胞类型多的话设置小点，不然图太大画不出来
                
                trackhight2 = 0.032, track.margin2 = c(0.01,0.12), DIY = FALSE)
dev.off()


a=grep('Epi',colnames(n))
b=grep('Mac/Mono',colnames(n))

rt=n[,c(b,a)]
rt=rt[rowSums(rt) >1,]

dev.off()
pheatmap::pheatmap(rt,cluster_rows = F,cluster_cols = F,border_color = NA)+scale_colour_gradientn(colours =  c("#6894B9","grey90","#f1707d"))

mt <- LR2TF(object = mt, sender_cell="Mac/Mono", recevier_cell="Epi",
            slot="expr_l_r_log2_scale", org="Homo sapiens")
head(mt@reductions$sankey)


# if(!require(networkD3)){
#   BiocManager::install("networkD3")
# }

sank <- LRT.Dimplot(mt, 
                    fontSize = 8, 
                    nodeWidth = 30, 
                    height = NULL, 
                    width = 1200,    
                    sinksRight=FALSE, 
                    DIY.color = FALSE)
sank
networkD3::saveNetwork(sank, 
                       "./Mac/Mono-Epi.html")



################################################################################
################################################################################
################################################################################
# step3 scMetabolism
################################################################################
{
  library(Seurat)
  library(scMetabolism)
  library(ggplot2)
  library(rsvd)
}
rm(list = ls())
gc()
load(file = "./rdata/scRNA_estrous.rda")
library(Seurat)
scRNA@meta.data$condition = factor(scRNA@meta.data$condition,levels = c("P","E","M","D"))
table(scRNA@meta.data$condition)
DimPlot(scRNA, group.by = "Level0", 
        split.by = "condition", ncol = 4,pt.size = 0.2)


# devtools::install_github("YosefLab/VISION@v2.1.0")

countexp.Seurat <-scRNA
source("utools.R")
mouse_results <- Mouse.sc.metabolism(countexp.Seurat, metabolism.type = 'KEGG')
# countexp.Seurat@meta.data$cell_type <- countexp.Seurat@active.ident

# countexp.Seurat<-sc.metabolism.Seurat(obj = countexp.Seurat, method = "AUCell",rna_slot = "RNA",
#                                       imputation = F, ncores = 2, metabolism.type = "KEGG")
save(mouse_results,file="rdata/ov_scmetabolism.rda")
load("rdata/ov_scmetabolism.rda")
countexp.Seurat = mouse_results
metabolism.matrix <- countexp.Seurat@assays$METABOLISM$score
metabolism.matrix


library(scMetabolism)
source("utools.R")
library(ggplot2)
p.dimplot<-DimPlot.metabolism(obj = countexp.Seurat, 
                              pathway = "Arachidonic acid metabolism", 
                              dimention.reduction.type = "umap", 
                              dimention.reduction.run = F, size = 1)
p.dimplot

input.pathway <- rownames(countexp.Seurat@assays[["METABOLISM"]][["score"]])[1:85] # 通路序号更改
input.pathway 

gg_table_median_norm = DotPlot_metabolism(obj = countexp.Seurat, 
                             pathway = input.pathway, 
                             phenotype = "condition", 
                             norm = "y")

save(gg_table_median_norm,file = "./rdata/scmetabolism_template.rda")
# head(gg_table_median_norm)
load("./rdata/scmetabolism_template.rda")
{
  gg_table_median_norm$X1 <- factor(gg_table_median_norm$X1, levels = c("P","E","M","D"))

  ####
  ####
  # This function checks if all the records in a subset of df (where subset is based on unique values of 'X2') have 0 on 'X3'
  is_all_zero <- function(x){ all(x == 0) }

  # We group the data by 'X2', apply our function and negate the result
  # This gives us a logical vector that is TRUE for all rows we want to keep and FALSE for those we want to remove
  keep <- !ave(gg_table_median_norm$X3, gg_table_median_norm$X2, FUN = is_all_zero)
  
  # We apply the logical vector to our dataframe
  gg_table_median_norm <- gg_table_median_norm[keep, ]
  gg_table_median_norm = na.omit(gg_table_median_norm)
  library(wesanderson)
  pal <- wes_palette("Zissou1", 100, type = "continuous")
  
  ggplot(data=gg_table_median_norm, aes(x=gg_table_median_norm[,1], y=gg_table_median_norm[,2], color = gg_table_median_norm[,3])) +
    geom_point(data=gg_table_median_norm, aes(size = gg_table_median_norm[,3])) + #geom_line() +
    #theme_bw()+theme(aspect.ratio=0.5, axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Metabolic Pathway")+ xlab(input.parameter)+
    theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1), #aspect.ratio=1,
                     panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    scale_color_gradientn(colours = pal) +
    labs(color = "Value", size = "Value") +
    #facet_wrap(~tissueunique, ncol = 1) +
    #theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    NULL
  # 定义一个新的颜色向量
  # Glycolysis/Gluconeogenesis, Citrate cycle, Pentose phosphate pathway
  gg_table_median_norm$X2
  gg_select = gg_table_median_norm[c(1:12,17:24,37:76,81:112,137:152,165:216),]
  mycolors <- c("#6894B9","grey90","#f1707d")
  
  # 创建一个色板函数
  mypalette <- colorRampPalette(mycolors)
  
  # 使用色板函数生成11个颜色
  pal <- mypalette(100)
  
  # 创建你的散点图
  library(ggplot2)
  ggplot(data=gg_select, aes(x=gg_select[,1], y=gg_select[,2], color = gg_select[,3])) +
    geom_point(data=gg_select, aes(size = gg_select[,3])) +
    ylab('Metabolic Pathway') + xlab("condition")+
    theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1),
                     panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    scale_color_gradientn(colours = pal) +
    labs(color = 'Value', size = 'Value')
  
  


}





p_dotplot <- DotPlot.metabolism(obj = countexp.Seurat, 
                                pathway = input.pathway, 
                                phenotype = "condition", 
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




