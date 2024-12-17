
library(Seurat)
library(SeuratData)
library(clusterProfiler)
library(org.Mm.eg.db)
rm(list = ls())
gc()
load("./rdata/ov_scRNA_myelo.rda")
# FeaturePlot(scRNA_myelo,features = c("MUC1"))
# 
# scRNA_myelo=subset(scRNA_myelo, Myeloid_subtype != c("DC")) 
# scRNA_myelo=subset(scRNA_myelo, Myeloid_subtype != c("Monocyte")) 
scRNA_myelo@meta.data$Myeloid_subtype = droplevels(scRNA_myelo@meta.data$Myeloid_subtype, exclude = setdiff(levels(scRNA_myelo@meta.data$Myeloid_subtype),unique(scRNA_myelo@meta.data$Myeloid_subtype)))

table(scRNA_myelo@meta.data$Myeloid_subtype)

scRNA_myelo$group=ifelse(as.numeric(scRNA_myelo@assays$RNA@counts['Gpnmb',])>0,'Gpnmb+Mac','Gpnmb-Mac')
table(scRNA_myelo@meta.data$group)

# 使用原有的细胞注释
Idents(scRNA_myelo) <- "group"

# 得到每个细胞亚群的高表达基因
obj.markers <- FindAllMarkers(scRNA_myelo, only.pos = TRUE)
head(obj.markers)

# 使用clusterProfiler进行KEGG Pathway通路富集分析，
# KEGG通路中基因为基因 ENTREZID，需要进行转换；当然如果是GO，直接设置keytype类型就可以了。

group <- data.frame(gene=obj.markers$gene,group=obj.markers$cluster)
Gene_ID <- bitr(obj.markers$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
data  <- merge(Gene_ID,group,by.x='SYMBOL',by.y='gene')

data_GO <- compareCluster(ENTREZID~group, 
                          data=data, 
                          fun="enrichGO", 
                          OrgDb="org.Mm.eg.db",
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 1)

# xx <- compareCluster(gcSample, fun="enrichKEGG", organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
res <- data_GO@compareClusterResult

## 将富集结果中的 ENTREZID 重新转为 SYMBOL
for (i in 1:dim(res)[1]) {
  arr = unlist(strsplit(as.character(res[i,"geneID"]), split="/"))
  gene_names = paste(unique(names(Symbol[Symbol %in% arr])), collapse="/")
  res[i,"geneID"] = gene_names
}

head(res)



## 通路筛选
enrich <- res %>% 
  group_by(Cluster) %>% 
  top_n(n = 5, wt = -pvalue) %>% 
  filter(Cluster %in% c('Gpnmb+Mac','Gpnmb-Mac'))

dt <- enrich
dt <- dt[order(dt$Cluster), ]
dt$Description <- factor(dt$Description, levels = dt$Description)
colnames(dt)


# 先自定义主题：
mytheme <- theme(
  axis.title = element_text(size = 13),
  axis.text = element_text(size = 11),
  axis.text.y = element_blank(), # 在自定义主题中去掉 y 轴通路标签:
  axis.ticks.length.y = unit(0,"cm"),
  plot.title = element_text(size = 13, hjust = 0.5, face = "bold"),
  legend.title = element_text(size = 13),
  legend.text = element_text(size = 11),
  plot.margin = margin(t = 5.5, r = 10, l = 5.5, b = 5.5)
)

p <- ggplot(data = dt, aes(x = -log10(pvalue), y = rev(Description), fill = Cluster)) +
  scale_fill_manual(values =c('#6bb9d2', '#d55640')) +
  geom_bar(stat = "identity", width = 0.5, alpha = 0.8) +
  scale_x_continuous(expand = c(0,0)) + # 调整柱子底部与y轴紧贴
  labs(x = "-Log10(pvalue)", y = "Gpnmb+Mac       Gpnmb-Mac", title = "KEGG Pathway enrichment") +
  # x = 0.61 用数值向量控制文本标签起始位置
  geom_text(size=3.8, aes(x = 0.05, label = Description), hjust = 0) + # hjust = 0,左对齐
  geom_text(size=3, aes(x = 0.05, label = geneID), hjust = 0, vjust = 2.5, color=rep(c('#6bb9d2', '#d55640'),each=5)) + # hjust = 0,左对齐
  theme_classic() + 
  mytheme +
  NoLegend()

p

# 保存，这里的保存宽和高进行了调整，可以使得结果比较美观
ggsave(filename = "test.png", width = 5.2, height = 5.1, plot = p)