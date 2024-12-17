
rm(list = ls())
gc()
# https://saezlab.github.io/liana/articles/liana_nichenet.html
remotes::install_github('saezlab/liana')
# source("ktplots.R")
library(data.table)
library(Seurat)

options(stringsAsFactors = F)
load(file = "../scRNA_estrous/rdata/scRNA_all_with_mac.rda")
table(scRNA@meta.data$Level0)


load(file = "./rdata/scRNA_myelo_gpnmb.rda")
table(scRNA_myelo@meta.data$Level0)
table(scRNA_myelo@meta.data$group)

scRNA_M_other=subset(scRNA,Level0 !='Mac')
table(scRNA_M_other$Level0)

## 记得同等地给一列，与celltype保持一致即可。
scRNA_M_other$group=scRNA_M_other$Level0

table(scRNA_M_other$group)

scRNA_cellcall=merge(scRNA_M_other,scRNA_myelo)
table(scRNA_cellcall$group)
save(scRNA_cellcall,file = "./rdata/scRNA_cellcall.rda")




# cellcall
load("./rdata/scRNA_cellcall.rda")

table(scRNA_cellcall$group)

options()
library(liana)

show_methods()
# [1] "connectome"      "logfc"           "natmi"           "sca"             "cellphonedb"     "cytotalk"        "call_squidpy"    "call_cellchat"  
# [9] "call_connectome" "call_sca"        "call_italk"      "call_natmi" 

# Run liana
liana_test <- liana_wrap(scRNA_cellcall, 
                         method = c("call_cellchat", "cellphonedb","connectome","natmi",  "logfc", "sca","cytotalk"),
                         resource = c("MouseConsensus"),
                         idents_col = "group")
save(liana_test,file = "./rdata/liana_test_raw.rda")
library(liana)
load("./rdata/liana_test_raw.rda")
# Liana returns a list of results, each element of which corresponds to a method
liana_test %>% dplyr::glimpse()
# We can aggregate these results into a tibble with consensus ranks
library(dplyr)
names(liana_test)

liana_test1 <- liana_test[c(2:7)]
# names(liana_test1) = c("call_cellchat", "cellphonedb" ,  "connectome"  ,  "natmi"    ,     "logfc"   ,      "sca"    ,       "cytotalk")

liana_test1 <- liana_test1 %>%
  liana_aggregate(
    aggregate_how = NULL,
    resource = NULL,
    set_cap = "max",
    cap = NULL,
    get_ranks = TRUE,
    get_agrank = TRUE,
    verbose = TRUE,
    join_cols = NULL
  )
table(liana_test1$source)

# dotplot
liana_test2 = liana_test1[liana_test1$receptor.complex == "Gpnmb",]
liana_test2
liana_test1 %>%
  liana_dotplot(source_groups = c("Granulosa",   "Mesenchyme"),
                target_groups = c("Gpnmb-Mac","Gpnmb+Mac"),
                ntop = 50)+scale_colour_gradientn(colours =  c("#6894B9","grey90","#f1707d"))

liana_test1 %>%
  liana_dotplot(source_groups = c("Gpnmb-Mac","Gpnmb+Mac"),
                target_groups = c( "Granulosa",   "Mesenchyme"),
                ntop = 50)+scale_colour_gradientn(colours =  c("#6894B9","grey90","#f1707d"))

liana_test1 %>%
  liana_dotplot(source_groups = c( "Mesenchyme"),
                target_groups = c("Gpnmb+Mac"),
                ntop = 40)+scale_colour_gradientn(colours =  c("#6894B9","grey90","#f1707d"))


liana_test1 %>%
  liana_dotplot(source_groups = c("Gpnmb+Mac"),
                target_groups = c( "Granulosa"),
                ntop = 40)+scale_colour_gradientn(colours =  c("#6894B9","grey90","#f1707d"))


# heatmap
liana_trunc <- liana_test1 %>%
  # only keep interactions concordant between methods
  filter(aggregate_rank <= 0.01) # note that these pvals are already corrected

heat_freq(liana_trunc)

library(circlize)
p <- chord_freq(liana_trunc,
                source_groups = c("Gpnmb-Mac","Gpnmb+Mac"),
                target_groups = c("Granulosa"))
chord_freq(liana_trunc)


###################################################################################
# mac -> Grnulosa & Mesenchyme
# 自定义绘制dotplot
# Modify for the plot
liana_mod <- liana_test1 %>%
  # Filter to only the cells of interest
  `if`(!is.null(c("Gpnmb-Mac","Gpnmb+Mac")),
       filter(., source %in% c("Gpnmb-Mac","Gpnmb+Mac")),
       .) %>%
  `if`(!is.null(c( "Granulosa",   "Mesenchyme")),
       filter(., target %in% c( "Granulosa",   "Mesenchyme")),
       .)
top_int <- liana_mod %>% distinct_at(c("ligand.complex", "receptor.complex")) %>% head(50)
top_int1 = top_int[c(1:4,10,14,17,20,33,35),]
entities <- c("ligand.complex", "receptor.complex")
liana_mod %<>% inner_join(top_int1, by=entities)
magnitude = "sca.LRscore"
specificity = "natmi.edge_specificity"
liana_mod %<>%
  rename(magnitude = !!magnitude) %>%
  rename(specificity = !!specificity) %>%
  unite(entities, col = "interaction", sep = " -> ") %>%
  unite(c("source", "target"), col = "source_target", remove = FALSE)



# ensure levels & order is kept the plot
interactions_order <- liana_mod %>% pull("interaction") %>% unique()
liana_mod %<>%
  mutate(interaction = factor(interaction, levels=rev(interactions_order))) %>%
  mutate(across(where(is.character), as.factor))

# colour blind palette from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
cbPalette <- c("#E69F00", "#56B4E9",
               "#009E73", "#F0E442", "#0072B2",
               "#D55E00", "#CC79A7", "#DF69A7")

# plot
size_range = c(2, 10)
y.label = "Interactions (Ligand -> Receptor)"
colour.label = "Expression\nMagnitude"
size.label = "Interaction\nSpecificity"

p = ggplot(liana_mod,
         aes(x = target,
             y = interaction,
             colour = magnitude,
             size = specificity,
             group = target
         )) +
    geom_point() +
    # scale_color_gradientn(colours = viridis::viridis(20)) +
    scale_colour_gradientn(colours =  c("#6894B9","grey90","#f1707d"))+
    scale_size_continuous(range = size_range) +
    facet_grid(. ~ source,
               space = "free",
               scales ="free",
               switch = "y")  +
    # scale_x_discrete(position = "right") +
    labs(y = y.label,
         colour = colour.label,
         size = size.label,
         x = "Target",
         title= "Source"
    ) +
    theme_bw(base_size = 20) +
    theme(
      legend.text = element_text(size = 16),
      axis.text.x = element_text(colour =
                                   cbPalette[1:length(
                                     unique(liana_mod$source)
                                   )],
                                 face = "bold",
                                 size = 23),
      axis.title.x = element_text(colour = "gray6"),
      axis.text.y = element_text(size = 18,
                                 vjust = 0.5),
      legend.title = element_text(size = 18),
      panel.spacing = unit(0.1, "lines"),
      strip.background = element_rect(fill = NA),
      plot.title = element_text(vjust = 0, hjust=0.5, colour = "gray6"),
      strip.text = element_text(size = 24, colour = "gray6") #,
      # strip.text.y.left = element_text(angle = 0)
    )
  

p


################################################################################
# mac -> Grnulosa & Mesenchyme
# 自定义绘制dotplot
# Modify for the plot
liana_mod <- liana_test1 %>%
  # Filter to only the cells of interest
  `if`(!is.null(c("Granulosa",   "Mesenchyme")),
       filter(., source %in% c("Granulosa",   "Mesenchyme")),
       .) %>%
  `if`(!is.null(c("Gpnmb-Mac","Gpnmb+Mac")),
       filter(., target %in% c( "Gpnmb-Mac","Gpnmb+Mac")),
       .)
top_int <- liana_mod %>% distinct_at(c("ligand.complex", "receptor.complex")) %>% head(50)
top_int1 = top_int[c(2:3,8, 19:21,31,36,38,42),]
entities <- c("ligand.complex", "receptor.complex")
liana_mod %<>% inner_join(top_int1, by=entities)
magnitude = "sca.LRscore"
specificity = "natmi.edge_specificity"
liana_mod %<>%
  rename(magnitude = !!magnitude) %>%
  rename(specificity = !!specificity) %>%
  unite(entities, col = "interaction", sep = " -> ") %>%
  unite(c("source", "target"), col = "source_target", remove = FALSE)



# ensure levels & order is kept the plot
interactions_order <- liana_mod %>% pull("interaction") %>% unique()
liana_mod %<>%
  mutate(interaction = factor(interaction, levels=rev(interactions_order))) %>%
  mutate(across(where(is.character), as.factor))

# colour blind palette from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
cbPalette <- c("#E69F00", "#56B4E9",
               "#009E73", "#F0E442", "#0072B2",
               "#D55E00", "#CC79A7", "#DF69A7")

# plot
size_range = c(2, 10)
y.label = "Interactions (Ligand -> Receptor)"
colour.label = "Expression\nMagnitude"
size.label = "Interaction\nSpecificity"

p = ggplot(liana_mod,
           aes(x = target,
               y = interaction,
               colour = magnitude,
               size = specificity,
               group = target
           )) +
  geom_point() +
  # scale_color_gradientn(colours = viridis::viridis(20)) +
  scale_colour_gradientn(colours =  c("#6894B9","grey90","#f1707d"))+
  scale_size_continuous(range = size_range) +
  facet_grid(. ~ source,
             space = "free",
             scales ="free",
             switch = "y")  +
  # scale_x_discrete(position = "right") +
  labs(y = y.label,
       colour = colour.label,
       size = size.label,
       x = "Target",
       title= "Source"
  ) +
  theme_bw(base_size = 20) +
  theme(
    legend.text = element_text(size = 16),
    axis.text.x = element_text(colour =
                                 cbPalette[1:length(
                                   unique(liana_mod$source)
                                 )],
                               face = "bold",
                               size = 23),
    axis.title.x = element_text(colour = "gray6"),
    axis.text.y = element_text(size = 18,
                               vjust = 0.5),
    legend.title = element_text(size = 18),
    panel.spacing = unit(0.1, "lines"),
    strip.background = element_rect(fill = NA),
    plot.title = element_text(vjust = 0, hjust=0.5, colour = "gray6"),
    strip.text = element_text(size = 24, colour = "gray6") #,
    # strip.text.y.left = element_text(angle = 0)
  )


p


if(!require('nichenetr')) remotes::install_github("saeyslab/nichenetr", quiet = TRUE)
library(tidyverse)
library(liana)
library(nichenetr)
library(Seurat)
library(ggrepel)
library(cowplot)
options(timeout=600) # required to download expression data /w slow connection

load("./rdata/liana_test_raw.rda")
liana_results <- liana_test[c(2:7)] %>%
  liana_aggregate()
# filter results to cell types of interest
cam_tumor_results <- liana_results %>%
  subset(source == "Gpnmb+Mac" & target == "Granulosa") %>%
  dplyr::rename(ligand=ligand.complex, receptor=receptor.complex)

# filter results to top N interactions
n <- 30
top_n_caf_tumor <- cam_tumor_results %>%
  arrange(aggregate_rank) %>%
  slice_head(n = n) %>%
  mutate(id = fct_inorder(paste0(ligand, " -> ", receptor)))

# visualize median rank
top_n_caf_tumor %>%
  ggplot(aes(y = aggregate_rank, x = id)) +
  geom_bar(stat = "identity", fill = "#f1707d") + # <--- 将颜色修改为 "#FF5733",这是一个橙色的十六进制颜色代码，你可以修改为你喜欢的颜色
  xlab("Interaction") + ylab("LIANA's aggregate rank") +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 8, angle = 60, hjust = 1, vjust = 1))






load("./rdata/scRNA_cellcall.rda")
scRNA_cellcall <- scRNA_cellcall %>% subset(group %in% c("Gpnmb+Mac", "Granulosa"))




expression <- t(as.matrix(scRNA_cellcall@assays$RNA@data))
sample_info <- scRNA_cellcall@meta.data
colnames(sample_info) <- make.names(colnames(sample_info))

# filter samples based on vignette's information and add cell type
sample_info <- sample_info %>%
  # subset( !(tumor %in% tumors_remove) & Lymph.node == 0) %>%
  # fix some cell type identity names
  # mutate(cell_type = ifelse(classified..as.cancer.cell == 1, "Tumor", non.cancer.cell.type)) %>%
  subset(group %in% c("Gpnmb+Mac", "Granulosa"))

# cell ID as rownames
sample_info$cell <- rownames(sample_info) 

# subset expression to selected cells
expression <- expression[sample_info$cell, ]

# gene set of interest
library(data.table)
library(dplyr)
# 读取配体靶基因矩阵 https://zenodo.org/record/3260758/files/ligand_target_matrix.rds
ligand_target_matrix = readRDS("./ligand_target_matrix.rds")
# Convert the ligand-target model from human to mouse symbols.

# Because not all human genes have a mouse one-to-one ortholog, these genes will be removed from the mouse model.

colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols() 
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols() 

ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]

dim(ligand_target_matrix)




# geneset_oi <- fread("./scenscence_marker.csv")  %>%
#   pull(V1) %>%
#   .[. %in% rownames(ligand_target_matrix)]
geneset_oi <- fread("./files/gmt/GOBP_OVARIAN_FOLLICLE_DEVELOPMENT.v2023.2.Hs.gmt",header = F)
geneset_oi = geneset_oi[,3:ncol(geneset_oi)]
geneset_oi = t(geneset_oi)
geneset_oi = as.vector(geneset_oi)
geneset_oi = geneset_oi  %>% convert_human_to_mouse_symbols()
geneset_oi <-  geneset_oi %>%
  .[. %in% rownames(ligand_target_matrix)]
background_genes <- expression %>%
  apply(2,function(x){10*(2**x - 1)}) %>%
  apply(2,function(x){log2(mean(x) + 1)}) %>%
  .[. >= 4] %>%
  names()

# get ligands and filter to those included in NicheNet's ligand-target matrix
ligands <- unique(top_n_caf_tumor$ligand)
ligands <- ligands[ligands %in% colnames(ligand_target_matrix)]
ligands

nichenet_activities <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = background_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = ligands
)


# prepare data for visualization
vis_liana_nichenet <- top_n_caf_tumor %>%
  inner_join(nichenet_activities, by = c("ligand" = "test_ligand")) %>%
  arrange(pearson) %>%
  mutate(ligand = fct_inorder(ligand))

# prepare NicheNet figure
nichenet_scores_plot <- vis_liana_nichenet %>%
  group_by(ligand) %>%
  summarize(pearson = mean(pearson)) %>%
  # mutate(pearson_discrete = cut(pearson, breaks = 5)) %>%  # 离散化pearson值
  ggplot(aes(y = ligand, x = pearson, fill = pearson)) +  # 将fill映射到新的离散变量
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "white", high = "#f1707d") +  # 设置颜色渐变
  ggtitle("NicheNet") +
  xlab("Pearson's score") +
  theme_cowplot() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_line(color = "white"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))
# prepare LIANA figure
liana_receptor_heatmap <- vis_liana_nichenet %>%
  ggplot(aes(y = ligand, x = receptor, fill = aggregate_rank)) +
  geom_tile() +
  scale_fill_gradient(low = "grey90", high = "#f1707d") +  # 这行代码定义了颜色渐变
  theme_cowplot() +
  ggtitle("LIANA") +
  ylab("Ligand") + xlab("Receptor") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_line(colour = "gray", linetype = 2),
        legend.position = "left")

# combine plots
plot_grid(liana_receptor_heatmap, nichenet_scores_plot,
          align = "h", nrow = 1, rel_widths = c(0.8,0.3))
