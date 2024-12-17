

################################################################################
# 1. Anndata -> Seurat (V4)
################################################################################
library(scPDtools)

# AnnData to Seurat
h5ad_file = './pydata/adata_ovary_combined_raw_sample_follicle.h5ad'
convertFormat(h5ad_file, from="anndata", to="seurat",main_layer = "counts_log1p",
              outFile='rdata/adata_ovary_combined_raw_sample_follicle.rds')

library(Seurat)
adata = readRDS("rdata/adata_ovary_combined_raw_sample_follicle.rds")

save(adata,file="rdata/adata_ovary_combined_raw_sample_follicle.rda")

table(adata@meta.data$celltype_annotation)
adata_fibro = subset(adata,celltype_annotation %in% c("M_Fibroblast-like_cells"))
save(adata_fibro,file="rdata/adata_fibro.rda")
adata_oocyte = subset(adata,celltype_annotation %in% c("Oocyte"))
save(adata_oocyte,file="rdata/adata_oocyte.rda")
adata_immune = subset(adata,celltype_annotation %in% c("Immune"))
save(adata_immune,file="rdata/adata_immune.rda")
# adata_GC = subset(adata,celltype_annotation %in% c("GC_Antral",""))
# save(adata_immune,file="rdata/adata_immune.rda")


library(scPDtools)
# AnnData to Seurat
h5ad_file = './pydata/adata_antral_gpnmb_to_r.h5ad'
convertFormat(h5ad_file, from="anndata", to="seurat",main_layer = "counts_log1p",
              outFile='rdata/adata_antral_gpnmb_to_r.rds')

library(Seurat)
adata = readRDS("rdata/adata_antral_gpnmb_to_r.rds")

save(adata,file="rdata/adata_antral_gpnmb_to_r.rda")



library(scPDtools)
# AnnData to Seurat
h5ad_file = './pydata/granulosa_pseudotime_true.h5ad'
convertFormat(h5ad_file, from="anndata", to="seurat",main_layer = "counts_log1p",
              outFile='rdata/granulosa_pseudotime_true.rds')

library(Seurat)
adata = readRDS("rdata/granulosa_pseudotime_true.rds")

save(adata,file="rdata/granulosa_pseudotime_true.rda")




