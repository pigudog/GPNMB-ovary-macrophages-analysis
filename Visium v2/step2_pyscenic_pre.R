
load("./rdata/scRNA_myelo_gpnmb.rda")
library(SeuratDisk)
library(Seurat) 
# seurat2h5seurat中间过渡        
SaveH5Seurat(scRNA_myelo,filename="pyscenic.h5seurat", overwrite = TRUE)

# 数据转为最终h5ad格式
Convert("pyscenic.h5seurat", dest = "h5ad", overwrite = TRUE) 
