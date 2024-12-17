
library(Seurat)
# library(nichenetr)
library(dplyr)
# library(nichenetr)
# library(dplyr)
mapper = function(df, value_col, name_col) setNames(df[[value_col]], df[[name_col]])
# 基因转化
convert_mouse_to_human_symbols = function(symbols, version = 1){
  
  if(!is.character(symbols))
    stop("symbols should be a character vector of mouse gene symbols")
  
  requireNamespace("dplyr")
  requireNamespace("tidyverse")
  load("./important_tmp/geneinfo_human.rda")
  if(version == 1){
    unambiguous_mouse_genes = geneinfo_human %>% 
      na.omit() %>% 
      group_by(symbol_mouse) %>% 
      count() %>% 
      filter(n < 2) %>% 
      .$symbol_mouse
    
    ambiguous_mouse_genes = geneinfo_human %>% 
      na.omit() %>% 
      group_by(symbol_mouse) %>% 
      count() %>% 
      filter(n >= 2) %>% 
      .$symbol_mouse
    
    geneinfo_ambiguous_solved = geneinfo_human %>% 
      filter(symbol_mouse %in% ambiguous_mouse_genes) %>% 
      filter(symbol == toupper(symbol_mouse))
    
    geneinfo_human = geneinfo_human %>% 
      filter(symbol_mouse %in% unambiguous_mouse_genes) %>% 
      bind_rows(geneinfo_ambiguous_solved) %>% 
      na.omit()
    
    mousesymbol2humansymbol = mapper(geneinfo_human, "symbol", "symbol_mouse")
    
    converted_symbols = symbols %>% 
      mousesymbol2humansymbol[.]
  } else if(version == 2) {
    unambiguous_mouse_genes = geneinfo_2022 %>% 
      na.omit() %>% 
      group_by(symbol_mouse) %>% 
      count() %>% 
      filter(n < 2) %>% 
      .$symbol_mouse
    
    ambiguous_mouse_genes = geneinfo_2022 %>% 
      na.omit() %>% 
      group_by(symbol_mouse) %>% 
      count() %>% 
      filter(n >= 2) %>% 
      .$symbol_mouse
    
    geneinfo_ambiguous_solved = geneinfo_2022 %>% 
      filter(symbol_mouse %in% ambiguous_mouse_genes) %>% 
      filter(symbol == toupper(symbol_mouse))
    
    geneinfo_2022 = geneinfo_2022 %>% 
      filter(symbol_mouse %in% unambiguous_mouse_genes) %>% 
      bind_rows(geneinfo_ambiguous_solved) %>% 
      na.omit()
    
    mousesymbol2humansymbol = mapper(geneinfo_2022, "symbol", "symbol_mouse")
    
    converted_symbols = symbols %>% 
      mousesymbol2humansymbol[.]
  }
  
  
  return(converted_symbols)
}

# 小鼠的scmetabolism
Mouse.sc.metabolism <- function(obj,
                                metabolism.type=c("KEGG","REACTOME")){
  #将鼠的基因名转化为人的
  # obj = countexp.Seurat
  # load("./important_tmp/biomaRt_tmp.rda")
  # library(biomaRt)
  
  # human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
  # mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  
  # obj = countexp.Seurat
  load("./important_tmp/geneinfo_human.rda")
  # source("convert.R")
  # mice.gene <- c("Sem1" ,"Gm42418" ,"Gm9843" ,"Rack1", "Rpl13a" ,"mt-Nd1" ,"Selenok",        "Rpl23a", "Rpl23a-ps3",'Kap')
  # human.gene <- c("HSBP1","INPP5F","ABCD2","AURKAIP1","ANP32E","AL138778.1","FEZ2",      
  #                 "PCYOX1","CHROMR","RPL24")
  gene_trans <- convert_mouse_to_human_symbols(symbols = rownames(obj))
  # result2 <- rownames(obj)
  gene_trans <- cbind(names(gene_trans),as.character(gene_trans)) %>% na.omit
  
  library(dplyr)
  # gene_trans$MGI.symbol %>% unique()
  # gene_trans$HGNC.symbol %>% unique()
  # colnames(gene_trans) <- c('mouse','human')
  gene_trans = as.data.frame(gene_trans)
  colnames(gene_trans) <- c('mouse','human')
  gene_trans = gene_trans[!duplicated(gene_trans$mouse),]
  gene_trans = gene_trans[!duplicated(gene_trans$human),]
  mouse_data_trans <- subset(obj,features=gene_trans$mouse)
  #函数参考
  #https://www.jianshu.com/p/6495706bac53
  RenameGenesSeurat <- function(obj,newnames,gene.use=NULL,de.assay) {
    # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration.
    # It only changes obj@assays$RNA@counts, @data and @scale.data.
    print("Run this before integration. It only changes obj@assays$*@counts, @data and @scale.data, @var.features,@reductions$pca@feature.loadings")
    # obj = mouse_data_trans
    lassays <- Assays(obj)
    #names(obj@assays)
    # assay.use <- obj@reductions$pca@assay.used
    DefaultAssay(obj) <- de.assay
    if (is.null(gene.use)) {
      all_genenames <- rownames(obj)
    }else{
      all_genenames <- gene.use
      obj <- subset(obj,features=gene.use)
    }
    
    order_name <- function(v1,v2,ref){
      v2 <- make.names(v2,unique=T)
      df1 <- data.frame(v1,v2)
      rownames(df1) <- df1$v1
      df1 <- df1[ref,]
      return(df1)
    }
    
    df1 <- order_name(v1=all_genenames,v2=newnames,ref=rownames(obj))
    all_genenames <- df1$v1
    newnames <- df1$v2
    
    if ('SCT' %in% lassays) {
      if ('SCTModel.list' %in%  slotNames(obj@assays$SCT)) {
        obj@assays$SCT@SCTModel.list$model1@feature.attributes <- obj@assays$SCT@SCTModel.list$model1@feature.attributes[all_genenames,]
        rownames(obj@assays$SCT@SCTModel.list$model1@feature.attributes) <- newnames
      }
    }
    change_assay <- function(a1=de.assay,obj,newnames=NULL,all_genenames=NULL){
      RNA <- obj@assays[a1][[1]]
      if (nrow(RNA) == length(newnames)) {
        if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
        if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
        if (length(RNA@var.features)) {
          df1 <- order_name(v1=all_genenames,v2=newnames,ref=RNA@var.features)
          all_genenames1 <- df1$v1
          newnames1 <- df1$v2
          RNA@var.features            <- newnames1
        }
        if (length(RNA@scale.data)){
          df1 <- order_name(v1=all_genenames,v2=newnames,ref=rownames(RNA@scale.data))
          all_genenames1 <- df1$v1
          newnames1 <- df1$v2
          rownames(RNA@scale.data)    <- newnames1
        }
        
      } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
      obj@assays[a1][[1]] <- RNA
      return(obj)
    }
    
    for (a in lassays) {
      DefaultAssay(obj) <- a
      df1 <- order_name(v1=all_genenames,v2=newnames,ref=rownames(obj))
      all_genenames1 <- df1$v1
      newnames1 <- df1$v2
      obj <- change_assay(obj=obj,a1=a,newnames=newnames1,all_genenames=all_genenames1)
    }
    assay.use = de.assay
    hvg <- VariableFeatures(obj,assay=assay.use)
    if (length(obj@reductions$pca)){
      df1 <- order_name(v1=all_genenames,v2=newnames,ref=hvg)
      df1 <- df1[rownames(obj@reductions$pca@feature.loadings),]
      all_genenames1 <- df1$v1
      newnames1 <- df1$v2
      rownames(obj@reductions$pca@feature.loadings) <- newnames1
    }
    try(obj[[de.assay]]@meta.features <- data.frame(row.names = rownames(obj[[de.assay]])))
    return(obj)
  }
  #转化
  mouse_data_trans <- RenameGenesSeurat(mouse_data_trans, 
                                        newnames = gene_trans$human,
                                        gene.use = gene_trans$mouse,
                                        de.assay = 'RNA')
  

  return(mouse_data_trans)
}




# mouse_results <- Mouse.sc.metabolism(mouse_data, metabolism.type = 'KEGG')


sc.metabolism.Seurat <- function(obj, method = "VISION", imputation = F, ncores = 2, metabolism.type = "KEGG",rna_slot="SCT") {
  if (rna_slot=="SCT"){
    countexp<-obj@assays$SCT@counts # RNA or SCT
  }else{
    countexp<-obj@assays$RNA@counts 
  }
  
  
  countexp<-data.frame(as.matrix(countexp))
  
  signatures_KEGG_metab <- "./important_tmp//KEGG_metabolism_nc.gmt"
  signatures_REACTOME_metab <- "./important_tmp/REACTOME_metabolism.gmt"
  
  # signatures_KEGG_metab <- system.file("data", "KEGG_metabolism_nc.gmt", package = "scMetabolism")
  # signatures_REACTOME_metab <- system.file("data", "REACTOME_metabolism.gmt", package = "scMetabolism")
  
  
  if (metabolism.type == "KEGG")  {gmtFile<-signatures_KEGG_metab; cat("Your choice is: KEGG\n")}
  if (metabolism.type == "REACTOME")  {gmtFile<-signatures_REACTOME_metab; cat("Your choice is: REACTOME\n")}
  
  #imputation
  if (imputation == F) {
    countexp2<-countexp
  }
  if (imputation == T) {
    
    cat("Start imputation...\n")
    
    #Citation: George C. Linderman, Jun Zhao, Yuval Kluger. Zero-preserving imputation of scRNA-seq data using low-rank approximation. bioRxiv. doi: https://doi.org/10.1101/397588
    #Github: https://github.com/KlugerLab/ALRA
    
    cat("Citation: George C. Linderman, Jun Zhao, Yuval Kluger. Zero-preserving imputation of scRNA-seq data using low-rank approximation. bioRxiv. doi: https://doi.org/10.1101/397588 \n")
    
    
    result.completed <- alra(as.matrix(countexp))
    countexp2 <- result.completed[[3]]; row.names(countexp2) <- row.names(countexp)
  }
  
  #signature method
  cat("Start quantify the metabolism activity...\n")
  
  #VISION
  if (method == "VISION") {
    library(VISION)
    n.umi <- colSums(countexp2)
    scaled_counts <- t(t(countexp2) / n.umi) * median(n.umi)
    vis <- Vision(scaled_counts, signatures = gmtFile)
    
    options(mc.cores = ncores)
    
    vis <- analyze(vis)
    
    signature_exp<-data.frame(t(vis@SigScores))
  }
  
  #AUCell
  if (method == "AUCell") {
    library(AUCell)
    library(GSEABase)
    cells_rankings <- AUCell_buildRankings(as.matrix(countexp2), nCores=ncores, plotStats=F) #rank
    geneSets <- getGmt(gmtFile) #signature read
    cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings) #calc
    signature_exp <- data.frame(getAUC(cells_AUC))
  }
  
  #ssGSEA
  if (method == "ssGSEA") {
    library(GSVA)
    library(GSEABase)
    geneSets <- getGmt(gmtFile) #signature read
    gsva_es <- gsva(as.matrix(countexp2), geneSets, method=c("ssgsea"), kcdf=c("Poisson"), parallel.sz=ncores) #
    signature_exp<-data.frame(gsva_es)
  }
  
  #GSVA
  if (method == "ssGSEA") {
    library(GSVA)
    library(GSEABase)
    geneSets <- getGmt(gmtFile) #signature read
    gsva_es <- gsva(as.matrix(countexp2), geneSets, method=c("gsva"), kcdf=c("Poisson"), parallel.sz=ncores) #
    signature_exp<-data.frame(gsva_es)
  }
  
  cat("\nPlease Cite: \nYingcheng Wu, Qiang Gao, et al. Cancer Discovery. 2021. \nhttps://pubmed.ncbi.nlm.nih.gov/34417225/   \n\n")
  
  obj@assays$METABOLISM$score<-signature_exp
  obj
}


DotPlot_metabolism <- function(obj, 
                               pathway, 
                               phenotype, 
                               norm = "y", 
                               calc_type = "median") {
  input.norm = norm
  input.pathway <- pathway
  input.parameter <- phenotype
  
  metadata <- obj@meta.data
  metabolism.matrix <- obj@assays$METABOLISM$score
  
  cat("\nPlease Cite: \nYingcheng Wu, Qiang Gao, et al. Cancer Discovery. 2021. \nhttps://pubmed.ncbi.nlm.nih.gov/34417225/   \n\n")
  
  metadata[, input.parameter] <- as.character(metadata[, input.parameter])
  metabolism.matrix_sub <- t(metabolism.matrix[input.pathway, ])
  
  # Arrange large table
  gg_table <- c()
  for (i in 1:length(input.pathway)) {
    gg_table <- rbind(gg_table, cbind(metadata[, input.parameter], input.pathway[i], metabolism.matrix_sub[, i]))
  }
  gg_table <- data.frame(gg_table)
  
  # Calculate median or average
  gg_table_stat <- c()
  input.group.x <- unique(as.character(gg_table[, 1])) # condition
  input.group.y <- unique(as.character(gg_table[, 2])) # pathway
  
  for (x in 1:length(input.group.x)) {
    for (y in 1:length(input.group.y)) {
      gg_table_sub <- subset(gg_table, gg_table[, 1] == input.group.x[x] & gg_table[, 2] == input.group.y[y])
      
      if (calc_type == "median") {
        stat_value <- median(as.numeric(as.character(gg_table_sub[, 3])), na.rm = TRUE)
      } else if (calc_type == "average") {
        stat_value <- mean(as.numeric(as.character(gg_table_sub[, 3])), na.rm = TRUE)
      }
      
      gg_table_stat <- rbind(gg_table_stat, cbind(input.group.x[x], input.group.y[y], stat_value))
    }
  }
  # print(gg_table_stat)
  gg_table_stat <- data.frame(gg_table_stat)
  gg_table_stat[, 3] <- as.numeric(as.character(gg_table_stat[, 3]))
  
  #normalize
  gg_table_stat_norm<-c()
  # input.group.x<-unique(as.character(gg_table[,1]))
  # input.group.y<-unique(as.character(gg_table[,2]))
  
  
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  
  if (input.norm == "y") {
    for (y in 1:length(input.group.y)) {
      gg_table_stat_sub <- subset(gg_table_stat, gg_table_stat[, 2] == input.group.y[y]) # pathway
      norm_value <- range01(as.numeric(as.character(gg_table_stat_sub[, 3])))
      gg_table_stat_sub[, 3] <- norm_value
      gg_table_stat_norm <- rbind(gg_table_stat_norm, gg_table_stat_sub)
    }
  } else if (input.norm == "x") {
    for (x in 1:length(input.group.x)) {
      gg_table_stat_sub <- subset(gg_table_stat, gg_table_stat[, 1] == input.group.x[x])
      norm_value <- range01(as.numeric(as.character(gg_table_stat_sub[, 3])))
      gg_table_stat_sub[, 3] <- norm_value
      gg_table_stat_norm <- rbind(gg_table_stat_norm, gg_table_stat_sub)
    }
  } else if (input.norm == "na") {
    gg_table_stat_norm <- gg_table_stat
  }
  # print(gg_table_stat_norm)
  gg_table_stat_norm <- data.frame(gg_table_stat_norm)
  gg_table_stat_norm[, 3] <- as.numeric(as.character(gg_table_stat_norm[, 3]))
  
  colnames(gg_table_stat_norm) = c("Condition","Pathway","Value")
  # 假设你的数据框为 df，并且包含 Condition, Pathway, 和 value1 列
  gg_table_stat_norm <- gg_table_stat_norm %>%
    group_by(Pathway) %>%
    # 如果某个 Pathway 中有任意一个 NA，就过滤掉整个 Pathway
    filter(!any(is.na(Value))) %>%
    ungroup()
  
  
  return(gg_table_stat_norm)
}


