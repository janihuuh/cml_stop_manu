

getClusterPhenotypesNK <- function(clusters){
  
  plyr::revalue(clusters, replace = c("0"	 = "0 CD56dim 1\nHAVCR2+ TIGIT+",
                                      "1"	 = "1 CD56dim 2\nLAG3+ FCGR3A+ PRF1+",
                                      "2"	 = "2 CD56dim 3\nNKG7+",
                                      "3"	 = "3 CD56dim 4\nKLRB1+ KLRG1+",
                                      "4"	 = "4 Adaptive",
                                      "5"	 = "5 CD56 bright",
                                      "6"	 = "6 CD56dim 5\nIFNg+",
                                      "7"	 = "7 CD56dim 6\nMALAT1"))
}

getQC <- function(seurat_object){
  
  ###################
  
  min_mito     <- 0
  max_mito     <- 15
  
  min_ribo     <- 10
  max_ribo     <- 50
  
  min_features <- 250
  max_features <- 4500
  
  min_counts   <- 1000
  max_counts   <- 20e3
  
  cycle.genes  <- c("ANLN", "ASPM","BIRC5","CCNA2","CCNB1","CCNB2","CCND1","CD63","CDC20","CDCA8","CDKN3","CENPE","CENPF",
                    "CEP55","CKAP2L","DLGAP5","FOXM1","GTSE1","H2AFZ","HIST1H1B", "HIST1H1C", "HIST1H1D", "HIST1H1E", "HIST1H2AJ",
                    "HIST1H4C", "HJURP", "HMGB1", "HMGB2", "HMMR", "KIF11", "KIF14", "KIF15", "KIF2C", "LMNA",
                    "MCM3", "MKI67", "NCAPG", "NUSAP1", "PCNA", "PLK1", "PRC1", "RRM2", "SMC4", "STMN1", "TK1", "TOP2A", "TPX2", "TUBA1B",
                    "TUBB", "TYMS", "UBE2C")
  
  ###################
  
  seurat_object  <- PercentageFeatureSet(seurat_object, pattern = "^MT-", col.name = "percent.mt")
  seurat_object  <- PercentageFeatureSet(seurat_object, pattern = "^RP", col.name = "percent.ribo")
  # seurat_object  <- PercentageFeatureSet(seurat_object, features = cycle.genes, col.name = "percent.cycle")
  seurat_object@meta.data$barcode <- colnames(seurat_object)
  
  ## In total, we remove with the following conditions:
  qc_df <- seurat_object@meta.data %>% as.data.frame()
  
  # VlnPlot(seurat_object, group.by = "Patient_ID", features = "percent.mt") # + geom_hline(yintercept = 4500) + geom_hline(yintercept = 100)
  # VlnPlot(seurat_object, group.by = "Patient_ID", features = "percent.ribo") # + geom_hline(yintercept = 4500) + geom_hline(yintercept = 100)
  # VlnPlot(seurat_object, group.by = "Patient_ID", features = "nFeature_RNA") # + geom_hline(yintercept = 4500) + geom_hline(yintercept = 100)
  # VlnPlot(seurat_object, group.by = "Patient_ID", features = "nCount_RNA") + geom_hline(yintercept = 30e3) + geom_hline(yintercept = 5e2)
  
  
  percent_mito_outlier <- qc_df %>% dplyr::filter(percent.mt   > max_mito     | percent.mt   < min_mito)     %>% pull(barcode) %>% as.character()
  percent_ribo_outlier <- qc_df %>% dplyr::filter(percent.ribo > max_ribo     | percent.ribo < min_ribo)     %>% pull(barcode) %>% as.character()
  features_outlier     <- qc_df %>% dplyr::filter(nFeature_RNA < min_features | nFeature_RNA > max_features) %>% pull(barcode) %>% as.character()
  umis_outlier         <- qc_df %>% dplyr::filter(nCount_RNA   > max_counts   | nCount_RNA   < min_counts)   %>% pull(barcode) %>% as.character()
  
  outlier_cells        <- c(percent_mito_outlier,
                            percent_ribo_outlier,
                            features_outlier,
                            umis_outlier)
  
  reason               <- c(rep("percent_mito_outlier", length(percent_mito_outlier)),
                            rep("percent_ribo_outlier", length(percent_ribo_outlier)),
                            rep("features_outlier",     length(features_outlier)),
                            rep("umis_outlier",         length(umis_outlier)))
  
  outlier_df <- data.frame(barcode = outlier_cells, reason = reason) %>% dplyr::mutate(from = extractName(barcode)) #, 1, 10))
  
  ## Remove the cells from Seurat-object and save a new seurat-object
  cells.to.use  <- colnames(seurat_object)[!colnames(seurat_object) %in% outlier_df$barcode]
  seurat_object <- subset(seurat_object, cells = cells.to.use)
  return(seurat_object)
  
}



getSingler <- function(seurat_object, folder){

  ### Run singler on public
  ## Try SingleR
  dir.create(folder, showWarnings = F)

  for(patient in unique(seurat_object$orig.ident)){

    require(SingleR)
    message(patient)

    ## Cells to select
    cells.to.keep <- seurat_object@meta.data[seurat_object$orig.ident == patient, ] %>% pull(barcode)
    public_temp      = subset(seurat_object, cells = cells.to.keep)

    annot = data.frame(public_temp@meta.data)
    rownames(annot) <- colnames(public_temp)

    singler = SingleR::CreateSinglerObject(public_temp@assays$RNA@counts,
                                           project.name = patient,
                                           min.genes = 0,
                                           technology = "10X",
                                           species = "Human",
                                           do.signatures = F,
                                           fine.tune = F,
                                           clusters = Idents(public_temp))

    saveRDS(singler, file = paste0(folder, patient, '.rds'))

  }


  ## Combine the data sets
  singler.objects.file <- list.files(folder,  full.names = T)
  singler.objects.file <- grep(".rds", singler.objects.file, value = T)
  singler.objects      <- lapply(singler.objects.file, FUN = function(x) {message(x); readRDS(x)})

  singler = SingleR.Combine(singler.objects,
                            order    = colnames(seurat_object),
                            clusters = Idents(seurat_object),
                            xy       = seurat_object@reductions$umap@cell.embeddings)

  metadata <- data.frame(barcode        = names(singler$meta.data$orig.ident),
                         hpca_pred      = singler$singler[[1]]$SingleR.single$labels,
                         blueprint_pred = singler$singler[[2]]$SingleR.single$labels) %>% add_rownames(var = "barcode")

  public_metadata <- merge(seurat_object@meta.data, metadata, by = "barcode")


  ## Add into Seurat
  public_metadata                      <- public_metadata[match(colnames(seurat_object), public_metadata$barcode), ]
  seurat_object$singler_hpca_pred      <- public_metadata$hpca_pred
  seurat_object$singler_blueprint_pred <- public_metadata$blueprint_pred

  return(seurat_object)
}


## Run seurat
preprocessSeurat <- function(orig_object, cells.to.use){

  ## Subset object
  object <- subset(orig_object, cells = cells.to.use)

  # orig_object@meta.data$barcode
  temp_meta <- orig_object@meta.data[as.character(orig_object@meta.data$barcode) %in% cells.to.use, ]
  temp_meta <- temp_meta[match(colnames(object), temp_meta$barcode), ]
  temp_meta$barcode == colnames(object)
  object@meta.data <- temp_meta

  ## Normalize and find HVGs
  object  <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
  object  <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000, clip.max = 10)

  ## Remove clonality genes
  hvg     <- VariableFeatures(object)
  too_hvg <- HVFInfo(object = object) %>% add_rownames(var = "gene") %>% filter(variance.standardized > 10) %>% pull("gene") %>% as.character()
  hvg     <- hvg[!hvg %in% too_hvg]
  hvg     <- hvg[!hvg %in% clonality_genes]
  hvg     <- hvg[!hvg %in% unwanted_genes]

  VariableFeatures(object) <- hvg
  # plotHVG(object, 30) #+ ylim(values = c(0,10))

  ## Scale data
  object <- ScaleData(object, features = hvg)

  ## PCA data
  object <- RunPCA(object, features = hvg, npcs = 50)
  nPCs   <- sum(object[["pca"]]@stdev > 2)
  print(paste("nPCs:", nPCs))

  ## RunUMAP does not work
  object <- RunUMAP(object, dims = 1:nPCs, learning.rate = 1)

  # Meanwhile try something hacky-ish
  # umap_df <- object[["pca"]]@cell.embeddings[,1:nPCs] %>% umapr::umap() %>% select(UMAP1:UMAP2)
  # umap_df <- CreateDimReducObject(key = "umap", embeddings = as.matrix(x = umap_df))
  # object[["umap"]] <- umap_df

  return(object)

}


plotSlingshot <- function(slingshot_object, reducedDim){

  colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)

  cololors                <- slingshot_object$orig.clusters %>% extractClusterNumber() %>% as.numeric()
  cololors[cololors == 0] <- getPalette5(8)[1]
  cololors[cololors == 1] <- getPalette5(8)[2]
  cololors[cololors == 2] <- getPalette5(8)[3]
  cololors[cololors == 3] <- getPalette5(8)[4]
  cololors[cololors == 4] <- getPalette5(8)[5]
  cololors[cololors == 5] <- getPalette5(8)[6]
  cololors[cololors == 6] <- getPalette5(8)[7]
  cololors[cololors == 7] <- getPalette5(8)[8]


  curve_cols <- c("black", "darkred", "darkblue", "darkolivegreen4")

  colors2 <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(9)

  sling_curve <- SlingshotDataSet(slingshot_object)

  par(mfrow = c(1,2))
  plot(reducedDims(slingshot_object)[[reducedDim]], col = cololors, pch = 16, asp = 1)
  for(i in 1:2){
    lines(sling_curve@curves[[i]], lwd = 4, col = curve_cols[i])
  }

  plot(reducedDims(slingshot_object)[[reducedDim]], col = colors[cut(slingshot_object$slingPseudotime_1,breaks=100)], pch=16, asp = 1)
  # lines(SlingshotDataSet(slingshot_object), lwd = 2, col = "black")
  for(i in 1:2){
    lines(sling_curve@curves[[i]], lwd = 4, col = curve_cols[i])
  }

  par(mfrow = c(1,1))


}


runSlingshot <- function(sce_object, cluster_with_earliset_timepoint, reducedDim){

  require(slingshot); require(SummarizedExperiment)
  sce <- slingshot(data = sce_object, clusterLabels = 'orig.clusters', reducedDim = reducedDim, start.clus = cluster_with_earliset_timepoint)

}

getNewClusters <- function(clusters){
  clusters %>% extractClusterNumber() %>% as.numeric() %>% as.factor() %>% getClusterPhenotypesCml()
}

extractCoarsePhenotype <- function(strs){

  p <- NULL
  i <- 1
  for(str1 in strs){
    p[[i]] <- strsplit(str1, "[ ]")[[1]][2]
    i <- i + 1
  }

  return(p)

}

facets_nice <- theme(strip.background = element_rect(fill="grey96"), strip.text = element_text(colour = 'black'))

reorderClusters <- function(cluster_vec){

  ## Get clusters in order
  clusters <- cluster_vec %>% unique()
  cluster_vec <- factor(as.character(cluster_vec), levels = clusters[order(as.numeric(extractClusterNumber(clusters)))])
  return(cluster_vec)

}


extractClusterNumber <- function(strs){

  p <- NULL
  i <- 1
  for(str1 in strs){
    p[[i]] <- strsplit(str1, "[ ]")[[1]][1]
    i <- i + 1
  }

  return(p)

}


getClusterPhenotypesCml <- function(clusters){

  clusters <- plyr::revalue(clusters, replace   = c("0"  = "0 NK CD56dim",
                                                    "1"  = "1 CD4 CM/naive",
                                                    "2"  = "2 CD4 CM",
                                                    "3"  = "3 CD8 effectory",
                                                    "4"  = "4 CD4 Th1",
                                                    "5"  = "5 CD8 effectory/exhausted",
                                                    "6"  = "6 B-cells naive/resting",
                                                    "7"  = "7 T cells activated",
                                                    "8"  = "8 Monocytes CD16-",
                                                    "9"  = "9 CD4 CM/naive",
                                                    "10" = "10 CD4 treg",

                                                    "11" = "11 CD8 CM/naive",
                                                    "12" = "12 NK adaptive",
                                                    "13" = "13 low quality",
                                                    "14" = "14 CD8 EM",
                                                    "15" = "15 innate lymphocytes",
                                                    "16" = "16 CD4",
                                                    "17" = "17 NK CD56bright",
                                                    "18" = "18 CD4 acute activated",
                                                    "19" = "19 Monocytes CD16+",
                                                    "20" = "20 pDC",
                                                    "21" = "21 CD8 low quality",
                                                    "22" = "22 CD4 low quality",
                                                    "23" = "23 Monocytes low quality"))

  return(clusters)

}




getClusterPhenotypesPublic <- function(clusters){

  plyr::revalue(clusters, replace = c("0"	= "0 CD8 effectory/exhausted",
                                     "1"	= "1 CD4 naive/CM",
                                     "2"	= "2 CD4 CM",
                                     "3"	= "3 NK CD56dim",
                                     "4"	= "4 Monocytes CD16-",

                                     "5"	= "5 CLL cells",
                                     "6"	= "6 B-cells immature",
                                     "7"	= "7 T cells activated",
                                     "8"	= "8 innate lymphocytes \n/ low quality",
                                     "9"	= "9 CD8 CM/naive",

                                     "10" =	"10 B-cells memory",
                                     "11" =	"11 CD4 treg",
                                     "12" =	"12 Monocytes CD16+",
                                     "13" =	"13 CD8 EM",
                                     "14" =	"14 ILC-like",

                                     "15" =	"15 NK adaptive",
                                     "16" =	"16 cDC",
                                     "17" =	"17 CD8+",
                                     "18" =	"18 NK CD56bright",
                                     "19" =	"19 pDC",

                                     "20" =	"20 CD4 Acute Act",
                                     "21" =	"21 NK cycling",
                                     "22" =	"22 Monocytes CD16-",
                                     "23" =	"23 low quality",
                                     "24" =	"24 B-cell plasma"))
}


getClusterPhenotypesTKI <- function(clusters){

  plyr::revalue(clusters, replace = c("0"	 = "0 CD8 effectory/exhausted",
                                      "1"	 = "1 CD4 naive/CM",
                                      "2"	 = "2 NK CD56dim",
                                      "3"	 = "3 CD4 CM",
                                      "4"	 = "4 NK CD56dim",
                                      "5"	 = "5 Monocytes CD16-",
                                      "6"	 = "6 B-cells immature",
                                      "7"	 = "7 innate lymphocytes / low quality",
                                      "8"	 = "8 T cells activated",
                                      "9"	 = "9 low quality",
                                      "10" = "10 CD4 treg",
                                      "11" = "11 CD8 CM/naive",
                                      "12" = "12 B-cells memory",
                                      "13" = "13 CD8 EM",
                                      "14" = "14 ILC-like",
                                      "15" = "15 NK CD56bright",
                                      "16" = "16 Monocytes CD16+",
                                      "17" = "17 Monocytes CD16-",
                                      "18" = "18 pDC",
                                      "19" = "19 NK cycling"))
}

getClusterPhenotypesHealthyTKI <- function(clusters){

  plyr::revalue(clusters, replace = c("0"	 = "0 CD4 naive/CM",
                                      "1"	 = "1 CD4 CM",
                                      "2"	 = "2 CD8 effectory",
                                      "3"	 = "3 NK CD56dim",
                                      "4"	 = "4 NK adaptive",

                                      "5"	 = "5 Monocytes CD16-",
                                      "6"	 = "6 CD8 CM/naive",
                                      "7"	 = "7 T cells activated",
                                      "8"	 = "8 B-cells immature",
                                      "9"	 = "9 CD4 treg",

                                      "10" = "10 B-cells memory",
                                      "11" = "11 CD8 EM",
                                      "12" = "12 innate lymphocytes / low quality",
                                      "13" = "13 B-cells immature",
                                      "14" = "14 NK CD56bright",

                                      "15" = "15 Monocytes CD16+",
                                      "16" = "16 Monocytes",
                                      "17" = "17 Monocytes CD16-",
                                      "18" = "18 pDC",
                                      "19" = "19 NK cycling"))
}


getClusterPhenotypesRelapse <- function(clusters){

  plyr::revalue(clusters, replace = c("0"	 = "0 CD8 effectory/exhausted",
                                      "1"	 = "1 CD4 CM/naive",
                                      "2"	 = "2 NK CD56dim",
                                      "3"	 = "3 CD4 Th1",
                                      "4"	 = "4 B-cells naive/resting",
                                      "5"	 = "5 CD8+ effector",
                                      "6"	 = "6 Monocytes CD16-",
                                      "7"	 = "7 innate lymphocytes / low quality",
                                      "8"	 = "8 CD8+ effector",
                                      "9"	 = "9 T cells activated",
                                      "10" = "10 CD4 treg",
                                      "11" = "11 B-cells naive/resting",
                                      "12" = "12 CD8 CM/naive",
                                      "13" = "13 NK adaptive",
                                      "14" = "14 CD8 EM",
                                      "15" = "15 Monocytes CD16-",
                                      "16" = "16 Monocytes CD16+",
                                      "17" = "17 B-cells naive/resting",
                                      "18" = "18 pDC",
                                      "19" = "19 unspesific"))
}




getClusterColors <- function(){

  nCD4_clusters   <- 10
  nCD8_clusters   <- 7
  nNK_clusters    <- 3
  nOther_clusters <- 11


  cd4_cols   <- colorRampPalette(brewer.pal(nCD4_clusters, "Blues"))(nCD4_clusters)
  cd8_cols   <- colorRampPalette(brewer.pal(nCD8_clusters, "Reds"))(nCD8_clusters)
  nk_cols    <- colorRampPalette(brewer.pal(nNK_clusters, "Greys"))(nNK_clusters)
  other_cols <- colorRampPalette(brewer.pal(nOther_clusters, "Greens"))(nOther_clusters)

  col_vec <- c("cd4" = cd4_cols, "cd8" = cd8_cols, "nk" = nk_cols, "other" = other_cols)
  n_vec   <- c(rep("cd4", nCD4_clusters), rep("cd8", nCD8_clusters),rep("nk", nNK_clusters),  rep("other", nOther_clusters))

  col_vec_fin <- c(nk_cols[1:2],
    cd4_cols[1:2], cd8_cols[1], cd4_cols[3:6], cd8_cols[2:3], cd4_cols[7:9],
    other_cols[1:4], cd8_cols[4], nk_cols[3], other_cols[5:6], cd8_cols[5],
    other_cols[7], cd8_cols[6:7], cd4_cols[10], other_cols[8:11])


  return(col_vec_fin)

}





getClusterCoarsePhenotype <- function(clusters){

  clusters <- plyr::revalue(clusters, replace = c(

    "1" = "nk",
    "2" = "nk",
    "3" = "cd4",
    "4" = "cd4",
    "5" = "cd8",
    "6/emra" = "cd4",
    "7" = "cd4",
    "8" = "cd4",
    "9" = "cd4",
    "10" = "cd8",

    "11" = "cd8",
    "12" = "cd4",
    "13" = "cd4",
    "14" = "cd4",
    "15" = "other",
    "16?" = "other",
    "17" = "other",
    "18" = "other",
    "19/exh" = "cd8",
    "20" = "nk",
    "21_junk" = "other",
    "22_mait" = "other",
    "23_cd8_effectory" = "cd8",
    "24_junk" = "other",
    "25_cd8_eff/exh" = "cd8",
    "26_cd8_eff/exh" ="cd8",
    "27_cd4_acute" = "cd4",
    "28_pDC_1" = "other",
    "29_monocytes" = "other",
    "30_pDC_2" = "other",
    "31_junk" = "other"))

  return(clusters)

}




## For GSEA
forGSEA <- function(de_df){

  ## Create .rnk files for GSEA
  rnk  <- de_df %>% arrange(desc(avg_logFC)) %>% select(gene, avg_logFC)
  colnames(rnk) <- c('Name', 'metric')
  return(rnk)

}



extractName = function(str1){
  sub("\\_.*", "", str1)
}

extractFileName = function(str1){
  sub(".*\\/", "", str1)
}



getLatentUMAP <- function(seurat_object){

  umap_df           <- seurat_object[["latent"]]@cell.embeddings %>% uwot::umap()
  colnames(umap_df) <- c("latent_umap1", "latent_umap2")
  rownames(umap_df) <- colnames(seurat_object)
  umap_df           <- CreateDimReducObject(key = "latent_umap", embeddings = as.matrix(x = umap_df))
  seurat_object[['latent_umap']] <- umap_df

  return(seurat_object)

}

fixSeurat <- function(seurat_object){

  ## Fix meta data if it brokes

  meta.data           <- seurat_object@meta.data
  count.data          <- seurat_object@assays$RNA@counts
  scale.data          <- seurat_object@assays$RNA@scale.data
  # hvg                 <- VariableFeatures(seurat_object)

  # pca_dimred          <- seurat_object[["pca"]]
  # umap_dimred         <- seurat_object[["umap"]]
  latent_dimred       <- seurat_object[["latent"]]
  latent_umap_dimred  <- seurat_object[["latent_umap"]]

  rownames(meta.data) <- meta.data$barcode

  old_idents <- Idents(seurat_object)
  new_seurat <- CreateSeuratObject(counts = count.data)

  new_seurat@meta.data             <- meta.data
  new_seurat@assays$RNA@counts     <- count.data
  new_seurat@assays$RNA@scale.data <- scale.data
  # VariableFeatures(seurat_object)  <- hvg

  # new_seurat[["pca"]]              <- pca_dimred
  # new_seurat[["umap"]]             <- umap_dimred
  new_seurat[["latent"]]           <- latent_dimred
  new_seurat[["latent_umap"]]      <- latent_umap_dimred
  Idents(new_seurat) <- old_idents
  return(new_seurat)

}


getLatentClustering <- function(seurat_object){

  ## Clustering
  res        <- c(seq(0.1, 1, 0.1), seq(1.2, 2, 0.2), 2.5, 3)
  seurat_object <- FindNeighbors(seurat_object, reduction = "latent", dims = c(1:ncol(seurat_object@reductions$latent@cell.embeddings)))
  seurat_object <- FindClusters(object = seurat_object, resolution = res, verbose = F)
  return(seurat_object)

}
