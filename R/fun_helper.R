
extractName = function(str1){
  sub("\\_.*", "", str1)
}

extractFileName = function(str1){
  sub(".*\\/", "", str1)
}


filter_contig_files <- function(contig_file){

  require(dplyr)
  contig_file <- contig_file %>%
    filter(high_confidence == "True") %>%
    # filter(chain %in% c("TRA", "TRB", "TRD", "TRG", "multi")) %>%
    filter(productive != "False") %>%
    filter(cdr3 != "None")

  return(contig_file)

}

get_vdj_for_file <- function(file_temp){

  # @ params: file_temp = contig file, which is hopefully filtered
  # @ output: file with consensus V, D and J annotations

  get_vdj <- function(file_temp, clonotype_id, chain){

    file_temp       <- filter(file_temp, chain == chain_temp)
    clonotype_index <- which(as.character(file_temp$clonotype_id) %in% clonotype_id)

    conv_temp       <- file_temp[clonotype_index, ]

    ## Get consensus cdr3aa, v, d and j genes and their respective frequencies
    # conv_temp$cdr3s_nt     <- droplevels(conv_temp$cdr3s_nt)
    # conv_temp$cdr3s_aa     <- droplevels(conv_temp$cdr3s_aa)
    # conv_temp$clonotype_id <- droplevels(conv_temp$clonotype_id)

    cdr3s_nt_all    <- names(sort(table(conv_temp$cdr3s_nt), decreasing = T))
    cdr3s_nt        <- cdr3s_nt_all[1]
    cdr3s_nt_amount <- as.numeric(sort(table(conv_temp$cdr3s_nt), decreasing = T)[1])
    cdr3s_nt_freq   <- cdr3s_nt_amount / sum(table(conv_temp$cdr3s_nt))

    cdr3s_aa_all    <- names(sort(table(conv_temp$cdr3s_aa), decreasing = T))
    cdr3s_aa        <- cdr3s_aa_all[1]
    cdr3s_aa_amouaa <- as.numeric(sort(table(conv_temp$cdr3s_aa), decreasing = T)[1])
    cdr3s_aa_freq   <- cdr3s_aa_amouaa / sum(table(conv_temp$cdr3s_aa))

    v_all           <- names(sort(table(conv_temp$v_gene), decreasing = T))
    v               <- v_all[1]
    v_amount        <- as.numeric(sort(table(conv_temp$v), decreasing = T)[1])
    v_freq          <- v_amount / sum(table(conv_temp$v))

    d               <- names(sort(table(conv_temp$d_gene), decreasing = T)[1])
    d_amount        <- as.numeric(sort(table(conv_temp$d), decreasing = T)[1])
    d_freq          <- d_amount / sum(table(conv_temp$d))

    j_all           <- names(sort(table(conv_temp$j_gene), decreasing = T))
    j               <- j_all[1]
    j_amount        <- as.numeric(sort(table(conv_temp$j), decreasing = T)[1])
    j_freq          <- j_amount / sum(table(conv_temp$j))


    # Combine the results
    tot <- data.frame(clonotype_id, cdr3s_nt, cdr3s_aa, chain, v, d, j, cdr3s_nt_freq, cdr3s_aa_freq, v_freq, d_freq, j_freq)
    return(tot)

  }

  ## Analyse the data
  # file_temp$chain            <- droplevels(file_temp$chain)
  # file_temp$clonotype_id     <- droplevels(file_temp$clonotype_id)

  tcr_list <- list()
  i <- 1

  for(chain_temp in unique(file_temp$chain)){

    print(paste("Analysing chain", chain_temp, "..."))

    ## Analyse only spesific chain at once
    file_temp2              <- filter(file_temp, chain == chain_temp)
    # file_temp2$clonotype_id <- droplevels(file_temp2$clonotype_id)

    ## Go through every clonotype
    vdj <- lapply(unique(file_temp2$clonotype_id), get_vdj, chain = chain_temp, file_temp = file_temp2)
    vdj <- vdj %>% rbindlist(fill = T)
    tcr_list[[i]] <- vdj
    i <- i + 1

  }

  ## Combine TRA & TRB or TRD & TRG; add NAs if only TRA or TRB is found

  tra_ind = NA
  trb_ind = NA
  trg_ind = NA
  trd_ind = NA
  tot_tcrab = NA
  tot_tcrgd = NA

  for(i in 1:length(tcr_list)){

    if(tcr_list[[i]]$chain == "TRA"){tra_ind = i}
    if(tcr_list[[i]]$chain == "TRB"){trb_ind = i}
    if(tcr_list[[i]]$chain == "TRG"){trd_ind = i}
    if(tcr_list[[i]]$chain == "TRD"){trg_ind = i}

  }

  if(!is.na(tra_ind) & !is.na(trb_ind)){
    tot_tcrab <- merge(tcr_list[[tra_ind]], tcr_list[[trb_ind]], by = "clonotype_id", all = T)
  }

  if(!is.na(trg_ind) & !is.na(trd_ind)){
    tot_tcrgd <- merge(tcr_list[[trg_ind]], tcr_list[[trd_ind]], by = "clonotype_id", all = T)
  }

  tot_tcr <- tot_tcrab

  return(tot_tcr)

}



splitCDR3nt <- function(cdr3s_nt){

  require(seqinr)
  tra_cdr3s_nt = ""
  trb_cdr3s_nt = ""

  # print(as.character(cdr3s_nt))
  temp_nt <- s2c(as.character(cdr3s_nt))

  if(sum(temp_nt == ";") > 0){

    tra_cdr3s_nt <- c2s(temp_nt[6:which(temp_nt == ";") - 1])
    trb_cdr3s_nt <- c2s(temp_nt[c(which(temp_nt == ";")+5):length(temp_nt)])

  }

  else{

    name = c2s(temp_nt[1:3])

    if(name == "TRA"){
      tra_cdr3s_nt = c2s(temp_nt[6:length(temp_nt)])

    }

    if(name == "TRB"){
      trb_cdr3s_nt = c2s(temp_nt[5:length(temp_nt)])

    }

  }

  total_cdr3 <- data.frame(cdr3s_nt, tra_cdr3s_nt, trb_cdr3s_nt)
  return(total_cdr3)

}

splitCDR3aa <- function(cdr3s_aa){


  # cdr3s_aa = as.character(s1a1_tcr_clonotype$cdr3s_aa[996])

  require(seqinr)
  tra_cdr3s_aa = ""
  trb_cdr3s_aa = ""

  # priaa(as.character(cdr3s_aa))
  temp_aa <- s2c(as.character(cdr3s_aa))

  if(sum(temp_aa == ";") == 1){

    tra_cdr3s_aa <- c2s(temp_aa[6:which(temp_aa == ";") - 1])
    trb_cdr3s_aa <- c2s(temp_aa[c(which(temp_aa == ";")+5):length(temp_aa)])

  }

  if(sum(temp_aa == ";") == 0){

    name = c2s(temp_aa[1:3])

    if(name == "TRA"){
      tra_cdr3s_aa = c2s(temp_aa[6:length(temp_aa)])
    }

    if(name == "TRB"){
      trb_cdr3s_aa = c2s(temp_aa[5:length(temp_aa)])
    }

  }

  total_cdr3 <- data.frame(cdr3s_aa, tra_cdr3s_aa, trb_cdr3s_aa)
  return(total_cdr3)

}





generateNewClonotypes <- function(old_clonotypes_df){

  temp <- data.frame(table(old_clonotypes_df$cdr3s_nt))
  temp <- temp[order(temp$Freq, decreasing = T), ]
  temp$clonotype_new <- paste0("clonotype", 1:nrow(temp))

  colnames(temp) <- c("cdr3s_nt", "freq", "new_clonotypes_id")
  return(temp)

}

newClonotype_df <- function(old_clonotypes_df){

  # basically just merge but faster (and neater)

  # @ params: old clonotypes df
  # @ output: df with new clonotypes, merged by the same cdr3s_net

  old_clonotypes_df$new_clonotypes_id <- "NA"
  new_clonotypes_df <- generateNewClonotypes(old_clonotypes_df)

  for(cdr3 in unique(old_clonotypes_df$cdr3s_nt)){

    old_clonotypes_df[old_clonotypes_df$cdr3s_nt %in% cdr3, "new_clonotypes_id"] <- new_clonotypes_df[new_clonotypes_df$cdr3s_nt %in% cdr3, "new_clonotypes_id"]

  }

  return(old_clonotypes_df)

}

newClonotype_trb_df <- function(old_clonotypes_df){

  # basically just merge but faster (and neater)

  # @ params: old clonotypes df
  # @ output: df with new clonotypes, merged by the same cdr3s_net



  old_clonotypes_df$new_clonotypes_id <- "NA"
  new_clonotypes_df <- generateNewClonotypes(old_clonotypes_df)

  for(cdr3 in levels(old_clonotypes_df$cdr3_trb)){

    old_clonotypes_df[old_clonotypes_df$cdr3_trb %in% cdr3, "new_clonotypes_id"] <- new_clonotypes_df[new_clonotypes_df$cdr3_trb %in% cdr3, "new_clonotypes_id"]

  }

  return(old_clonotypes_df)

}


## Add conservative cysteine into every cdr3 determined
add_C <- function(cdr3){

  cdr3 <- as.character(cdr3)
  if(nchar(cdr3) > 0){

    if(substr(cdr3, 1, 1) != "C"){
      cdr3 <- paste0("C", cdr3)
    }
  }

  return(cdr3)

}

preprocess_10X_TCR <- function(contig_file, clonotype_file, prefix){

  # contig_file = s1a1_contig
  # clonotype_file = s1a1_clonotype

  # 1) Filter
  contig_file                     <- filter_contig_files(contig_file)
  contig_file                     <- contig_file[,c("barcode", "raw_clonotype_id", "chain", "v_gene", "d_gene", "j_gene", "c_gene")]

  master_tcr_raw                  <- merge(clonotype_file, contig_file, by.x = "clonotype_id", by.y = "raw_clonotype_id")
  master_tcr_raw                  <- master_tcr_raw[!duplicated(master_tcr_raw), ]


  # 2) Get the consensus TCRab and TCRgd files
  master_tcr_raw_detailed             <- get_vdj_for_file(file_temp = master_tcr_raw)
  colnames(master_tcr_raw_detailed)   <- gsub("\\.x", "_tra", colnames(master_tcr_raw_detailed))
  colnames(master_tcr_raw_detailed)   <- gsub("\\.y", "_trb", colnames(master_tcr_raw_detailed))
  master_tcr_raw_detailed$tcr_type    <- ifelse(master_tcr_raw_detailed$chain_tra %in% c("TRA", "TRB"), "TCRab", "TCRgd")

  # 3) Combine the data
  master_tcr_barcoded             <- merge(master_tcr_raw, master_tcr_raw_detailed, by.x = "clonotype_id", by.y = "clonotype_id")
  master_tcr_barcoded             <- dplyr::select(master_tcr_barcoded, -chain, -v_gene, -d_gene, -j_gene, -c_gene)

  # 4) Split the CDR3 information into pieces
  nt <- lapply(master_tcr_barcoded$cdr3s_nt, FUN = splitCDR3nt); nt <- do.call(rbind, nt)
  aa <- lapply(master_tcr_barcoded$cdr3s_aa, FUN = splitCDR3aa); aa <- do.call(rbind, aa)

  master_tcr_barcoded <- data.frame(nt, aa, dplyr::select(master_tcr_barcoded, -cdr3s_nt_tra, -cdr3s_nt_trb, -cdr3s_aa_tra, -cdr3s_aa_trb))

  # 5) Create file just by clonotypes
  master_tcr_clonotype            <- master_tcr_barcoded %>% dplyr::select(-clonotype_id, -barcode)
  master_tcr_clonotype            <- master_tcr_clonotype[!duplicated(master_tcr_clonotype), ]
  master_tcr_clonotype$proportion <- master_tcr_clonotype$frequency / sum(master_tcr_clonotype$frequency)

  # 6) Remove duplicate entries by barcode
  master_tcr_barcoded <- master_tcr_barcoded[!duplicated(master_tcr_barcoded$barcode), ]

  ## Add conservative C if it is missing
  master_tcr_barcoded$tra_cdr3s_aa <- lapply(master_tcr_barcoded$tra_cdr3s_aa, add_C)
  master_tcr_barcoded$trb_cdr3s_aa <- lapply(master_tcr_barcoded$trb_cdr3s_aa, add_C)

  master_tcr_clonotype$tra_cdr3s_aa <- lapply(master_tcr_clonotype$tra_cdr3s_aa, add_C)
  master_tcr_clonotype$trb_cdr3s_aa <- lapply(master_tcr_clonotype$trb_cdr3s_aa, add_C)


  ### ==== Write down ====
  data.table::fwrite(master_tcr_barcoded,  paste0(prefix, "_barcoded.txt"),   sep = "\t", row.names = F)
  data.table::fwrite(master_tcr_clonotype, paste0(prefix, "_clonotyped.txt"), sep = "\t", row.names = F)

}





getSlingNk <- function(project_temp){

  cells.to.keep       <- nk_seurat@meta.data %>% filter(project == project_temp) %>% pull(barcode)
  project_seurat      <- subset(nk_seurat, cells = cells.to.keep) %>% as.SingleCellExperiment()
  proje_nk_sling      <- runSlingshot(project_seurat, cluster_with_earliset_timepoint = "5 CD56 bright", reducedDim = "LATENT_UMAP")
  return(proje_nk_sling)

}


mergeTCRtoSeurat <- function(seurat_object, tcr_df){

    ## Merge with TCRab data with seurat-object metadata
    metadata_temp <- merge(seurat_object@meta.data, dplyr::select(tcr_df, -barcode), all.x = T, by.x = "barcode", by.y = "barcode_uniq")
    metadata_temp <- metadata_temp[match(colnames(seurat_object), metadata_temp$barcode), ]

    ## Add some meta data;

    ## 1) Major: over 10 cells in clonotype
    major_clonotypes               <- unique(subset(metadata_temp, metadata_temp$frequency > 10)$new_clonotypes_id)
    metadata_temp$major_clonotypes <- metadata_temp$new_clonotypes_id
    metadata_temp$major_clonotypes[!metadata_temp$new_clonotypes_id %in% major_clonotypes] <- "minor"

    ## 2) Expanded: over 2 cells in clonotype
    expanded_clonotypes <- unique(subset(metadata_temp, metadata_temp$frequency > 2)$new_clonotypes_id)
    metadata_temp$expanded_clonotypes <- metadata_temp$new_clonotypes_id
    metadata_temp$expanded_clonotypes[!metadata_temp$new_clonotypes_id %in% expanded_clonotypes] <- "unexpanded"

    ## Add metadata into Seurat object; make sure that the colnames match
    rownames(metadata_temp) <- metadata_temp$barcode
    colnames(seurat_object) == rownames(metadata_temp)
    seurat_object@meta.data <- metadata_temp
    return(seurat_object)

}



require(clusterProfiler)
require(org.Hs.eg.db)

if(me == "hru"){

  hallmark   <- read.gmt("/Users/hru/Dropbox/applications/GSEA/h.all.v6.2.symbols.gmt")
  tf         <- read.gmt("/Users/hru/Dropbox/applications/GSEA/c3.tft.v7.0.symbols.gmt")
  go         <- read.gmt("/Users/hru/Dropbox/applications/GSEA/c5.bp.v7.0.symbols.gmt")
  immunology <- read.gmt("/Users/hru/Dropbox/applications/GSEA/c7.all.v7.0.symbols.gmt")


}


if(me == "janihuuh"){

  hallmark   <- read.gmt("/Users/janihuuh/Dropbox/applications/GSEA/h.all.v6.2.symbols.gmt")
  tf         <- read.gmt("/Users/janihuuh/Dropbox/applications/GSEA/c3.tft.v7.0.symbols.gmt")
  go         <- read.gmt("/Users/janihuuh/Dropbox/applications/GSEA/c5.bp.v7.0.symbols.gmt")
  immunology <- read.gmt("/Users/janihuuh/Dropbox/applications/GSEA/c7.all.v7.0.symbols.gmt")

}





plotHypergeometric <- function(genes_df, universe_df, term_df){

  require(clusterProfiler)

  if(nrow(genes_df) == 0) return(NULL)

  # Enrichment in hypergeometric test
  # in: de_df with BM-annotations, direction (Up/Down) and universe to count the enrichment
  # from, e.g. all genes available (~60 000) or all genes expressed in the data (usually ~20 000)

  # out: df with enrichment results

  enrich <- enricher(genes_df$gene, universe = universe_df, TERM2GENE = term_df)

  if(table(enrich@result$p.adjust < 0.05) %>% length() > 1){
    heatplot(enrich)
  }

  else(NULL)

}

getHypergeometric <- function(genes_df, universe_df, term_df){

  require(clusterProfiler)

  if(nrow(genes_df) == 0) return(NULL)

  # Enrichment in hypergeometric test
  # in: de_df with BM-annotations, direction (Up/Down) and universe to count the enrichment
  # from, e.g. all genes available (~60 000) or all genes expressed in the data (usually ~20 000)

  # out: df with enrichment results

  enrich <- enricher(genes_df$gene, universe = universe_df, TERM2GENE = term_df)
  enrich <- do.call(rbind, enrich@result) %>% t %>% as.data.frame()
  enrich[,c(5:7, 9)] <- sapply(enrich[,c(5:7, 9)], function(x) {as.numeric(as.character(x))})
  return(enrich)

}






getSingler <- function(seurat_object, cluster = NULL, method = NULL, sample = NULL){

  message(unique(seurat_object$orig.ident))
  if(!is.null(sample)){

    set.seed(123)
    seurat_object <- subset(seurat_object, cells = colnames(seurat_object)[sample(1:ncol(seurat_object), sample)])

  }

  sce       <- as.SingleCellExperiment(seurat_object)

  ## Predictions
  if(is.null(method)){
    pred.hca <- SingleR::SingleR(test = sce, ref = hpca.se, assay.type.test = 1,   labels = hpca.se$label.fine)
    pred.blu <- SingleR::SingleR(test = sce, ref = blueprint, assay.type.test = 1, labels = blueprint$label.fine)

    if(is.null(sample)){
      seurat_object$singler_hpca_pred      <- pred.hca$first.labels
      seurat_object$singler_blueprint_pred <- pred.blu$first.labels
      return(seurat_object)
    }

    else{
      df <- data.frame(barcode = rownames(pred.hca), cluster = seurat_object$cluster, singler_hpca_pred = pred.hca$labels, singler_blueprint_pred = pred.blu$labels)
      return(df)
    }

  }


  if(method == "cluster"){
    if(is.null(cluster)){
      cluster=seurat_object$cluster
    }
    pred.hca <- SingleR::SingleR(test = sce, ref = hpca.se, assay.type.test = 1,   labels = hpca.se$label.fine, method = "cluster", clusters = cluster)
    pred.blu <- SingleR::SingleR(test = sce, ref = blueprint, assay.type.test = 1, labels = blueprint$label.fine, method = "cluster", clusters = cluster)
    df <- data.frame(cluster = rownames(pred.hca), singler_hpca_pred = pred.hca$labels, singler_blueprint_pred = pred.blu$labels)
    return(df)
  }
}



plotClustering <- function(seurat_object){

  res <- c(seq(0.1, 1, 0.1), seq(1.2, 2, 0.2), 2.5, 3)
  clustering_columns <- grep("res", colnames(seurat_object@meta.data), value = T)
  clustering_columns <- clustering_columns[order(substr(clustering_columns, 13, nchar(clustering_columns)) %>% as.numeric())]

  q <- NULL; i <- 1

  for(clustering_column in clustering_columns){
    q[[i]] <- seurat_object@meta.data[,clustering_column] %>% levels %>% length
    i <- i + 1
  }

  data.frame(resolution = res, nClusters = do.call(q, what="c")) %>%
    ggplot(aes((resolution),nClusters), label = nClusters) + geom_point(shape = 21) + theme_bw()

}

putLatentsSeurat <- function(seurat_object, latent){

  latent_umap <- uwot::umap(latent) %>% as.data.frame() %>% dplyr::rename(UMAP1=V1, UMAP2=V2)

  latent      <- as.matrix(latent)
  latent_umap <- as.matrix(latent_umap)

  rownames(latent)      <- colnames(seurat_object)
  rownames(latent_umap) <- colnames(seurat_object)

  latent_dim_red            <- CreateDimReducObject(key = "latent", embeddings = as.matrix(x = latent))
  latent_umap_dim_red       <- CreateDimReducObject(key = "latent_umap", embeddings = as.matrix(x = latent_umap))

  seurat_object[['latent']]      <- latent_dim_red
  seurat_object[['latent_umap']] <- latent_umap_dim_red
  return(seurat_object)
}

getScviInput <- function(seurat_object, folder){

  dir.create(folder, showWarnings = F)
  genes_to_keep <- list()
  i <- 1
  Idents(seurat_object) <- seurat_object$orig.ident

  for(patient in unique(seurat_object$orig.ident)){

    message(patient)
    seurat_temp <- subset(seurat_object, idents = patient)
    counts_temp <- seurat_temp@assays$RNA@counts %>% as.data.frame

    genes_to_keep[[i]] <- rownames(counts_temp)
    i <- i + 1

    counts_temp <- seurat_object@assays$RNA@counts[ ,seurat_object$orig.ident == patient] %>% as.data.frame
    counts_temp <- counts_temp[!rownames(counts_temp) %in% clonality_genes, ]
    data.table::fwrite(counts_temp, paste0(folder, patient, ".csv"), sep = ",", quote = F, row.names = T, col.names = T)

  }
}


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








getQC <- function(seurat_object){

  ###################

  min_mito     <- 0
  max_mito     <- 10

  min_ribo     <- 5
  max_ribo     <- 50

  min_features <- 500
  max_features <- 3000
  # max_features <- 4500

  min_counts   <- 1e3
  max_counts   <- 30e3

  min_pct50    <- 25
  max_pct50    <- 60

  ###################

  seurat_object@meta.data$barcode <- colnames(seurat_object)

  ## In total, we remove with the following conditions:
  qc_df <- seurat_object@meta.data %>% as.data.frame()

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





getLatentClustering <- function(seurat_object){

  ## Clustering
  res        <- c(seq(0.1, 1, 0.1), seq(1.2, 2, 0.2), 2.5, 3)
  seurat_object <- FindNeighbors(seurat_object, reduction = "latent", dims = c(1:ncol(seurat_object@reductions$latent@cell.embeddings)))
  seurat_object <- FindClusters(object = seurat_object, resolution = res, verbose = F)
  return(seurat_object)

}







getReducedNames <- function(strs){

  p <- NULL
  i <- 1
  for(str1 in strs){
    x <- strsplit(str1, "[ ]")[[1]][c(1,3)]
    p[[i]] <- paste(x, collapse = " ")
    i <- i + 1
  }

  return(p)


}

getNewClusters <- function(clusters){

  clusters %>% extractClusterNumber() %>% as.numeric %>% as.factor() %>% extractClusterNumber() %>% getClusterPhenotypes()

}

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

facets_nice <- theme(strip.background = element_rect(fill="grey96"), strip.text = element_text(colour = 'black'))

extractName = function(str1){
  # strsplit(str1, "[_]")[[1]][1]
  sub("\\_.*", "", str1)
}

extractFileName = function(str1){
  # strsplit(str1, "[_]")[[1]][1]
  sub(".*\\/", "", str1)
}

extractSeuratName <- function(str1){

  str1 <- substr(str1, 1, nchar(str1) - 27)
  extractFileName(str1)
}

extractTimepoint <- function(strs){

  strs2 <-NULL
  i <- 1
  for(str1 in strs){
    strs2[[i]] <- strsplit(str1, "[_]")[[1]][2]
    i <- i + 1
  }

  # return(strs2)
  return(factor(strs2, levels = c("311019", "141119", "021219", "301219", "130120", "280120")))

}


plotQcViolin <- function(viz_df, var_to_plot, grouping, min, max){

  ## Plot univariate violin plots with filter thresholds

  # @ params:
  # viz_df = df that contains qc-analysis results and covariates of interest
  # var_to_plot = char, a column name that contains the variable to plot
  # grouping = char, a column name that contains the x-axis grouping
  # min = num, min value for variable
  # max = num, max value for variable

  viz_df_temp <- viz_df %>% dplyr::select(var_to_plot)

  label_df_min <- ifelse(viz_df_temp > min, "above", "below") %>% table
  label_df_max <- ifelse(viz_df_temp < max, "above", "below") %>% table

  ggplot(data = viz_df, aes_string(x = grouping, y = var_to_plot, fill = grouping)) +
    geom_violin(alpha = 0.5) +
    # geom_boxplot(width = 0.1, alpha = 0.8, outlier.shape = NA) +

    geom_hline(yintercept = min, linetype = "dotted") +
    geom_hline(yintercept = max, linetype = "dotted") +

    annotate(geom = "text", x = 2.5, y = min, label = paste("Below the line:\n", label_df_min[2]), fontface = "italic") +
    annotate(geom = "text", x = 2.5, y = max, label = paste("Above the line:\n", label_df_max[2]), fontface = "italic") +

    labs(x = "", title = var_to_plot) + theme(legend.position = "none")

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






getLatentUMAP <- function(seurat_object){

  umap_df           <- seurat_object[["latent"]]@cell.embeddings %>% uwot::umap()
  colnames(umap_df) <- c("latent_umap1", "latent_umap2")
  rownames(umap_df) <- colnames(seurat_object)
  umap_df           <- CreateDimReducObject(key = "latent_umap", embeddings = as.matrix(x = umap_df))
  seurat_object[['latent_umap']] <- umap_df

  return(seurat_object)

}









initCellphonedb <- function(df, seurat_object, name = "", sample_n = 50, folder = "results/cellphonedb/input_files/", seed = 123){

  set.seed(seed)

  ## If one wants to use subsampling; remove clusters with less than sample_n cells
  if(!is.null(sample_n)){

    clusters.to.rm  <- df %>% group_by(cluster) %>% summarise(n = n()) %>% filter(n < sample_n) %>% pull(cluster)
    message(paste("remove", clusters.to.rm))
    clusters.to.use <- df %>% group_by(cluster) %>% summarise(n = n()) %>% filter(n >= sample_n) %>% pull(cluster)
    df              <- df %>% filter(cluster %in% clusters.to.use) %>% group_by(cluster)
    df              <- df %>% sample_n(size = sample_n) %>% ungroup()
  }

  ## Remove empty clusters
  df <- df %>% filter(cluster != "")

  counts <- seurat_object@assays$RNA@counts[,df$barcode]
  counts[Matrix::rowSums(counts) != 0, ] %>% as.data.frame %>% add_rownames(var = "Gene") %>%
    dplyr::select(Gene, everything()) %>%
    fwrite(paste0(folder, name, "_counts.txt"), sep = "\t", quote = F, row.names = F)

  df %>% dplyr::select(barcode, cluster) %>% dplyr::rename(Cell = barcode, cell_type = cluster) %>%
    fwrite(paste0(folder, name, "_meta.txt"), sep = "\t", quote = F, row.names = F)

  return(NULL)

}
