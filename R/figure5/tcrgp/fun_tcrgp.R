


## Load the thresholds
source("src/R/tcrgp/load_thresholds.R")

## Generate a function that understands the output
getPreds <- function(tcrgp_df, thresholds_df = thresholds, fdr = 0.05){

  ## This function removes all predictions that are under the selected threshold, and leaves anything else behind
  ## summarisePreds() is adviced to use after this

  # @ param
  ## input:
  # tcrgp_df: df which contains the epitopes as column names. Can contain other meta-info as well
  # thresholds_df: df which contains the thresholds for predictions for the epitoeps in tcrgp_df
  # fdr: or fpr; the criterion to use for thresholds

  require(dplyr)

  ## Select the fdr for the look-up table, thresholds
  if(fdr == 0){fdr = "FPRS_0"}
  if(fdr == 0.05){fdr = "FPRS_005"}
  if(fdr == 0.1){fdr = "FPRS_01"}
  if(fdr == 0.2){fdr = "FPRS_02"}

  ## Read in the prediction data
  tcrgp_df <- tcrgp_df %>% mutate(clonotype_name = cdr3aa)
  tcrgp_df[tcrgp_df == "NaN"] <- NA

  thresholds_df$model
  colnames(tcrgp_df)
  ## Select only the epitopes found in this TCRGP prediction
  thresholds_temp   <- thresholds %>% dplyr::slice(which(thresholds_df$model %in% colnames(tcrgp_df)))
  thresholds_to_use <- thresholds_temp %>% pull(fdr)

  ## In the TCRGP-tcrgp_df, get rid of meta data
  tcrgp_df_meta <- as.data.frame(tcrgp_df)[ , which(!colnames(tcrgp_df) %in% thresholds_temp$model)]

  ## Select only the columns with the model and reorder accodingly
  tcrgp_df_temp <- as.data.frame(tcrgp_df)[ , which(colnames(tcrgp_df) %in% thresholds_temp$model)]
  tcrgp_df_temp <- tcrgp_df_temp[ ,match(colnames(tcrgp_df_temp), thresholds_temp$model)]

  if(colnames(tcrgp_df_temp) == thresholds_temp$model){

    ## Main part: filter the cdr3s with low predictions
    antigen_names <- colnames(tcrgp_df_temp)
    tcrgp_df_temp <- pbapply::pbapply(tcrgp_df_temp, 1, function(x){ x[x < thresholds_to_use] <- 0; return(data.frame(matrix(x,nrow=1)))})
    tcrgp_df_temp <- do.call(rbind, tcrgp_df_temp)

    colnames(tcrgp_df_temp) <- antigen_names
    tcrgp_df_temp           <- cbind(tcrgp_df_meta, tcrgp_df_temp)

    return(tcrgp_df_temp)

  }

  else{print("No matching models")}

}

summarisePreds <- function(pred_file, thresholds_df = thresholds){

  ## This function removes all predictions that are under the selected threshold, and leaves anything else behind
  ## summarisePreds() is adviced to use after this

  # @ param
  ## input:
  # pred_file: df from getPreds()
  # thresholds: df which contains the thresholds for predictions for the epitoeps in tcrgp_df


  ## Get the models from the threshold look-up table
  pred_file_meta <- pred_file %>% dplyr::select(which(!colnames(pred_file) %in% thresholds_df$model))
  pred_file      <- pred_file %>% dplyr::select(which(colnames(pred_file) %in% thresholds_df$model))
  pred_file[is.na(pred_file)] <- 0

  ## To find predictions, first define which cells deviate from 0
  getNonzeros <- function(row){

    ## If no prediction is found:
    if(sum(row) == 0){
      pred_epitope    <- ""
      n_pred_epitopes <- 0
      return(data.frame(pred_epitope = pred_epitope, n_pred_epitopes = n_pred_epitopes))
    }

    antigen_names     <- names(row)
    specific_antigens <- which(row > 0)
    pred_epitope      <- paste(antigen_names[specific_antigens], collapse = ", ")
    n_pred_epitopes   <- length(specific_antigens)

    return(data.frame(pred_epitope = pred_epitope, n_pred_epitopes = n_pred_epitopes))

  }

  summ_file      <- do.call(rbind, pbapply::pbapply(pred_file, 1, FUN = getNonzeros))

  message("Summarisation done")

  ## Summarise the targets
  summ_file$target                                 <- "anti-viral"
  summ_file$target[summ_file$n_pred_epitopes > 1]  <- "multi"
  summ_file$target[summ_file$n_pred_epitopes == 0] <- "none"

  df <- cbind(pred_file_meta, pred_file, summ_file)
  return(df)

}



fisherEpitope <- function(seurat_object, epitope){


  ## Calculate Fisher's exact for epitope enrichment to clusters
  ## Clusters are Idents from seurat object

  # @ params
  # seurat_object = seurat_object
  # epitope = char, epitope to predict


  cell_df <- cbind(seurat_object@meta.data, cluster = Idents(seurat_object))
  target_df  <- cell_df %>% filter(pred_epitope %in% epitope)

  target_cells             = nrow(target_df)
  cells                    = nrow(cell_df) - target_cells

  target_cells_in_cluster  = table(target_df$cluster) %>% data.frame()
  cells_in_cluster         = table(cell_df$cluster) %>% data.frame()

  cell_df <- merge(cells_in_cluster, target_cells_in_cluster, by = "Var1", all.x = T)
  cell_df$Freq.y[is.na(cell_df$Freq.y)] <- 0

  colnames(cell_df) <- c("cluster", "cells_in_cluster", "target_cells_in_cluster")
  cell_df$cells        <- cells        - cell_df$cells_in_cluster
  cell_df$target_cells <- target_cells - cell_df$target_cells_in_cluster

  cell_df <- cell_df[,c(1,3,5,2,4)]

  if(length(epitope) > 1){epitope = paste(epitope, collapse = ", ")}

  fisher_df <- do.call(rbind, apply(cell_df[,-1], 1, function(x) broom::tidy(fisher.test(matrix(as.numeric(x), byrow = T, nrow = 2), alternative = "greater"))))
  total <- cbind(cell_df, fisher_df) %>%
    mutate(p.adj = p.adjust(p.value, method = "BH"), epitope = epitope) %>%
    mutate(sigf = ifelse(p.adj < 0.05, "Sigf", "Not_sigf")) %>% dplyr::select(-method) %>%
    arrange(p.adj)

  return(total)

}





