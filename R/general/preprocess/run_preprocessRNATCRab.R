
message("Merge seurat-objects with TCRab-info...")

mergeTCRtoSeurat <- function(seurat_object, tcr_df){

    ## Merge with TCRab data with seurat-object metadata
    metadata_temp <- merge(seurat_object@meta.data, select(tcr_df, -barcode), all.x = T, by.x = "barcode", by.y = "barcode_uniq")
    metadata_temp <- metadata_temp[match(colnames(seurat_object), metadata_temp$barcode), ]

    ## Add some meta data;

    ## 1) Major: over 10 cells in clonotype
    major_clonotypes               <- unique(subset(metadata_temp, metadata_temp$frequency.y > 10)$new_clonotypes_id.y)
    metadata_temp$major_clonotypes <- metadata_temp$new_clonotypes_id.y
    metadata_temp$major_clonotypes[!metadata_temp$new_clonotypes_id.y %in% major_clonotypes] <- "minor"

    ## 2) Expanded: over 2 cells in clonotype
    expanded_clonotypes <- unique(subset(metadata_temp, metadata_temp$frequency.y > 2)$new_clonotypes_id.y)
    metadata_temp$expanded_clonotypes <- metadata_temp$new_clonotypes_id.y
    metadata_temp$expanded_clonotypes[!metadata_temp$new_clonotypes_id.y %in% expanded_clonotypes] <- "unexpanded"

    ## Add metadata into Seurat object; make sure that the colnames match
    rownames(metadata_temp) <- metadata_temp$barcode
    colnames(seurat_object) == rownames(metadata_temp)
    seurat_object@meta.data <- metadata_temp
    return(seurat_object)

}

seurat_object$barcode <- colnames(seurat_object)
seurat_object <- mergeTCRtoSeurat(seurat_object = seurat_object, tcr_df = tot_barcode)
unique(subset(tot_barcode, tot_barcode$frequency > 10)$new_clonotypes_id)
saveRDS(cml_seurat, "results/cml_seurat_fin.rds")

metadata_temp <- merge(seurat_object@meta.data, select(tot_barcode, -barcode), all.x = T, by.x = "barcode", by.y = "barcode_uniq")
metadata_temp <- metadata_temp[match(colnames(seurat_object), metadata_temp$barcode), ]
major_clonotypes               <- unique(subset(metadata_temp, metadata_temp$frequency > 10)$new_clonotypes_id)

metadata <- seurat_object[[]]
table(metadata$patient.y)

# loop and check that everything with .x goes to null and .y gets the .y removed

i <- 1
for (cname in colnames(metadata)){
  if (is.na(strsplit(cname, '[.]')[[1]][2]) == FALSE){
    if (strsplit(cname, '[.]')[[1]][2] == 'y'){colnames(metadata)[i] <- strsplit(cname, '[.]')[[1]][1]}
  }
  i <- i + 1
}

# Removing all the old stuff with 'x'
colnames(metadata)
colnames(metadata)[12]
colnames(metadata)[43]
metadata[12:43] <- NULL

# And adding it to the Seurat object
seurat_object@meta.data <- metadata

metadata2 <- seurat_object[[]]
head(metadata2) # All seems fine!

setwd('/Users/theodjas/Dropbox (Aalto)/cml_stop/latest/data')
saveRDS(seurat_object, 'cml_fin_28jan.RDS')

seurat_object <- readRDS('cml_fin_28jan.RDS')




