
## Make TCRab file to be used in TCRGP predictions
cml_tcrgp <- data.frame(barcode = cml_seurat$barcode, v = cml_seurat$v_trb, cdr3aa = cml_seurat$trb_cdr3s_aa)
write.table(cml_tcrgp, "results/tcrgp/cml_input.txt", sep = "\t", quote = F, row.names = F)

## Load the original file
message("Read in the files...")
df            <- fread(paste0("results/tcrgp/cml_input.txt"))

## Load the prediction files
pred_files    <- list.files(paste0("results/tcrgp/raw/"), full.names = T, pattern = "csv")
pred_df       <- lapply(pred_files, function(x){fread(x) %>% dplyr::select(prediction)})
pred_df       <- do.call(cbind.data.frame, pred_df)
# tcrb          <- fread(pred_files[1]) %>% pull("CDR3B")

## Add names
pred_epitopes <- extractFileName(pred_files)
# pred_epitopes <- substr(pred_epitopes, 5, nchar(pred_epitopes))
pred_epitopes <- gsub(".csv", "", pred_epitopes)

colnames(pred_df)    <- pred_epitopes
# pred_df$tcrb_padded  <- tcrb
# pred_df$cdr3aa       <- gsub("\\-", "", tcrb)
pred_df$cdr3aa       <- cml_tcrgp$cdr3aa

## Summarise the predictions
message("Getting predictions...")
gotPreds <- getPreds(tcrgp_df = pred_df, fdr = 0.05)
gotPreds <- gotPreds %>% dplyr::select(cdr3aa:clonotype_name, GILGFVFTL_cdr3b:YVLDHLIVV_cdr3b, GLCTLVAML_cdr3b:RPRGEVRFL_cdr3b)
message("Summarising predictions...")
predSummary <- summarisePreds(gotPreds)


## Cbind the two df:s
cml_tcrgp <- cbind(df, predSummary)
write.table(cml_tcrgp, "results/tcrgp/tcrgp_predictions.txt", sep = "\t", quote = F, row.names = F)

## Add TCRGP predictions into Seurat-object
rownames(cml_tcrgp) <- cml_tcrgp$barcode
# cml_seurat@meta.data <- cml_tcrgp

cml_meta  <- cml_seurat@meta.data
cml_tcrgp <- cml_tcrgp[cml_tcrgp$barcode %in% cml_meta$barcode, ]

table(cml_tcrgp$barcode == cml_seurat$barcode)
cml_tcrgp <- cml_tcrgp[,-3]
cml_tcrgp <- cml_tcrgp %>% dplyr::select(CMV_p65_IPS:target)
rownames(cml_tcrgp) <- colnames(cml_seurat)

cml_seurat <- AddMetaData(cml_seurat, metadata = cml_tcrgp)
saveRDS(cml_seurat, "results/cml_fin_19feb_jh.rds")
