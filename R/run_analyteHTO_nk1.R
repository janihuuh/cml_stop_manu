
cml1_seurat <- cml_seurat_new_filt1
cml1_seurat$is_nk <- ifelse(cml1_seurat$RNA_snn_res.0.1 == 2 & grepl("NK", cml1_seurat$my.demux), "NK", "CML")

## Select only tumor cells
cells.to.keep     <- cml1_seurat_tumor@meta.data %>% filter(is_nk != "NK" & !my.demux %in% c("E-NK-only", "NE-NK-only")) %>% pull(barcode)
cml1_seurat_tumor <- subset(cml1_seurat_tumor, cells = cells.to.keep)
cml1_seurat_tumor <- cml1_seurat_tumor %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml1_seurat_tumor), nPCs = 14)

## Select only tumor cells
cells.to.keep  <- cml1_seurat_tumor@meta.data %>% filter(is_nk == "NK") %>% pull(barcode)
cml1_seurat_nk <- subset(cml1_seurat, cells = cells.to.keep)
cml1_seurat_nk <- cml1_seurat_nk %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml1_seurat_nk), nPCs = 12)

saveRDS(cml1_seurat_tumor, "results/functional/cml1_seurat_tumor.rds")
saveRDS(cml1_seurat_nk, "results/functional/cml1_seurat_nk.rds")
