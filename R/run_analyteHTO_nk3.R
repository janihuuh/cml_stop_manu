
cml3_seurat <- cml_seurat_new_filt3

## Select only tumor cells
cml3_seurat_tumor <- subset(cml3_seurat, is_nk != "NK")
cml3_seurat_tumor <- subset(cml3_seurat_tumor, my.demux != "E-NK-only")
cml3_seurat_tumor <- cml3_seurat_tumor %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml3_seurat_tumor), nPCs = 14)

## Select only NK cells
cml3_seurat$is_nk <- ifelse(cml3_seurat$RNA_snn_res.0.1 == 1 & grepl("NK", cml3_seurat$my.demux), "NK", "CML")

cml3_seurat_nk <- subset(cml3_seurat, is_nk == "NK")
cml3_seurat_nk <- cml3_seurat_nk %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml3_seurat_nk), nPCs = 10)

saveRDS(cml3_seurat_tumor, "results/cml3_seurat_tumor.rds")
saveRDS(cml3_seurat_nk, "results/cml3_seurat_nk.rds")
