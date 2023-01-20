
dense_umis <- fread("data/functional/citeseq/CML_NK_1_HTO/uncorrected_cells/dense_umis.tsv", gene.column = 1)

cml_nk_1_hto_seurat <- Read10X(data.dir = "data/functional/citeseq/CML_NK_1_HTO/umi_count/", gene.column = 1)
cml_nk_1_hto_seurat <- CreateSeuratObject(cml_nk_1_hto_seurat, project = "CML_NK_1")

cml_nk_2_hto_seurat <- Read10X(data.dir = "data/functional/citeseq/CML_NK_2_HTO/umi_count/", gene.column = 1)
cml_nk_2_hto_seurat <- CreateSeuratObject(cml_nk_2_hto_seurat, project = "CML_NK_2")

cml_nk_3_hto_seurat <- Read10X(data.dir = "data/functional/citeseq/CML_NK_3_HTO/umi_count/", gene.column = 1)
cml_nk_3_hto_seurat <- CreateSeuratObject(cml_nk_3_hto_seurat, project = "CML_NK_3")

colnames(cml_nk_1_hto_seurat) <- paste0("CML_NK_1_", colnames(cml_nk_1_hto_seurat), "-1")
colnames(cml_nk_2_hto_seurat) <- paste0("CML_NK_2_", colnames(cml_nk_2_hto_seurat), "-1")
colnames(cml_nk_3_hto_seurat) <- paste0("CML_NK_3_", colnames(cml_nk_3_hto_seurat), "-1")

hto_names <- c(paste0("CML_NK_1_", colnames(cml_nk_1_hto_seurat), "-1"), 
               paste0("CML_NK_2_", colnames(cml_nk_2_hto_seurat), "-1"), 
               paste0("CML_NK_3_", colnames(cml_nk_3_hto_seurat), "-1"))

joint.bcs <- intersect(colnames(cml_seurat), hto_names)

cml_nk_1_hto_seurat_new <- cml_nk_1_hto_seurat[, c(paste0("CML_NK_1_", colnames(cml_nk_1_hto_seurat), "-1") %in% joint.bcs)]
cml_nk_2_hto_seurat_new <- cml_nk_2_hto_seurat[, c(paste0("CML_NK_2_", colnames(cml_nk_2_hto_seurat), "-1") %in% joint.bcs)]
cml_nk_3_hto_seurat_new <- cml_nk_3_hto_seurat[, c(paste0("CML_NK_3_", colnames(cml_nk_3_hto_seurat), "-1") %in% joint.bcs)]

cml_nk_1_hto_mtx <- cml_nk_1_hto_seurat_new@assays$RNA@counts
colnames(cml_nk_1_hto_mtx) <- c(paste0("CML_NK_1_", colnames(cml_nk_1_hto_mtx), "-1"))
rownames(cml_nk_1_hto_mtx) <- c(hto_id %>% filter(experiment == "CML_NK_1") %>% pull(toname), "unmapped")

cml_nk_2_hto_mtx <- cml_nk_2_hto_seurat_new@assays$RNA@counts
cml_nk_2_hto_mtx <- cml_nk_2_hto_mtx[c(1:10,15), ]
colnames(cml_nk_2_hto_mtx) <- c(paste0("CML_NK_2_", colnames(cml_nk_2_hto_mtx), "-1"))
rownames(cml_nk_2_hto_mtx) <- c(hto_id %>% filter(experiment == "CML_NK_2") %>% pull(toname), "unmapped")

cml_nk_3_hto_mtx <- cml_nk_3_hto_seurat_new@assays$RNA@counts
cml_nk_3_hto_mtx <- cml_nk_3_hto_mtx[c(1:10,15), ]
colnames(cml_nk_3_hto_mtx) <- c(paste0("CML_NK_3_", colnames(cml_nk_3_hto_mtx), "-1"))
rownames(cml_nk_3_hto_mtx) <- c(hto_id %>% filter(experiment == "CML_NK_3") %>% pull(toname), "unmapped")

cml_seurat_new1 <- subset(cml_seurat_new, orig.ident == "CML_NK_1")
cml_seurat_new2 <- subset(cml_seurat_new, orig.ident == "CML_NK_2")
cml_seurat_new3 <- subset(cml_seurat_new, orig.ident == "CML_NK_3")

cml_seurat_new1[["HTO"]] <- CreateAssayObject(counts = cml_nk_1_hto_mtx)
cml_seurat_new2[["HTO"]] <- CreateAssayObject(counts = cml_nk_2_hto_mtx)
cml_seurat_new3[["HTO"]] <- CreateAssayObject(counts = cml_nk_3_hto_mtx)


#### HTO 1
cml_seurat_new1 <- cml_seurat_new1 %>% NormalizeData(assay = "HTO", normalization.method = "CLR") %>% HTODemux(assay = "HTO", positive.quantile = 0.99)

DefaultAssay(cml_seurat_new1) <- "HTO"
cml_seurat_new1 <- cml_seurat_new1 %>% ScaleData(features = rownames(cml_seurat_new1)) %>% RunPCA(features = rownames(cml_seurat_new1), approx = F) %>% RunTSNE(dims = 1:8, perplexity = 100)
cml_seurat_new1 <- cml_seurat_new1 %>% getDemuxEssi
cml_seurat_new_filt1 <- subset(cml_seurat_new1, my.demux2 == "Singlet")
cml_seurat_new_filt1 <- cml_seurat_new_filt1 %>% ScaleData(features = rownames(cml_seurat_new_filt1)) %>% RunPCA(features = rownames(cml_seurat_new_filt1), approx = F) %>% RunTSNE(dims = 1:8, perplexity = 100)

DimPlot(cml_seurat_new_filt1, label = T, repel = T, group.by = "HTO_classification.global")
DimPlot(cml_seurat_new_filt1, label = T, repel = T, group.by = "my.demux")
DimPlot(cml_seurat_new_filt1, label = T, repel = T, group.by = "my.demux2")
HTOHeatmap(cml_seurat_new_filt1, assay = "HTO", ncells = 5000)

cml_seurat_new_filt1 <- subset(cml_seurat_hto, orig.ident == "CML_NK_1")

DefaultAssay(cml_seurat_new_filt1) <- "RNA"
cml_seurat_new_filt1 <- cml_seurat_new_filt1 %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml_seurat_new_filt1))
cml_seurat_new_filt1 <- cml_seurat_new_filt1 %>% getClustering()

a <- DimPlot(cml_seurat_new_filt1, label = F, repel = T, cols = getPalette(14)) + labs(title = "Seurat")
b <- DimPlot(cml_seurat_new_filt1, label = F, repel = T, cols = getPalette(14), group.by = "my.demux") + labs(title = "Essi")

c <- DimPlot(cml_seurat_new_filt1, label = F, repel = T, cols = getPalette(14), group.by = "HTO_classification.global") + labs(title = "Seurat")
d <- DimPlot(cml_seurat_new_filt1, label = F, repel = T, cols = getPalette(14), group.by = "my.demux2") + labs(title = "Seurat")

e <- a + b + c + d
ggsave(plot = e, "results/functional/demux/umap_hto1_rna.png", width = 12, height = 8)

a <- DimPlot(cml_seurat_new_filt1, label = T, repel = T, cols = getPalette(14), group.by = "RNA_snn_res.0.1") 
b <- DimPlot(cml_seurat_new_filt1, label = T, repel = T, cols = getPalette(14), group.by = "singler_blueprint_pred") 
# c <- DimPlot(cml_seurat_new_filt1, label = T, repel = T, cols = getPalette(14), group.by = "singler_hpca_pred") 

cml_seurat_new_filt1$is_nk <- ifelse(cml_seurat_new_filt1$RNA_snn_res.0.1 == 2, "NK", "CML")

Idents(cml_seurat_new_filt1) <- cml_seurat_new_filt1$RNA_snn_res.0.1
cml_seurat_new_markers <- FindAllMarkers(cml_seurat_new_filt1, test.use = "t", max.cells.per.ident = 1e3) %>% filter(avg_log2FC > 0)
cml_seurat_new_markers <- cml_seurat_new_markers %>% filter(avg_log2FC > 0)

FeaturePlot(cml_seurat_new_filt1, features = big_markers)
ggsave("results/functional/demux/feature_hto1.png", width = 12, height = 18)


DefaultAssay(cml_seurat_new_filt1) <- "HTO"
cml_seurat_new_filt1 <- ScaleData(cml_seurat_new_filt1, features = rownames(cml_seurat_new_filt1), verbose = FALSE)
cml_seurat_new_filt1 <- RunPCA(cml_seurat_new_filt1, features = rownames(cml_seurat_new_filt1), approx = FALSE)
cml_seurat_new_filt1 <- RunTSNE(cml_seurat_new_filt1, dims = 1:8, perplexity = 100)

table(cml_seurat_new_filt1$HTO_classification.global)

a <- DimPlot(cml_seurat_new_filt1, label = F, repel = T, cols = getPalette(14)) + labs(title = "Seurat")
b <- DimPlot(cml_seurat_new_filt1, label = F, repel = T, cols = getPalette(14), group.by = "my.demux") + labs(title = "Essi")

c <- DimPlot(cml_seurat_new_filt1, label = F, repel = T, cols = getPalette(14), group.by = "HTO_classification.global") + labs(title = "Seurat")
d <- DimPlot(cml_seurat_new_filt1, label = F, repel = T, cols = getPalette(14), group.by = "my.demux2") + labs(title = "Seurat")

e <- DimPlot(cml_seurat_new_filt1, label = F, repel = T, cols = getPalette(2), group.by = "is_nk") + labs(title = "RNA")

f <- a + b + c + d + e
ggsave(plot = f, "results/functional/demux/umap_hto1_hto.png", width = 12, height = 6)



##### HTO 2

cml_seurat_new2 <- cml_seurat_new2 %>% NormalizeData(assay = "HTO", normalization.method = "CLR") %>% HTODemux(assay = "HTO", positive.quantile = 0.99)

DefaultAssay(cml_seurat_new2) <- "HTO"
cml_seurat_new2 <- cml_seurat_new2 %>% ScaleData(features = rownames(cml_seurat_new2)) %>% RunPCA(features = rownames(cml_seurat_new2), approx = F) %>% RunTSNE(dims = 1:8, perplexity = 100)
cml_seurat_new2 <- cml_seurat_new2 %>% getDemuxEssi
cml_seurat_new_filt2 <- subset(cml_seurat_new2, my.demux2 == "Singlet")
cml_seurat_new_filt2 <- cml_seurat_new_filt2 %>% ScaleData(features = rownames(cml_seurat_new_filt2)) %>% RunPCA(features = rownames(cml_seurat_new_filt2), approx = F) %>% RunTSNE(dims = 1:8, perplexity = 100)

DimPlot(cml_seurat_new_filt2, label = T, repel = T, group.by = "HTO_classification.global")
DimPlot(cml_seurat_new_filt2, label = T, repel = T, group.by = "my.demux")
DimPlot(cml_seurat_new_filt2, label = T, repel = T, group.by = "my.demux2")
HTOHeatmap(cml_seurat_new_filt2, assay = "HTO", ncells = 5000)


cml_seurat_new_filt2$hto_classification_single <- Idents(cml_seurat_new_filt2)


DefaultAssay(cml_seurat_new_filt2) <- "RNA"
cml_seurat_new_filt2 <- cml_seurat_new_filt2 %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml_seurat_new_filt2))
cml_seurat_new_filt2 <- cml_seurat_new_filt2 %>% getClustering()
Idents(cml_seurat_new_filt2) <- cml_seurat_new_filt2$RNA_snn_res.0.1

a <- DimPlot(cml_seurat_new_filt2, label = F, repel = T, cols = getPalette(14)) + labs(title = "RNA_snn_res.0.1")
b <- DimPlot(cml_seurat_new_filt2, label = F, repel = T, cols = getPalette(14), group.by = "singler_blueprint_pred") + labs(title = "SingleR")
c <- DimPlot(cml_seurat_new_filt2, label = F, repel = T, cols = getPalette(14), group.by = "HTO_classification.global") + labs(title = "Seurat")
d <- DimPlot(cml_seurat_new_filt2, label = F, repel = T, cols = getPalette(14), group.by = "my.demux") + labs(title = "Essi")

e <- a + b + c + d
ggsave(plot = e, "results/functional/demux/umap_hto2_rna.png", width = 12, height = 8)






DefaultAssay(cml_seurat_new_filt2) <- "HTO"
cml_seurat_new_filt2 <- ScaleData(cml_seurat_new_filt2, features = rownames(cml_seurat_new_filt2), verbose = FALSE)
cml_seurat_new_filt2 <- RunPCA(cml_seurat_new_filt2, features = rownames(cml_seurat_new_filt2), approx = FALSE)
cml_seurat_new_filt2 <- RunTSNE(cml_seurat_new_filt2, dims = 1:8, perplexity = 100)




###### HTO 3
cml_seurat_new3 <- cml_seurat_new3 %>% NormalizeData(assay = "HTO", normalization.method = "CLR") %>% HTODemux(assay = "HTO", positive.quantile = 0.99)

DefaultAssay(cml_seurat_new3) <- "HTO"
cml_seurat_new3 <- cml_seurat_new3 %>% ScaleData(features = rownames(cml_seurat_new3)) %>% RunPCA(features = rownames(cml_seurat_new3), approx = F) %>% RunTSNE(dims = 1:8, perplexity = 100)
cml_seurat_new3 <- cml_seurat_new3 %>% getDemuxEssi
cml_seurat_new_filt3 <- subset(cml_seurat_new3, my.demux2 == "Singlet")
cml_seurat_new_filt3 <- cml_seurat_new_filt3 %>% ScaleData(features = rownames(cml_seurat_new_filt3)) %>% RunPCA(features = rownames(cml_seurat_new_filt3), approx = F) %>% RunTSNE(dims = 1:8, perplexity = 100)

DimPlot(cml_seurat_new_filt3, label = T, repel = T, group.by = "HTO_classification.global")
DimPlot(cml_seurat_new_filt3, label = T, repel = T, group.by = "my.demux")
DimPlot(cml_seurat_new_filt3, label = T, repel = T, group.by = "my.demux3")
HTOHeatmap(cml_seurat_new_filt3, assay = "HTO", ncells = 5000)

cml_seurat_new_filt <- merge(cml_seurat_new_filt1, list(cml_seurat_new_filt2, cml_seurat_new_filt3))
cml_seurat_hto      <- cml_seurat_new_filt
saveRDS(cml_seurat_hto, "results/functional/cml_seurat_hto.rds")

cml_seurat_hto <- readRDS("results/functional/cml_seurat_hto.rds")

