
###### Used in cml stop project

cml_seurat <- readRDS("results/cml_stop/cml_seurat_clustered.rds")

## Get the CML immune cells
public_immune_seurat <- readRDS("/Users/janihuuh/Dropbox/cml_stop/results/public_seurat2.rds")
cells.to.keep        <- public_immune_seurat@meta.data %>% filter(project %in% c("CML on TKI", "CML dg")) %>% pull(barcode)
cml_immune_seurat    <- subset(public_immune_seurat, cells = cells.to.keep)

Idents(cml_immune_seurat) <- cml_immune_seurat$RNA_snn_res.0.5 %>% extractClusterNumber() %>% as.numeric %>% as.factor() %>% getClusterPhenotypesPublic()
cml_immune_seurat$cluster <- Idents(cml_immune_seurat)

cells.to.keep     <- cml_immune_seurat@meta.data %>% filter(cluster != "5 CLL cells") %>% pull(barcode)
cml_immune_seurat <- subset(cml_immune_seurat, cells = cells.to.keep) %>% getLatentUMAP()

## Merge with the CML-object
comparison_cml_seurat <- merge(cml_seurat, cml_immune_seurat)

## Samples to scVI; write .csv files
clonality_genes <- getClonalityGenes(comparison_cml_seurat)

## Get features that are at least in every sample
a <- rownames(cml_seurat)
b <- rownames(cml_immune_seurat)

genes_to_use <- Reduce(intersect, list(a,b))
genes_to_use <- genes_to_use[!genes_to_use %in% clonality_genes]

diet_seurat         <- DietSeurat(comparison_cml_seurat, features = genes_to_use)
Idents(diet_seurat) <- diet_seurat$orig.ident

genes_to_keep <- list()
i <- 1
dir.create("results/scvi/input_files/cml_stop/", showWarnings=F)

for(patient in unique(comparison_cml_seurat$orig.ident)){
  
  message(patient)
  seurat_temp <- subset(diet_seurat, idents = patient)
  counts_temp <- seurat_temp@assays$RNA@counts %>% as.data.frame
  
  genes_to_keep[[i]] <- rownames(counts_temp)
  i <- i + 1
  
  counts_temp <- diet_seurat@assays$RNA@counts[ ,diet_seurat$orig.ident == patient] %>% as.data.frame
  counts_temp <- counts_temp[!rownames(counts_temp) %in% clonality_genes, ]
  
  data.table::fwrite(counts_temp, paste0("results/scvi/input_files/cml_stop/", patient, ".csv"), sep = ",", quote = F, row.names = T, col.names = T)
  
}








## Get latent representation; UMAP it
dir.create("results/scvi/results/cml_stop/", showWarnings = F)

cml_latent      <- fread("results/scvi/output/cml_stop_lsc_oneshot_latent.csv", header = F)
cml_batch       <- fread("results/scvi/output/cml_stop_lsc_oneshot_indices.csv", header = F)
cml_latent_umap <- uwot::umap(cml_latent)
cml_latent_umap <- cml_latent_umap %>% as.data.frame() %>% dplyr::rename(UMAP1=V1, UMAP2=V2)

data.frame(cml_latent_umap, batch = cml_batch$V1) %>% ggplot(aes(UMAP1, UMAP2, color = as.factor(batch))) + geom_point(size = 0.3) + add_guide + scale_color_manual(values = getPalette3(24))
ggsave("results/scvi/results/cml_stop/latent_umap_best.pdf", width = 8, height = 6)

write.table(select(cml_latent_umap, UMAP1:UMAP2), "results/scvi/output/cml_stop_lsc_oneshot_latent_umap.csv", sep = ",")

## Put embeddings into Seurat object
cml_latent      <- as.matrix(cml_latent)
cml_latent_umap <- as.matrix(dplyr::select(as.data.frame(cml_latent_umap), UMAP1:UMAP2))

rownames(cml_latent)      <- colnames(diet_seurat)
rownames(cml_latent_umap) <- colnames(diet_seurat)
rownames(cml_latent) == rownames(cml_latent_umap)

latent_dim_red            <- CreateDimReducObject(key = "latent", embeddings = as.matrix(x = cml_latent))
latent_umap_dim_red       <- CreateDimReducObject(key = "latent_umap", embeddings = as.matrix(x = cml_latent_umap))

diet_seurat[['latent']]      <- latent_dim_red
diet_seurat[['latent_umap']] <- latent_umap_dim_red

diet_seurat$cluster_new <- ifelse(!is.na(diet_seurat$cluster), as.character(diet_seurat$cluster), as.character(diet_seurat$nabo)) %>% as.character()
diet_seurat$cluster_new[diet_seurat$cluster_new == ""] <- "no_pred"

diet_seurat$cluster_new[which(diet_seurat$kit == "kitpos")] <- paste(as.character(diet_seurat$cluster_new)[which(diet_seurat$kit == "kitpos")], "kitpos")
diet_seurat$cluster_new[which(diet_seurat$kit == "kitneg")] <- paste(as.character(diet_seurat$cluster_new)[which(diet_seurat$kit == "kitneg")], "kitneg")

DimPlot(diet_seurat, group.by = "cluster_new", label = T, repel = T, cols = getPalette5(41)) + theme_bw(base_size = 12) + theme(legend.position = "none") + labs(x = "UMAP 1", y = "UMAP 2")
ggsave("results/cml_stop/umap_total.png", width = 6, height = 5)

saveRDS(diet_seurat, "results/cml_stop/diet_seurat.rds")
