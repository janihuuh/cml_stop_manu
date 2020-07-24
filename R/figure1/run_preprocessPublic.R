
## Samples to scVI; write .csv files
comparison_comparison_cml_seurat <- readRDS("results/comparison_comparison_cml_seurat.rds")

clonality_genes <- getClonalityGenes(comparison_comparison_cml_seurat)

## Get features that are at least in every sample
a <- rownames(wu_seurat)
b <- rownames(rendeiro_d0_seurat)
c <- rownames(healthy_seurat)
d <- rownames(comparison_cml_seurat)
e <- rownames(nordcomparison_cml_seurat)

genes_to_use <- Reduce(intersect, list(a,b,c,d,e))
genes_to_use <- genes_to_use[!genes_to_use %in% clonality_genes]

comparison_comparison_cml_seurat         <- DietSeurat(comparison_comparison_cml_seurat, features = genes_to_use)
Idents(comparison_comparison_cml_seurat) <- comparison_comparison_cml_seurat$orig.ident

genes_to_keep <- list()
i <- 1

for(patient in unique(comparison_comparison_cml_seurat$orig.ident)){
  
  message(patient)
  seurat_temp <- subset(comparison_comparison_cml_seurat, idents = patient)
  counts_temp <- seurat_temp@assays$RNA@counts %>% as.data.frame
  
  genes_to_keep[[i]] <- rownames(counts_temp)
  i <- i + 1
  
  counts_temp <- comparison_comparison_cml_seurat@assays$RNA@counts[ ,comparison_comparison_cml_seurat$orig.ident == patient] %>% as.data.frame
  counts_temp <- counts_temp[!rownames(counts_temp) %in% clonality_genes, ]
  
  data.table::fwrite(counts_temp, paste0("results/scvi/input_files/comparison/", patient, ".csv"), sep = ",", quote = F, row.names = T, col.names = T)
  
}


## From scVI; get latent representation; UMAP it
public_indices <- fread("results/scvi/output/public_oneshot_indices2.csv", header = F)
public_latent  <- fread("results/scvi/output/public_oneshot_latent2.csv", header = F)

public_latent_umap <- uwot::umap(public_latent)
public_latent_umap <- public_latent_umap %>% as.data.frame() %>% dplyr::rename(UMAP1=V1, UMAP2=V2)

data.frame(public_latent_umap, batch = public_indices$V1) %>%
  ggplot(aes(UMAP1, UMAP2, color = as.factor(batch))) + geom_point(size = 0.3) + add_guide + scale_color_manual(values = getPalette3(150)) + theme(legend.position = "none")
ggsave("results/scvi/results/latent_umap_immune_best.pdf", width = 8, height = 6)


## Add latent UMAPs to Seurat
public_latent                <- as.matrix(public_latent)
public_latent_umap           <- as.matrix(public_latent_umap)
rownames(public_latent)      <- colnames(comparison_comparison_cml_seurat)
rownames(public_latent_umap) <- colnames(comparison_comparison_cml_seurat)
latent_dim_red               <- CreateDimReducObject(key = "latent", embeddings = as.matrix(x = public_latent))
latent_umap_dim_red          <- CreateDimReducObject(key = "latent_umap", embeddings = as.matrix(x = public_latent_umap))
comparison_comparison_cml_seurat[['latent']]      <- latent_dim_red
comparison_comparison_cml_seurat[['latent_umap']] <- latent_umap_dim_red


## Cluster the cells in latent space
comparison_comparison_cml_seurat <- getLatentClustering(comparison_comparison_cml_seurat)

## Plot clustering results
p <- NULL
i <- 1

for(clustering_column in clustering_columns){
  
  nColors <- comparison_comparison_cml_seurat@meta.data[,clustering_column] %>% levels %>% length
  
  p[[i]] <- DimPlot(comparison_comparison_cml_seurat, reduction = "latent_umap", group.by = clustering_column, cols = getPalette(nColors), label = T) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(title = clustering_column)
  i <- i + 1
  
}

png("results/public/latent_umap_per_res.png", width = 1024, height = 1024)
do.call(grid.arrange, c(p, ncol = 4))
dev.off()


q <- NULL; i <- 1

for(clustering_column in clustering_columns){
  q[[i]] <- comparison_comparison_cml_seurat@meta.data[,clustering_column] %>% levels %>% length
  i <- i + 1
}

data.frame(resolution = res, nClusters = q) %>%
  ggplot(aes((resolution),nClusters), label=nClusters) + geom_point(shape = 21) + theme_bw()
ggsave("results/public/scatter_res_nClusters.png", width = 4, height = 3)

## Choose res 0.5 as optimal
p <- DimPlot(comparison_comparison_cml_seurat, reduction = "latent_umap", group.by = "RNA_snn_res.0.5", cols = getPalette(25), label = T, repel = T) + theme(legend.position = "none")
ggsave(plot = p, "results/public/latent_umap.png", width = 8, height = 7)




### Add SingleR
dir.create("results/figure1/singler/", showWarnings = F)
require(SingleR)

for(patient in unique(public_seurat$orig.ident)){
  
  message(patient)
  
  ## Cells to select
  cells.to.keep <- public_seurat@meta.data[public_seurat$orig.ident == patient, ] %>% pull(barcode)
  public_temp      = subset(public_seurat, cells = cells.to.keep)
  
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
  
  saveRDS(singler, file = paste0('results/figure1/singler/public_temp', patient, '.rds'))
  
}


## Combine the data sets
singler.objects.file <- list.files('results/figure1/singler', pattern = 'public_temp', full.names = T)
singler.objects.file <- grep(".rds", singler.objects.file, value = T)
singler.objects      <- lapply(singler.objects.file, FUN = function(x) {message(x); readRDS(x)})

singler = SingleR.Combine(singler.objects,
                          order    = colnames(comparison_comparison_cml_seurat),
                          clusters = Idents(comparison_comparison_cml_seurat),
                          xy       = comparison_comparison_cml_seurat@reductions$umap@cell.embeddings)
# saveRDS(singler, 'results/singler/cml_singler.rds')
# singler <- readRDS('results/singler/cml_singler.rds')

metadata <- data.frame(barcode        = names(singler$meta.data$orig.ident),
                       hpca_pred      = singler$singler[[1]]$SingleR.single$labels,
                       blueprint_pred = singler$singler[[2]]$SingleR.single$labels) %>% add_rownames(var = "barcode")

public_metadata <- merge(comparison_comparison_cml_seurat@meta.data, metadata, by = "barcode")


## Add into Seurat
public_metadata                      <- public_metadata[match(colnames(comparison_comparison_cml_seurat), public_metadata$barcode), ]
comparison_comparison_cml_seurat$singler_hpca_pred      <- public_metadata$hpca_pred
comparison_comparison_cml_seurat$singler_blueprint_pred <- public_metadata$blueprint_pred


## Detect doublets
require(scds)

# Annotate doublet using co-expression based doublet scoring:
cml_sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = comparison_cml_seurat@assays$RNA@counts), colData = comparison_cml_seurat@meta.data)
cml_sce <- cxds(cml_sce)
cml_sce <- bcds(cml_sce)
cml_sce <- cxds_bcds_hybrid(cml_sce)

## Add into Seurat
comparison_cml_seurat$cxds_doublet_score   <- SingleCellExperiment::colData(cml_sce)$cxds_score
comparison_cml_seurat$bcds_doublet_score   <- SingleCellExperiment::colData(cml_sce)$bcds_score

comparison_cml_seurat$hybrid_doublet_score <- SingleCellExperiment::colData(cml_sce)$hybrid_score
comparison_cml_seurat$cxds_doublet_score_norm <- c(SingleCellExperiment::colData(cml_sce)$cxds_score - min(SingleCellExperiment::colData(cml_sce)$cxds_score)) / max(SingleCellExperiment::colData(cml_sce)$cxds_score)
comparison_cml_seurat$bcds_doublet_score_norm <- c(SingleCellExperiment::colData(cml_sce)$bcds_score - min(SingleCellExperiment::colData(cml_sce)$bcds_score)) / max(SingleCellExperiment::colData(cml_sce)$bcds_score)

saveRDS(comparison_comparison_cml_seurat, "results/comparison_comparison_cml_seurat_scvi.rds")