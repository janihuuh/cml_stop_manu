
#### NK
public_seurat   <- readRDS("results/public_seurat.rds")
nk_seurat <- subset(public_seurat, celltype == "NK")

## Preprocess data
nk_seurat <- nk_seurat %>% preprocessSeurat(cells.to.use = colnames(nk_seurat))

## Do scVI
nk_seurat %>% getScviInput(folder = "results/scvi/input_files/")

## Get scVI results
latents    <- fread("results/scvi/results/nk_latents.csv")
nk_seurat <- nk_seurat %>% putLatentsSeurat(latent = latents) %>% getLatentClustering() %>% fixSeurat()

## Decide on clustering
nk_seurat %>% plotClustering()
Idents(nk_seurat) <- public_seurat$RNA_snn_res.0.5 %>% as.character() %>% as.numeric() %>% as.factor()
nk_seurat$cluster <- Idents(nk_seurat)
saveRDS("results/nk_seurat.rds")


## Find DE-genes between all the clusters
all_markers <- FindAllMarkers(nk_seurat, test.use = "t", verbose = T, max.cells.per.ident = 1e3)
fwrite(all_markers, "results/manuscript/nk/all_markers_1e3.txt", sep = "\t", quote = F, row.names = F)
