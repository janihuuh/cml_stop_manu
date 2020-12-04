
#### CML-stop
public_seurat   <- readRDS("results/public_seurat.rds")
cml_stop_seurat <- subset(public_seurat, project == "CML stop")

## Preprocess data
cml_stop_seurat <- cml_stop_seurat %>% preprocessSeurat(cells.to.use = colnames(cml_stop_seurat))

## Do scVI
cml_stop_seurat %>% getScviInput(folder = "results/scvi/input_files/")

## Get scVI results
latents    <- fread("results/scvi/results/cml_stop_latents.csv")
cml_stop_seurat <- cml_stop_seurat %>% putLatentsSeurat(latent = latents) %>% getLatentClustering() %>% fixSeurat()

## Decide on clustering
cml_stop_seurat %>% plotClustering()
Idents(cml_stop_seurat) <- public_seurat$RNA_snn_res.0.5 %>% as.character() %>% as.numeric() %>% as.factor()
cml_stop_seurat$cluster <- Idents(cml_stop_seurat)
saveRDS("results/cml_stop_seurat.rds")


## Find DE-genes between all the clusters
all_markers <- FindAllMarkers(cml_stop_seurat, test.use = "t", verbose = T, max.cells.per.ident = 1e3)
fwrite(all_markers, "results/manuscript/cml_stop/all_markers_1e3.txt", sep = "\t", quote = F, row.names = F)
