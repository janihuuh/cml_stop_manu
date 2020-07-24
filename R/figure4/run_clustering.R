
dir.create("results/scvi/input_files/")

## Run scvi
cml_seurat$ident <- paste0(cml_seurat$patient, "_", cml_seurat$type)
Idents(cml_seurat) <- cml_seurat$ident

genes_to_keep <- list()
i <- 1

Idents(cml_seurat)

for(patient in unique(cml_seurat$sample)){
  
  message(patient)
  seurat_temp <- subset(cml_seurat, idents = patient)
  counts_temp <- seurat_temp@assays$RNA@counts %>% as.data.frame
  
  data.table::fwrite(counts_temp, paste0("results/scvi/input_files/", patient, ".csv"), sep = ",", quote = F, row.names = T, col.names = T)
  
  genes_to_keep[[i]] <- rownames(counts_temp)
  i <- i + 1
}



## Get latent representation; UMAP it
lsc_indices <- fread("results/scvi/output/lsc_optim_batch_indices.csv", header = F)
lsc_latent  <- fread("results/scvi/output/lsc_optim_batch_latent.csv", header = F)

lsc_latent_umap <- uwot::umap(lsc_latent)
lsc_latent_umap <- lsc_latent_umap %>% as.data.frame() %>% dplyr::rename(UMAP1=V1, UMAP2=V2)

data.frame(lsc_latent_umap, batch = lsc_indices$V1) %>% ggplot(aes(UMAP1, UMAP2, color = as.factor(batch))) + geom_point(size = 1) + add_guide + scale_color_manual(values = getPalette3(5)) + theme(legend.position = "none") + facet_wrap(~as.factor(batch),ncol=3)
ggsave("results/scvi/results/latent_umap_immune_best.pdf", width = 8, height = 6)


## Add to Seurat
cml_seurat$sample <- cml_seurat$orig.ident
lsc_latent      <- as.matrix(lsc_latent)
lsc_latent_umap <- as.matrix(lsc_latent_umap)

rownames(lsc_latent)      <- colnames(cml_seurat)
rownames(lsc_latent_umap) <- colnames(cml_seurat)

table(rownames(lsc_latent) == rownames(lsc_latent_umap))

latent_dim_red            <- CreateDimReducObject(key = "latent", embeddings = as.matrix(x = lsc_latent))
latent_umap_dim_red       <- CreateDimReducObject(key = "latent_umap", embeddings = as.matrix(x = lsc_latent_umap))

cml_seurat[['latent']]      <- latent_dim_red
cml_seurat[['latent_umap']] <- latent_umap_dim_red

DimPlot(cml_seurat, reduction = "latent_umap", group.by = "RNA_snn_res.0.6")
DimPlot(cml_seurat, reduction = "latent_umap", group.by = "nabo", split.by = "nabo", ncol = 3) + theme_bw(base_size = 12)
DimPlot(cml_seurat, reduction = "latent_umap", group.by = "sample")
FeaturePlot(cml_seurat, features = "bcrabl1", reduction = "latent_umap", cols = c("lightgrey", "salmon"))



## Cluster the cells
res                <- c(seq(0.1, 2, 0.1), 2.5, 3)
cml_seurat         <- FindNeighbors(cml_seurat, reduction = "latent", dims = 1:30)
cml_seurat         <- FindClusters(object = cml_seurat, resolution = res)
clustering_columns <- grep("res", colnames(cml_seurat@meta.data), value = T)
clustering_columns <- clustering_columns[order(substr(clustering_columns, 13, nchar(clustering_columns)) %>% as.numeric())]


q <- NULL; i <- 1

for(clustering_column in clustering_columns){
  
  message(clustering_column)
  q[[i]] <- cml_seurat@meta.data[,clustering_column] %>% levels %>% length
  i <- i + 1
  
}
data.frame(resolution = res, nClusters = do.call(q,what="c")) %>% 
  ggplot(aes((resolution),nClusters), label=nClusters) + geom_point(shape = 21) + theme_bw()
ggsave("results/public/scatter_res_nClusters.png", width = 4, height = 3)

Idents(cml_seurat) <- cml_seurat$RNA_snn_res.0.3 %>% getClusterPhenotypesLSC
cml_seurat$cluster <- cml_seurat$RNA_snn_res.0.3 %>% getClusterPhenotypesLSC

getClusterPhenotypesLSC <- function(clusters){
  
  clusters <- plyr::revalue(clusters, replace = c(
    
    "0"  = "0 MEP CD38+",
    "1"  = "1 MEP CD38+ cycling" ,
    "2"  = "2 BCRABL1- CD38- HSC" ,
    "3"  = "3 BCRABL1+ GMP cycling" ,
    "4"  = "4 BCRABL1+ HSC/CMP cycling" ,
    "5"  = "5 BCRABL1- CD38- HSC" ,
    "6"  = "6 BCRABL1+ GMP cycling" ,
    "7"  = "7 MEP low quality" ))
  
  return(clusters)
  
}        

cml_seurat <- RunUMAP(cml_seurat, reduction = "latent", dims = 1:30, reduction.name = "latent_umap_2", reduction.key = "latent ", n.neighbors = 10)

DimPlot(cml_seurat, reduction = "latent_umap", group.by = "nabo")
DimPlot(cml_seurat, reduction = "latent_umap_2", group.by = "nabo")

DimPlot(cml_seurat, reduction = "latent_umap", group.by = "RNA_snn_res.0.3")
DimPlot(cml_seurat, reduction = "latent_umap_2", group.by = "RNA_snn_res.0.3", label = T, repel = T)

DimPlot(cml_seurat, reduction = "latent_umap_2", label = T, repel = T)

Idents(cml_seurat) <- cml_seurat$RNA_snn_res.0.3 %>% getClusterPhenotypesLSC
cml_seurat$cluster <- cml_seurat$RNA_snn_res.0.3 %>% getClusterPhenotypesLSC

markers <- FindAllMarkers(cml_seurat, test = "t") %>% filter(p_val_adj < 0.05 & avg_logFC > 0)
markers %>% filter(gene %in% do.call(van_galen_markers, what = "c"))  
cml_seurat@meta.data %>% group_by(nabo, cluster) %>% summarise(n=n()) %>% mutate(prop=n/sum(n)) %>% View







##### Find the BCR-ABL1 positive cells

## Giustacchini
gisutacchini_deg <- fread("/Users/janihuuh/Dropbox/cml_stop/data/bcr_abl_pathway/giustacchini_deg_ablpos_neg.txt") %>% filter(Log2FC > 0) %>% pull(Gene)
gisutacchini_deg <- gisutacchini_deg[gisutacchini_deg %in% rownames(cml_seurat)]

cml_seurat <- AddModuleScore(cml_seurat, features = list(gisutacchini_deg), name = "bcrabl")
median_df  <- cml_seurat@meta.data %>% group_by(cluster) %>% summarise(median = median(bcrabl1))

cml_seurat@meta.data %>% 
  left_join(median_df) %>% 
  ggplot(aes(reorder(cluster,median),bcrabl1,fill=cluster)) + geom_violin(draw_quantiles = c(0.25,0.5,0.75)) + ggpubr::rotate_x_text(angle = 90) + scale_fill_manual(values = getPalette3(13)) + theme(legend.position = "none") + labs(x = "", y = "BCR-ABL1 score") + ggpubr::stat_compare_means()
ggsave("results/cml_stop/vln_bcrabl.pdf", width = 5, height = 4)

cml_seurat$bcrabl1_pos_cells <- ifelse(cml_seurat$bcrabl1 > median(cml_seurat$bcrabl1), "above median", "below median")
DimPlot(cml_seurat, group.by = "bcrabl1_pos_cells", cols = c("salmon", "lightgrey"), reduction = "latent_umap_2") + theme_bw(base_size = 12) + labs(x = "UMAP 1", y = "UMAP 2")
ggsave("results/cml_stop/umap_bcrablmedian.png", width = 5, height = 4)


## Singler
pred.blue <- SingleR(test = cml_seurat@assays$RNA@counts, ref = blueprint, labels =  blueprint$label.fine)
pred.hesc <- SingleR(test = cml_seurat@assays$RNA@counts, ref = hpca.se, labels =  hpca.se$label.main)

cml_seurat$pred_hcpa <- pred.hesc$first.labels
cml_seurat$pred_blue <- pred.blue$first.labels

DimPlot(cml_seurat, group.by = "pred_hcpa", reduction = "latent_umap_2")
DimPlot(cml_seurat, group.by = "pred_blue", reduction = "latent_umap_2")
DimPlot(cml_seurat, group.by = "nabo", reduction = "latent_umap_2",label=T)




DimPlot(cml_seurat, reduction = "latent_umap_2", cols = getPalette3(14), label = T, repel = T)  + theme_bw(base_size = 12) + labs(x = "UMAP 1", y = "UMAP 2") + theme(legend.position = "none")

saveRDS(cml_seurat, "results/cml_stop/cml_seurat_clustered.rds")
# cml_seurat <- readRDS("results/cml_stop/cml_seurat_clustered.rds")



## PR1-genes
cml_seurat <- AddModuleScore(cml_seurat, features = c("PRTN3", "ELANE"), name = "PR1_score")

FeaturePlot(cml_seurat, reduction = "latent_umap_2", features = c("PR1_score1", "PRTN3", "ELANE"), order = T, cols = c("lightgrey", "tomato"), ncol = 3)
ggsave("results/cml_stop/umap_pr1.png", width = 9, height = 3)

FeaturePlot(cml_seurat, reduction = "latent_umap_2", label.size = 3, label = T, repel = T, features = c("PR1_score1"), order = T, cols = c("lightgrey", "tomato")) + theme_bw() + labs(title = "", x = "UMAP 1", y = "UMAP 2", color = "PR1 score")
ggsave("results/cml_stop/umap_pr1.png", width = 5, height = 4)

cml_seurat@meta.data %>% 
  ggplot(aes(bcrabl1,PR1_score1)) + geom_point() + scale_y_log10() + scale_x_log10() + ggpubr::stat_cor() + geom_smooth(method="lm")


cml_seurat@meta.data %>% ggplot(aes(bcrabl1_quantiles,PR1_score1)) + geom_boxplot(draw_quantiles = c(0.25,0.5,0.75), adjust = 3, outlier.shape = NA) + geom_jitter() + ggpubr::stat_compare_means()

cml_seurat@meta.data %>% ggplot(aes(bcrabl1_pos_cells,PR1_score1)) + geom_boxplot(draw_quantiles = c(0.25,0.5,0.75), adjust = 3, outlier.shape = NA) + geom_jitter() + ggpubr::stat_compare_means()

table(cml_seurat$PR1_score1<0.5,cml_seurat$bcrabl1_pos_cells) %>% fisher.test()




