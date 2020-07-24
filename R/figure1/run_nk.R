
## Focus on NK cells
cells.to.keep  <- comparison_cml_seurat@meta.data %>% filter(public_clusters %in% c("3 NK CD56dim", "15 NK adaptive", "18 NK CD56bright")) %>% pull(barcode)
diet_nk_seurat <- subset(comparison_cml_seurat, cells = cells.to.keep)
diet_nk_seurat <- diet_nk_seurat %>% getLatentUMAP() %>% fixSeurat()
diet_nk_seurat <- diet_nk_seurat %>% getLatentClustering()

q <- NULL; i <- 1
for(clustering_column in clustering_columns){
  q[[i]] <- diet_nk_seurat@meta.data[,clustering_column] %>% levels %>% length
  i <- i + 1
}

data.frame(resolution = res, nClusters = q) %>% ggplot(aes((resolution),nClusters), label=nClusters) + geom_point(shape = 21) + theme_bw()
ggsave("results/figure1/scatter_res_nCluster_nk.png", width = 4, height = 3)

Idents(diet_nk_seurat) <- diet_nk_seurat$RNA_snn_res.0.4 %>% getClusterPhenotypesNK

## Find markers
Idents(diet_nk_seurat) <- diet_nk_seurat$RNA_snn_res.0.4
nk_markers <- FindAllMarkers(diet_nk_seurat, test.use = "t")
nk_markers <- nk_markers %>% mutate(dir = ifelse(avg_logFC > 0, "up", "down"))
write.table(nk_markers, "results/figure1/deg_nk.txt", sep = "\t", quote = F, row.names = F)

## Visualize
diet_nk_seurat@meta.data %>% mutate(cluster = Idents(diet_nk_seurat)) %>%
  group_by(project, cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  ggplot(aes(project, prop, fill = cluster)) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette5(8)) + labs(x = "") + ggpubr::rotate_x_text(45)
ggsave("results/figure1/bar_nk_clusters.pdf", width = 5, height = 5)

q <- NULL
min_cells <- 50

for(i in 1:100){
  set.seed(i)
  q[[i]] <- diet_nk_seurat@meta.data %>% mutate(cluster = Idents(diet_nk_seurat)) %>%
    group_by(project) %>% sample_n(min_cells) %>%
    group_by(cluster, project) %>% summarise(n = n()) %>% mutate(prop = n / sum(n))
}

q <- q %>% rbindlist()

tot_df <- q %>% group_by(cluster, project) %>% summarise(prop = median(prop)) %>% group_by(cluster) %>% summarise(tot = sum(prop))
q %>% group_by(cluster, project) %>% summarise(prop = median(prop)) %>% left_join(tot_df, by = "cluster") %>% mutate(prop = prop / tot) %>%
  ggplot(aes(cluster, prop, fill = project)) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette3(6)) + labs(x = "") + ggpubr::rotate_x_text(45)
ggsave("results/figure1/bar_nk_clusters2.pdf", width = 5, height = 4)

DimPlot(diet_nk_seurat, group.by = "public_clusters", label = T, repel = T) + theme_bw(base_size = 12) + theme(legend.position = "none") + labs(x = "UMAP 1", y = "UMAP 2") + facets_nice +
  scale_color_manual(values = getPalette5(3))
ggsave("results/figure1/umap_nk.png", width = 5, height = 4)

DimPlot(diet_nk_seurat, group.by = "public_clusters", split.by = "project", ncol = 3, label = T, label.size = 2) + theme_bw(base_size = 12) + theme(legend.position = "none") + labs(x = "UMAP 1", y = "UMAP 2") + facets_nice +
  scale_color_manual(values = getPalette5(3))
ggsave("results/figure1/umap_nk_project.png", width = 6, height = 4)

diet_nk_seurat@reductions$latent_umap@cell.embeddings %>% as.data.frame %>% bind_cols(diet_nk_seurat@meta.data) %>%
  ggplot(aes(latent_umap_1,latent_umap_2)) + stat_density_2d(aes(fill = ..level..), geom = "polygon") + facet_wrap(~project) + scale_fill_distiller(direction=1) + theme(legend.position = "none") + theme_bw(base_size = 12) + theme(legend.position = "none") + labs(x = "UMAP 1", y = "UMAP 2") + facets_nice
ggsave("results/figure1/umap_nk_project_dens.png", width = 6, height = 4)

DimPlot(diet_nk_seurat, label = T, repel = T) + theme_bw(base_size = 12) + theme(legend.position = "none") + labs(x = "UMAP 1", y = "UMAP 2") + facets_nice +
  scale_color_manual(values = getPalette5(8))
ggsave("results/figure1/umap_nk_new_clusters.png", width = 5, height = 4)

DimPlot(diet_nk_seurat, split.by = "project", ncol = 3, label = T, label.size = 2) + theme_bw(base_size = 12) + theme(legend.position = "none") + labs(x = "UMAP 1", y = "UMAP 2") + facets_nice +
  scale_color_manual(values = getPalette5(8))
ggsave("results/figure1/umap_nk_project_new_clusters.png", width = 6, height = 4)

## Get idea of big clusters
diet_nk_seurat@meta.data %>% group_by(RNA_snn_res.0.4, public_clusters) %>% summarise(n = n()) %>% mutate(prop = n / sum(n))

## Visualize known markers
p <- DotPlot(diet_nk_seurat, features = rev(unique(c("CD3E", "NCAM1", "GZMK", "SELL", "XCL1", "XCL2", "KLRC1", "IL7R", "LTB", "FCGR3A", "GZMA", "GZMB", "GZMH", "GZMM", "KLRC2", "ZEB2", "KLF2", "PRDM1", "GZMH", "LAG3"))), cols = "RdYlBu") +
  ggpubr::rotate_x_text(angle = 90) + labs(x = "", y  = "")
ggsave(plot = p, "results/figure1/dotplot_nk_dufva_markers.png", width = 9, height = 6)

## Visualize known scores
diet_nk_seurat <- AddModuleScore(diet_nk_seurat, features = list(unique(inhibitory_long)), nbin = 10, name = "inhibitory")
diet_nk_seurat <- AddModuleScore(diet_nk_seurat, features = list( c("LAG3", "CTLA4", "PDCD1", "HAVCR2", "TIGIT")), nbin = 10, name = "exhaustion")
diet_nk_seurat <- AddModuleScore(diet_nk_seurat, features = list(cytotoxic_markers), nbin = 10, name = "cytotoxicity")

FeaturePlot(diet_nk_seurat, label = F, features = "inhibitory1", cols = c("lightgrey", "salmon"), order = T, min.cutoff = 0)
ggsave("results/figure1/umap_nk_inhibitory.png", width = 5, height = 4)

FeaturePlot(diet_nk_seurat, label = F, features = "exhaustion1", cols = c("lightgrey", "salmon"), order = T, min.cutoff = 0, max.cutoff = 1) #, max.cutoff = 0.2)
ggsave("results/figure1/umap_nk_exhaustion.png", width = 5, height = 4)

FeaturePlot(diet_nk_seurat, features = "cytotoxicity1", cols = c("lightgrey", "salmon"), order = T, min.cutoff = 0) #, max.cutoff = 0.2)
ggsave("results/figure1/umap_nk_cytotoxicity.png", width = 5, height = 4)

diet_nk_seurat@meta.data %>% mutate(cluster = Idents(diet_nk_seurat)) %>%
  ggplot(aes(cluster, cytotoxicity1, fill = cluster)) + geom_violin(draw_quantiles = 0.5) + scale_fill_manual(values = getPalette5(8)) + labs(x = "") + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) +
  ggpubr::stat_compare_means()
ggsave("results/figure1/vln_nk_cytotoxicity.pdf", width = 5, height = 4)

diet_nk_seurat@meta.data %>% mutate(cluster = Idents(diet_nk_seurat)) %>%
  ggplot(aes(cluster, inhibitory1, fill = cluster)) + geom_violin(draw_quantiles = 0.5) + scale_fill_manual(values = getPalette5(8)) + labs(x = "") + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) +
  ggpubr::stat_compare_means()
ggsave("results/figure1/vln_nk_inhibitory.pdf", width = 5, height = 4)

diet_nk_seurat@meta.data %>% mutate(cluster = Idents(diet_nk_seurat)) %>%
  ggplot(aes(cluster, exhaustion1, fill = cluster)) + geom_violin(draw_quantiles = 0.5) + scale_fill_manual(values = getPalette5(8)) + labs(x = "") + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) +
  ggpubr::stat_compare_means()
ggsave("results/figure1/vln_nk_exhaustion.pdf", width = 5, height = 4)

saveRDS(diet_nk_seurat, "results/diet_nk_seurat.rds")



