
dir.create("results/relapse/", showWarnings = F)

## Analyze relapse samples, compare them to dg-samples
diet_seurat    <- readRDS("results/public_seurat2.rds")
cells.to.keep  <- diet_seurat@meta.data %>% filter(timepoint == "relapse" | project == "CML dg") %>% pull(barcode)
relapse_seurat <- subset(diet_seurat, cells = cells.to.keep)
relapse_seurat <- relapse_seurat %>% getLatentUMAP() %>% fixSeurat()
relapse_seurat <- relapse_seurat %>% getLatentClustering()
saveRDS(diet_seurat, "results/relapse_seurat.rds")

relapse_seurat <- readRDS("results/relapse_seurat.rds")
DimPlot(relapse_seurat, group.by = "public_clusters", label = T, repel = T) + theme(legend.position = "none")

res_new <- c(seq(0.1, 1, 0.1), seq(1.2, 2, 0.2), 2.5, 3)
res_old <- c(seq(0.1, 2, 0.1), 2.5, 3)

clustering_columns <- grep("res", colnames(relapse_seurat@meta.data), value = T)
clustering_columns <- clustering_columns[substr(clustering_columns, 13, nchar(clustering_columns)) %in% res_new]
clustering_columns <- clustering_columns[order(substr(clustering_columns, 13, nchar(clustering_columns)) %>% as.numeric())]
# clustering_columns <- clustering_columns[res_old %in% res_new]


## Plot clustering results
p <- NULL
i <- 1

for(clustering_column in clustering_columns){

  message(clustering_column)
  nColors <- relapse_seurat@meta.data[,clustering_column] %>% levels %>% length

  p[[i]] <- DimPlot(relapse_seurat, reduction = "latent_umap", group.by = clustering_column, cols = getPalette(nColors), label = T) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(title = clustering_column)
  i <- i + 1

}

png("results/relapse/latent_umap_per_res.png", width = 1024, height = 1024)
do.call(grid.arrange, c(p, ncol = 4))
dev.off()


q <- NULL
i <- 1

for(clustering_column in clustering_columns){

  message(clustering_column)
  q[[i]] <- relapse_seurat@meta.data[,clustering_column] %>% levels %>% length
  i <- i + 1

}

data.frame(resolution = res_new, nClusters = q) %>%
  ggplot(aes((resolution),nClusters), label=nClusters) + geom_point(shape = 21) + theme_bw()
ggsave("results/relapse/scatter_res_nClusters.png", width = 4, height = 3)


## Choose 0.5
p <- DimPlot(relapse_seurat, reduction = "latent_umap", group.by = "RNA_snn_res.0.5", cols = getPalette4(20), label = T, repel = T) + theme(legend.position = "none")
ggsave(plot = p, "results/relapse/latent_umap.png", width = 8, height = 7)




## Add the clustering info
clusters <- unique(relapse_seurat$saved.idents) %>% extractClusterNumber() %>% as.numeric() %>% order
relapse_seurat$saved.idents <- factor(as.character(relapse_seurat$saved.idents), levels = unique(relapse_seurat$saved.idents)[clusters])
relapse_seurat$saved.idents <- relapse_seurat$saved.idents %>% extractClusterNumber() %>% as.numeric() %>% as.factor() %>% getClusterPhenotypesCml()

p <- DimPlot(relapse_seurat, group.by = "saved.idents", split.by = "saved.idents", ncol = 5, cols = getPalette3(25), pt.size = 0.0001) + theme_bw(base_size = 12) + theme(legend.position = "none") + labs(x = "UMAP 1", y = "UMAP 2")
ggsave(plot = p, "results/relapse/latent_umap_per_old_cluster.png", width = 7, height = 7)


## Coarsely rename the clusters
relapse_seurat@meta.data %>%
  group_by(RNA_snn_res.0.5, saved.idents) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% group_by(RNA_snn_res.0.5) %>% top_n(n = 3) %>% arrange(as.numeric(as.character(RNA_snn_res.0.5))) %>% View

relapse_seurat$public_clusters <- relapse_seurat$RNA_snn_res.0.5 %>% reorderClusters %>% getClusterPhenotypesRelapse()
DimPlot(relapse_seurat, reduction = "latent_umap", group.by = "public_clusters", cols = getPalette4(20), label = T, repel = T) + theme(legend.position = "none") + theme_bw(base_size = 12) + theme(legend.position = "none") + labs(x = "UMAP 1", y = "UMAP 2")
ggsave("results/relapse/latent_umap_new_cluster.png", width = 6, height = 5)

## Sample to same amount of cells
relapse_seurat@meta.data %>% group_by(project) %>% summarise(n = n())
cells.to.keep <- relapse_seurat@meta.data %>% group_by(project) %>% sample_n(7622) %>% pull(barcode)
temp_seurat <- subset(relapse_seurat, cells = cells.to.keep)
temp_seurat <- temp_seurat %>% getLatentUMAP() %>% fixSeurat()

DimPlot(temp_seurat, reduction = "latent_umap", group.by = "public_clusters", cols = getPalette4(20), split.by = "project", label = T, repel = T, label.size = 2.5) + theme(legend.position = "none") + theme_bw(base_size = 12) + theme(legend.position = "none") + labs(x = "UMAP 1", y = "UMAP 2")
ggsave("results/relapse/latent_umap_new_cluster_project.png", width = 12, height = 4)

DimPlot(temp_seurat, reduction = "latent_umap", group.by = "public_clusters", cols = getPalette4(20), label = T, repel = T, label.size = 2.5) + theme(legend.position = "none") + theme_bw(base_size = 12) + theme(legend.position = "none") + labs(x = "UMAP 1", y = "UMAP 2")
ggsave("results/relapse/latent_umap_new_cluster_project_solo.png", width = 5, height = 4)


temp_seurat@reductions$latent_umap@cell.embeddings %>% as.data.frame %>% bind_cols(temp_seurat@meta.data) %>%
  ggplot(aes(latent_umap_1,latent_umap_2)) + stat_density_2d(aes(fill = ..level..), geom = "polygon") + facet_wrap(~project) + scale_fill_distiller(direction=1) + theme(legend.position = "none") + theme_bw(base_size = 12) + theme(legend.position = "none") + labs(x = "UMAP 1", y = "UMAP 2")
ggsave("results/relapse/latent_umap_new_cluster_project_dens.png", width = 12, height = 4)

temp_seurat@reductions$latent_umap@cell.embeddings %>% as.data.frame %>% bind_cols(temp_seurat@meta.data) %>%
  ggplot(aes(latent_umap_1,latent_umap_2)) + stat_density_2d(aes(fill = ..level..), geom = "polygon") + facet_wrap(~project, ncol = 1) + scale_fill_distiller(direction=1) + theme(legend.position = "none") + theme_bw(base_size = 12) + theme(legend.position = "none") + labs(x = "UMAP 1", y = "UMAP 2")
ggsave("results/relapse/latent_umap_new_cluster_project_dens.png", width = 4, height = 10)



Idents(relapse_seurat) <- relapse_seurat$public_clusters
p <- DotPlot(relapse_seurat, features = rev(unique(big_markers)), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y = "")
ggsave(plot = p, "results/relapse/dotplot_big_markers.png", width = 14, height = 6)



## Find project specific clusters
relapse_seurat@meta.data$project <- droplevels(relapse_seurat@meta.data$project)
relapse_seurat$project <- as.character(relapse_seurat$project)
relapse_seurat$project <- ifelse(relapse_seurat$project == "CML dg", "diagnosis", as.character(relapse_seurat$relapse)) %>% as.character() # as.factor()
relapse_seurat$project[is.na(relapse_seurat$project)] <- "Diagnosis"
relapse_seurat$project <- as.factor(relapse_seurat$project)

set.seed(123)
min_cells <- min(table(relapse_seurat@meta.data$project))
min_cells <- 5e3
relapse_seurat@meta.data %>% group_by(project) %>% sample_n(min_cells) %>% group_by(project) %>%
  group_by(public_clusters, project) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  ggplot(aes(public_clusters, prop, fill = project)) + geom_bar(stat = "identity") + ggpubr::rotate_x_text(angle = 45) + scale_fill_manual(values = getPalette3(6)[c(3,2,1,4:6)]) + labs(x = "") +
  geom_hline(yintercept = 1/3, linetype = "dotted") +
  geom_hline(yintercept = 2/3, linetype = "dotted")
ggsave("results/relapse/bar_new_cluster_project_resampled.png", width = 8, height = 4)



relapse_seurat@meta.data %>% group_by(project) %>% sample_n(min_cells) %>% group_by(project) %>%
  group_by(project, public_clusters) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  ggplot(aes(project, prop, fill = public_clusters)) + geom_bar(stat = "identity") + ggpubr::rotate_x_text(angle = 45) + scale_fill_manual(values = getPalette3(20)) + labs(x = "") +
  geom_hline(yintercept = 1/3, linetype = "dotted") +
  geom_hline(yintercept = 2/3, linetype = "dotted")
ggsave("results/relapse/bar_new_cluster_project_resampled2.png", width = 7, height = 4)

relapse_df <- relapse_seurat@meta.data %>% group_by(project, orig.ident) %>% summarise(n = n()) %>% ungroup() %>% dplyr::select(-n)

relapse_seurat@meta.data %>% group_by(orig.ident) %>% sample_n(3e3) %>%
  group_by(orig.ident, public_clusters) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  group_by(orig.ident) %>% summarise(gini = ineq::Gini(prop), diversity = vegan::diversity(prop)) %>% left_join(relapse_df) %>%
  ggplot(aes(gini,diversity,label=orig.ident, fill = project)) + geom_point(shape = 21, size = 3) + ggrepel::geom_text_repel() + add_guide
ggsave("results/relapse/scatter_diversity_sample.pdf", width = 5, height = 4)

relapse_seurat@meta.data %>% group_by(orig.ident) %>% sample_n(3e3) %>%
  group_by(orig.ident, public_clusters) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = "", prop, fill = public_clusters)) + geom_bar(stat = "identity") + facet_wrap(~orig.ident, ncol = 4) + coord_polar(theta="y") + scale_fill_manual(values = getPalette3(20)) + labs(x = "") + theme_void() +
  theme(legend.position = "bottom")
ggsave("results/relapse/pie_new_cluster_project_resampled2.png", width = 7, height = 4)



## Find project specific clusters
q <- NULL

for(i in 1:100){
  set.seed(i)
  q[[i]] <- relapse_seurat@meta.data %>% group_by(project) %>% sample_n(min_cells) %>% group_by(project) %>%
    group_by(public_clusters, project) %>% summarise(n = n()) %>% mutate(prop = n / sum(n))
}

q <- q %>% rbindlist()

q %>%
  ggplot(aes(project, prop)) + geom_violin(aes(fill = project)) + ggpubr::rotate_x_text(angle = 45) + scale_fill_manual(values = getPalette3(6)[c(3,2,1,4:6)]) + labs(x = "") +
  geom_hline(yintercept = 0.75, linetype = "dotted") +
  geom_hline(yintercept = 0.50, linetype = "solid") +
  geom_hline(yintercept = 0.25, linetype = "dotted") +
  facet_wrap(~public_clusters) + ggpubr::stat_compare_means(ref = "healthy", label = "p.signif") + facets_nice + theme(legend.position = "none")
ggsave("results/relapse/vln_new_cluster_resampled_project_100.png", width = 7, height = 8)

df_dg     <- q %>% group_by(public_clusters, project) %>% summarise(median = median(prop)) %>% filter(project == "Healthy")
df_on_tki <- q %>% group_by(public_clusters, project) %>% summarise(median = median(prop)) %>% filter(project == "CML on TKI")

df_on_tki %>% left_join(df_dg, by = "public_clusters") %>% mutate(log2fc = log2(median.x / median.y)) %>%
  ggplot(aes(reorder(public_clusters, log2fc), log2fc, fill = ifelse(log2fc > 0, "y", "n"))) + geom_bar(stat = "identity") + coord_flip() + ylim(c(-8,8)) + labs(x = "") +
  geom_hline(yintercept = -1, linetype = "dotted") + geom_hline(yintercept = 1, linetype = "dotted") + scale_fill_manual(values = getPalette3(6)[c(3,2,1,4:6)]) + theme_classic(base_size = 12) + theme(legend.position = "none")
ggsave("results/relapse/bar_tki_fc.pdf", width = 5, height = 4)

patient_df <- relapse_seurat@meta.data %>% group_by(orig.ident, project) %>% summarise(n = n()) %>% dplyr::select(-n)

relapse_seurat@meta.data %>%
  group_by(orig.ident, public_clusters) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% left_join(patient_df) %>%
  # mutate(project = factor(project, levels = c("CML on TKI", "Healthy"))) %>%
  ggplot(aes(project,prop, fill = project)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.5) + facet_wrap(~public_clusters, scales = "free_y") + ggpubr::stat_compare_means(label = "p.format") + facets_nice +
  scale_fill_manual(values = getPalette3(6)[c(3,2,1,4:6)]) + theme_classic(base_size = 12) + theme(legend.position = "none") + labs(x = "") + ggpubr::rotate_x_text(45)
ggsave("results/relapse/box_tki.pdf", width = 7, height = 6)


## Get the DEGs
relapse_seurat$timepoint <- relapse_seurat$project
relapse_seurat$timepoint <- plyr::revalue(relapse_seurat$timepoint, replace = c("Slow" = "12m", "Fast" = "6m", "diagnosis" = "baseline"))
DEG_rp_cluster    <- lapply(levels(Idents(relapse_seurat)), getDEGbyCluster, seurat_object = relapse_seurat, min_cells = 5)
DEG_rp_cluster_df <- DEG_tki_cluster %>% rbindlist()
write.table(DEG_rp_cluster_df, "results/relapse/deg_tki.txt", sep = "\t", quote = F, row.names = F)


## Get enrichment
output_dir <- "results/relapse/enrichment/"
dir.create(output_dir, showWarnings = F)
universe_df <- rownames(relapse_seurat)

lapply(levels(Idents(relapse_seurat)), function(x){

  message(x)
  p <- deg_rp_cluster_df %>% filter(timepoint == "2v1" & cluster == x & direction == "up") %>% plotHypergeometric(universe_df = universe_df, term_df = hallmark)
  x <- gsub(" ", "_", x)
  x <- gsub("\\/", "slash", x)

  ggsave(plot = p, paste0(output_dir, x, "_hallmark_up_2v1", ".png"), width = 24, height = 12)

})

lapply(levels(Idents(relapse_seurat)), function(x){

  message(x)
  p <- deg_rp_cluster_df %>% filter(timepoint == "2v1" & cluster == x & direction == "down") %>% plotHypergeometric(universe_df = universe_df, term_df = hallmark)
  x <- gsub(" ", "_", x)
  x <- gsub("\\/", "slash", x)
  ggsave(plot = p, paste0(output_dir, x, "_hallmark_down_2v1", ".png"), width = 24, height = 12)

})


lapply(levels(Idents(relapse_seurat)), function(x){

  message(x)
  p <- deg_rp_cluster_df %>% filter(timepoint == "2v1" & cluster == x & direction == "up") %>% plotHypergeometric(universe_df = universe_df, term_df = go)
  x <- gsub(" ", "_", x)
  x <- gsub("\\/", "slash", x)

  ggsave(plot = p, paste0(output_dir, x, "_go_up_2v1", ".png"), width = 24, height = 12)

})

lapply(levels(Idents(relapse_seurat)), function(x){

  message(x)
  p <- deg_rp_cluster_df %>% filter(timepoint == "2v1" & cluster == x & direction == "down") %>% plotHypergeometric(universe_df = universe_df, term_df = go)
  x <- gsub(" ", "_", x)
  x <- gsub("\\/", "slash", x)
  ggsave(plot = p, paste0(output_dir, x, "_go_down_2v1", ".png"), width = 24, height = 12)

})


## Volcano
volcano_df <- deg_rp_cluster_df %>% mutate(cluster = reorderClusters(cluster))
ggplot(volcano_df, aes(avg_logFC, -log10(p_val), fill = ifelse(avg_logFC > 0, "y", "n"))) + geom_point(shape = 21, alpha = 0.5) + geom_vline(xintercept = -0.1, linetype = "dotted") + geom_vline(xintercept = 0.1, linetype = "dotted") +
  theme(legend.position = "None") + facet_wrap(~cluster) + xlim(values = c(-2,2)) +
  #  ggrepel::geom_text_repel(data = subset(volcano_df, gene %in% inhibitory_long), aes(label=gene, color = ifelse(avg_logFC > 0, "y", "n")), fontface = 3) +
  scale_fill_manual(values = c("dodgerblue", "salmon")) +
  scale_color_manual(values = c("dodgerblue", "salmon")) + facets_nice
ggsave("results/relapse/volcano_tki.png", width = 7, height = 6)

deg_rp_cluster_df %>% group_by(cluster) %>% summarise(n = n()) %>% mutate(celltype = extractCoarsePhenotype(cluster)) %>%
  mutate(celltype = ifelse(celltype %in% c("innate", "low", "T", "CML", "pDC"), "other", celltype)) %>%
  ggplot(aes(reorder(cluster, n),n,fill=celltype)) + geom_bar(color = "grey", stat = "identity") + coord_flip() + labs(x = "", y = "nDEGs") + scale_fill_manual(values = getPalette4(8)) + theme_classic(base_size = 12)
ggsave("results/relapse/bar_deg_tki.pdf", width = 6, height = 4)


deg_rp_cluster_df %>% filter(gene %in% cytotoxic_markers) %>% filter(direction == "down")


## Cytotoxicity
relapse_seurat <- AddModuleScore(relapse_seurat, features = cytotoxic_markers, name = "cytotoxic")

relapse_seurat@meta.data %>%
  filter(public_clusters %in% c("7 T cells activated", "4 NK adaptive", "3 NK CD56dim", "2 CD8 effectory", "2 CD8 effectory",
                                "11 CD8 EM", "6 CD8 CM/naive", "14 NK CD56bright", "12 innate lymphocytes / low quality")) %>%

  ggplot(aes(project, cytotoxic1, fill = project)) + geom_violin(draw_quantiles = 0.5, adjust = 15) + facet_wrap(~public_clusters, scales = "free_y") + ggpubr::stat_compare_means(label = "p.format") +
  theme_classic(base_size = 12) + ggpubr::rotate_x_text(45) + theme(legend.position = "none") + facets_nice + scale_fill_manual(values = getPalette4(6)) + labs(x = "", y = "cytotoxicity score")
ggsave("results/relapse/vln_cytotoxicity.pdf", width = 6, height = 5)




## calcineurin activitiy?
calcineurin_pathway <- fread("/Users/janihuuh/Dropbox/aplastic_anemia_sc/results/scrnaseq/pathways/go_calcineurin.txt")
calcineurin_genes   <- rownames(relapse_seurat)[calcineurin_pathway$V1 %in% rownames(relapse_seurat)]
relapse_seurat  <- AddModuleScore(relapse_seurat, features = list(calcineurin_genes), name = "calcineurin")

unique(relapse_seurat$public_clusters)

relapse_seurat@meta.data %>%
  filter(public_clusters %in% c("0 CD4 naive", "1 CD4 CM", "2 CD8 effectory", "6 CD8 CM/naive", "7 T cells activated", "6 CD8 CM/naive", "7 T cells activated", "9 CD4 treg", "11 CD8 EM")) %>%
  ggplot(aes(project, calcineurin1, fill = project)) + geom_violin(draw_quantiles = 0.5, adjust = 2) + facet_wrap(~public_clusters, scales = "free_y") + ggpubr::stat_compare_means(label = "p.format") +
  theme_classic(base_size = 12) + ggpubr::rotate_x_text(45) + theme(legend.position = "none") + facets_nice + scale_fill_manual(values = getPalette4(6)) + labs(x = "", y = "calcineurin pathway score")
ggsave("results/relapse/vln_calcineurin.pdf", width = 6, height = 5)





### Zoom in on the CD8+ effecto exhausted cells
cells.to.keep <- relapse_seurat@meta.data %>% filter(public_clusters == "0 CD8 effectory/exhausted") %>% pull(barcode)
relapse_cd8_seurat <- subset(relapse_seurat, cells = cells.to.keep) %>% getLatentUMAP() %>% fixSeurat()
relapse_cd8_seurat <- relapse_cd8_seurat %>% getLatentClustering()

q <- NULL
i <- 1

for(clustering_column in clustering_columns){

  message(clustering_column)
  q[[i]] <- relapse_cd8_seurat@meta.data[,clustering_column] %>% levels %>% length
  i <- i + 1

}

data.frame(resolution = res_new, nClusters = q) %>%
  ggplot(aes((resolution),nClusters), label=nClusters) + geom_point(shape = 21) + theme_bw()
ggsave("results/relapse/scatter_res_nClusters_cd8.png", width = 4, height = 3)



DimPlot(relapse_cd8_seurat, label = T, repel = T, group.by = "RNA_snn_res.0.2")
DimPlot(relapse_cd8_seurat, label = T, repel = T, group.by = "RNA_snn_res.0.5")
DimPlot(relapse_cd8_seurat, label = T, repel = T, group.by = "project")

relapse_cd8_seurat@meta.data %>%
  group_by(project, RNA_snn_res.0.2) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  ggplot(aes(project, prop, fill = RNA_snn_res.0.2)) + geom_bar(stat = "identity")

relapse_cd8_seurat@meta.data %>%
  group_by(project, RNA_snn_res.0.5) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  ggplot(aes(project, prop, fill = RNA_snn_res.0.5)) + geom_bar(stat = "identity")

Idents(relapse_cd8_seurat) <- relapse_cd8_seurat$RNA_snn_res.0.5
cd8_markers <- FindAllMarkers(relapse_cd8_seurat, test.use = "t")
write.table(cd8_markers, "results/relapse/deg_cd8.txt", sep = "\t", quote = F, row.names = F)


cd8_markers %>% filter(gene %in% inhibitory_long) %>% filter(avg_logFC > 0) %>% group_by(cluster) #%>% summarise(n = n())
cd8_markers %>% filter(gene %in% cytotoxic_markers) %>% filter(avg_logFC > 0) %>% group_by(cluster) #%>% summarise(n = n())
cd8_markers %>% filter(gene %in% big_markers) %>% filter(avg_logFC > 0) %>% group_by(cluster) #%>% summarise(n = n())

cd8_markers %>% filter(cluster == "3") %>% filter(avg_logFC > 0) %>% arrange((p_val)) #summarise(n = n())
cd8_markers %>% filter(cluster == "5") %>% filter(avg_logFC > 0) %>% arrange((p_val)) #summarise(n = n())
cd8_markers %>% filter(cluster == "6") %>% filter(avg_logFC > 0) %>% arrange((p_val)) #summarise(n = n())

getClusterPhenotypesRelapseCD8 <- function(clusters){

  plyr::revalue(clusters, replace = c("0"	 = "0 Exhausted TIGIT+, \nfast-specific",
                                      "1"	 = "1 IFNg producing",
                                      "2"	 = "2 Most cytotoxic",
                                      "3"	 = "3 NFKBIA+, slow-specific",
                                      "4"	 = "4 Exhausted LAIR1+",
                                      "5"	 = "5 transitional 1",
                                      "6"	 = "6 transitional 2"))
}

Idents(relapse_cd8_seurat) <- relapse_cd8_seurat$RNA_snn_res.0.5 %>% getClusterPhenotypesRelapseCD8
saveRDS(relapse_cd8_seurat, "results/relapse_cd8_seurat.rds")

DimPlot(relapse_cd8_seurat, label = T, repel = T, cols = getPalette4(7)) + theme_bw(base_size = 12) + theme(legend.position = "none") + labs(x = "UMAP 1", y = "UMAP 2") #, group.by = "project")
ggsave("results/relapse/umap_cd8_clusters.pdf", width = 5, height = 4)

DimPlot(relapse_cd8_seurat, group.by = "project", split.by = "project", cols = getPalette3(3)) + theme_bw(base_size = 12) + theme(legend.position = "none") + labs(x = "UMAP 1", y = "UMAP 2") #, group.by = "project")
ggsave("results/relapse/umap_cd8_project.pdf", width = 12, height = 4)

relapse_cd8_seurat@meta.data %>% mutate(cluster = Idents(relapse_cd8_seurat)) %>%
  group_by(project, cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  ggplot(aes(project, prop, fill = cluster)) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette4(7)) + labs(x = "")
ggsave("results/relapse/bar_cd8_clusters.pdf", width = 5, height = 4)



DotPlot(relapse_cd8_seurat, features = c("LAG3", "CTLA4", "PDCD1", "HAVCR2", "TIGIT")) + ggpubr::rotate_x_text()
DotPlot(relapse_cd8_seurat, features = unique(inhibitory_long)) + ggpubr::rotate_x_text()

relapse_cd8_seurat <- AddModuleScore(relapse_cd8_seurat, features = list(unique(inhibitory_long)), nbin = 10, name = "inhibitory")
relapse_cd8_seurat <- AddModuleScore(relapse_cd8_seurat, features = list( c("LAG3", "CTLA4", "PDCD1", "HAVCR2", "TIGIT")), nbin = 10, name = "exhaustion")
relapse_cd8_seurat <- AddModuleScore(relapse_cd8_seurat, features = list(cytotoxic_markers), nbin = 10, name = "cytotoxicity")

FeaturePlot(relapse_cd8_seurat, features = "inhibitory1", cols = c("lightgrey", "salmon"), order = T, min.cutoff = 0, max.cutoff = 0.2)
FeaturePlot(relapse_cd8_seurat, features = "exhaustion1", cols = c("lightgrey", "salmon"), order = T, min.cutoff = 0) #, max.cutoff = 0.2)
FeaturePlot(relapse_cd8_seurat, features = "cytotoxicity1", cols = c("lightgrey", "salmon"), order = T, min.cutoff = 0) #, max.cutoff = 0.2)
ggsave("results/relapse/umap_cd8_cytotoxicity.pdf", width = 5, height = 4)

relapse_cd8_seurat@meta.data %>% mutate(cluster = Idents(relapse_cd8_seurat)) %>%
  ggplot(aes(cluster, cytotoxicity1, fill = cluster)) + geom_violin(draw_quantiles = 0.5) + scale_fill_manual(values = getPalette4(7)) + labs(x = "") + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) +
  ggpubr::stat_compare_means()
ggsave("results/relapse/vln_cd8_cytotoxicity.pdf", width = 5, height = 4)
