
relapse_seurat <- readRDS("results/relapse_seurat.rds")

dir.create("results/manuscript/", showWarnings = F)
dir.create("results/manuscript/relapse/", showWarnings = F)

Idents(relapse_seurat) <- relapse_seurat$RNA_snn_res.0.5 %>% extractClusterNumber() %>% as.numeric %>% as.factor() %>% getClusterPhenotypesRelapse()
relapse_seurat$cluster <- Idents(relapse_seurat)

relapse_seurat@meta.data %>% ggplot(aes(cluster, nCount_RNA, fill = cluster)) + geom_violin() + scale_fill_manual(values = getPalette5(25)) + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) + labs(x = "")
ggsave("results/manuscript/relapse/vln_nCount_rna.pdf", width = 7, height = 4)

relapse_seurat@meta.data %>% ggplot(aes(cluster, nFeature_RNA, fill = cluster)) + geom_violin() + scale_fill_manual(values = getPalette5(25)) + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) + labs(x = "")
ggsave("results/manuscript/relapse/vln_nFeature_RNA.pdf", width = 7, height = 4)

relapse_seurat@meta.data %>% ggplot(aes(cluster, percent.mt, fill = cluster)) + geom_violin() + scale_fill_manual(values = getPalette5(25)) + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) + labs(x = "")
ggsave("results/manuscript/relapse/vln_percent.mt.pdf", width = 7, height = 4)

relapse_seurat@meta.data %>% ggplot(aes(cluster, percent.ribo, fill = cluster)) + geom_violin() + scale_fill_manual(values = getPalette5(25)) + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) + labs(x = "")
ggsave("results/manuscript/relapse/vln_percent.ribo.pdf", width = 7, height = 4)

## From public
df <- relapse_seurat@meta.data %>% left_join(dplyr::select(public_seurat@meta.data, barcode, cluster, singler_blueprint_pred, singler_hpca_pred) %>% dplyr::rename(public_cluster = cluster), by = "barcode")
relapse_seurat$public_cluster <- df$public_cluster
relapse_seurat$singler_blueprint_pred <- df$singler_blueprint_pred
relapse_seurat$singler_hpca_pred <- df$singler_hpca_pred

p <- DimPlot(relapse_seurat, cols = getPalette4(25), group.by = "public_cluster", label = T, repel = T) + theme_bw(base_size = 12) + theme(legend.position = "none") + labs(x = "UMAP1", y = "UMAP2")
ggsave(plot = p, "results/manuscript/relapse/latent_umap_public.png", width = 7, height = 6)


## Plot most notable markers
p <- DimPlot(relapse_seurat, cols = getPalette5(25), label = T, repel = T) + theme_bw(base_size = 12) + theme(legend.position = "none") + labs(x = "UMAP1", y = "UMAP2")
ggsave(plot = p, "results/manuscript/relapse/latent_umap.png", width = 7, height = 6)

p <- DotPlot(relapse_seurat, features = rev(unique(big_markers)), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y = "") + theme(legend.position = "top")
ggsave(plot = p, "results/manuscript/relapse/dotplot_big_markers.png", width = 11, height = 8)

p <- DotPlot(relapse_seurat, features = rev(unique(c("CD3E", "NCAM1", "GZMK", "SELL", "XCL1", "XCL2", "KLRC1", "IL7R", "LTB", "FCGR3A", "GZMA", "GZMB", "GZMH", "GZMM", "KLRC2", "ZEB2", "KLF2", "PRDM1", "GZMH", "LAG3"))), cols = "RdYlBu") +
  ggpubr::rotate_x_text(angle = 90) + labs(x = "", y  = "")
ggsave(plot = p, "results/manuscript/relapse/dotplot_dufva_markers.png", width = 9, height = 6)

guo_genes <- guo_markers %>% do.call(what = "c") %>% as.character() %>% rev %>% unique()
p <- DotPlot(relapse_seurat, features = guo_genes, cols = "RdYlBu") +
  ggpubr::rotate_x_text(angle = 90) + labs(x = "", y  = "")
ggsave(plot = p, "results/manuscript/relapse/dotplot_guo_markers.png", width = 12, height = 6)

p <- relapse_seurat@meta.data %>% group_by(cluster, singler_hpca_pred) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% filter(prop > 0.1) %>%
  ggplot(aes(singler_hpca_pred, prop, fill = singler_hpca_pred, label = singler_hpca_pred)) + geom_bar(stat = "identity") + add_guide + scale_fill_manual(values = getPalette(20)) + facet_wrap(~cluster, ncol = 6) + coord_flip() +
  theme_bw(base_size = 12) + labs(x = "") + theme(legend.position = "none") + geom_hline(yintercept = 0.5, linetype = "dotted") + ggrepel::geom_text_repel()
ggsave(plot = p, "results/manuscript/relapse/bar_predictions_hpca.png", width = 12, height = 8)

p <- relapse_seurat@meta.data %>% group_by(cluster, singler_blueprint_pred) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% filter(prop > 0.1) %>%
  ggplot(aes(singler_blueprint_pred, prop, fill = singler_blueprint_pred, label = singler_blueprint_pred)) + geom_bar(stat = "identity") + add_guide + scale_fill_manual(values = getPalette(18)) + facet_wrap(~cluster, ncol = 6) + coord_flip() +
  theme_bw(base_size = 12) + labs(x = "") + theme(legend.position = "none") + geom_hline(yintercept = 0.5, linetype = "dotted") + ggrepel::geom_text_repel() + ylim(values = c(0,1))
ggsave(plot = p, "results/manuscript/relapse/bar_predictions_blueprint.png", width = 12, height = 8)



clusters <- relapse_seurat@meta.data %>% mutate(clusters = Idents(relapse_seurat)) %>% group_by(clusters, !is.na(new_clonotypes_id)) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% filter(`!is.na(new_clonotypes_id)` == T) %>% filter(prop > 0.4) %>% pull(clusters) %>% as.character() %>% unique()

relapse_seurat@meta.data %>%
  mutate(cluster = Idents(relapse_seurat)) %>%
  # filter(cluster %in% clusters) %>%
  filter(!is.na(new_clonotypes_id)) %>%
  group_by(cluster, new_clonotypes_id) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  summarise(diversity = vegan::diversity(prop), gini = ineq::Gini(prop)) %>%

  ggplot(aes(gini,diversity,label=cluster)) + geom_point(shape = 21, fill = "lightgrey", size = 3) + ggrepel::geom_text_repel() + labs(x = "Gini index", y = "Shannon diveristy")
ggsave("results/manuscript/relapse/scatter_gini_clonality.pdf", width = 5, height = 4)




relapse_seurat@meta.data %>% mutate(clusters = Idents(relapse_seurat)) %>%
  group_by(clusters, !is.na(new_clonotypes_id)) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  filter(`!is.na(new_clonotypes_id)` == T) %>%

  ggplot(aes(reorder(clusters, prop), prop)) + geom_bar(stat = "identity") + ggpubr::rotate_x_text(45) +
  geom_hline(yintercept = 0.5, linetype = "dotted") + labs(x = "", y = "proportion of cells with TCR") + ylim(c(0,1))
ggsave("results/manuscript/relapse/bar_tcrab_cluster.pdf", width = 5, height = 4)

relapse_seurat@meta.data %>% mutate(clusters = Idents(relapse_seurat)) %>%
  group_by(clusters) %>% summarise(n = n()) %>%

  ggplot(aes(reorder(clusters, n), n, fill = clusters, label = n)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_manual(values = getPalette3(25))+ labs(x = "", y = "nCells") +
  geom_text() + theme_classic(base_size = 12)  + theme(legend.position = "none")
ggsave("results/manuscript/relapse/bar_n_cluster.pdf", width = 5, height = 4)




## Find DE-genes between all the clusters
all_markers <- FindAllMarkers(relapse_seurat, test.use = "t", verbose = T, max.cells.per.ident = 1e3)
fwrite(all_markers, "results/manuscript/relapse/all_markers_1e3.txt", sep = "\t", quote = F, row.names = F)
all_markers <- fread("results/manuscript/relapse/all_markers_1e3.txt") %>% mutate(cluster = cluster %>% extractClusterNumber() %>% as.numeric %>% as.factor() %>% getClusterPhenotypesRelapse())


set.seed(125)
top10  <- all_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

cells <- relapse_seurat@meta.data %>% mutate(cluster = Idents(relapse_seurat)) %>% group_by(cluster) %>% sample_n(29) %>% pull(barcode)
clonality_genes <- getClonalityGenes(relapse_seurat)
unwanted_genes  <- getUnwantedGenes(relapse_seurat)
relapse_seurat2  <- preprocessSeurat(relapse_seurat, cells.to.use = cells)

p <- DoHeatmap(relapse_seurat2, cells = cells, features = top10$gene, angle = 90, group.colors = getPalette(24)) + theme(legend.position = "none")
ggsave(plot = p, "results/manuscript/relapse/heatmap_top10_markers.png", width = 16, height = 20)

p <- DotPlot(relapse_seurat, features = rev(unique(top10$gene)), cols = "RdYlBu") + theme(legend.position = "none") + coord_flip() + ggpubr::rotate_x_text(45) + labs(x = "", y = "")
ggsave(plot = p, "results/manuscript/relapse/heatmap_top10_markers.png", width = 8, height = 27)









### Find relapse specific clusters
relapse_seurat$relapse[is.na(relapse_seurat$relapse)] <- "Diagnosis"
relapse_seurat$cluster <- Idents(relapse_seurat)

patient_df <- relapse_seurat@meta.data  %>% group_by(orig.ident, relapse) %>% summarise(n = n()) %>% dplyr::select(-n)
df <- relapse_seurat@meta.data  %>% group_by(orig.ident, cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% left_join(patient_df)

p.df <- lapply(unique(df$cluster), FUN = function(x){
  message(x)
  y <- df %>% filter(cluster == x)
  if(length(unique(y$relapse)) > 2){
    kruskal.test(prop~relapse, data = y) %>% broom::tidy() %>% mutate(cluster = x)
  }
}) %>% rbindlist() %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% arrange(p.adj)

fwrite(p.df, "results/manuscript/relapse/cml_kruskal_p_df.txt", sep = "\t", quote = F, row.names = F)

p.df %>% filter(p.adj < 0.05)


df %>% ggplot(aes(relapse, prop, fill = relapse)) +
  scale_fill_manual(values = getPalette2(7)) + labs(x = "") + theme_classic(base_size = 12) +
  geom_boxplot(outlier.shape = NA) + ggpubr::stat_compare_means(label = "p.format") + facet_wrap(~cluster, scales = "free_y") + ggpubr::rotate_x_text(45) + theme(legend.position = "none") + geom_jitter(size = 0.5)+ facets_nice
ggsave("results/manuscript/relapse/box_clusters_patient.pdf", width = 10, height = 7)


#### CMV specific clusters
relapse_seurat$cmv_status <- ifelse(relapse_seurat$orig.ident %in% c("706_dg", "716_dg", "720_dg",
                                                                   "Batch 2 baseline", "Batch 2 6m", "Batch 2 relapse",
                                                                   "Batch 3 baseline", "Batch 3 6m", "Batch 3 relapse",
                                                                   "Batch 4 baseline", "Batch 4 12m",
                                                                   "Batch 6 baseline", "Batch 6 relapse"
), "CMV pos", "CMV unknown")
relapse_seurat$cmv_status <- ifelse(relapse_seurat$orig.ident %in% c("Batch 5 baseline", "Batch 5 6m", "Batch 5 12m",
                                                                   "Batch 7 baseline", "Batch 7 relapse"), "CMV neg", relapse_seurat$cmv_status)
relapse_seurat$cmv_status <- factor(as.character(relapse_seurat$cmv_status), levels = c("CMV pos", "CMV neg", "CMV unknown"))

patient_df <- relapse_seurat@meta.data %>% filter(cluster != "5 CLL cells") %>% group_by(orig.ident, cmv_status) %>% summarise(n = n()) %>% dplyr::select(-n)
df <- relapse_seurat@meta.data %>% filter(cluster != "5 CLL cells") %>% group_by(orig.ident, cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% left_join(patient_df)

p.df <- lapply(unique(df$cluster), FUN = function(x){
  message(x)
  y <- df %>% filter(cluster == x)
  if(length(unique(y$cmv_status)) > 2){
    kruskal.test(prop~cmv_status, data = y) %>% broom::tidy() %>% mutate(cluster = x)
  }
}) %>% rbindlist() %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% arrange(p.adj)

fwrite(p.df, "results/manuscript/relapse/cml_cmv_p_df.txt", sep = "\t", quote = F, row.names = F)



df %>% ggplot(aes(cmv_status, prop, fill = cmv_status)) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "") + theme_classic(base_size = 12) +
  geom_boxplot(outlier.shape = NA) + ggpubr::stat_compare_means(label = "p.format") + facet_wrap(~cluster, scales = "free_y") + ggpubr::rotate_x_text(45) + theme(legend.position = "none") + geom_jitter(size = 0.5)+ facets_nice
ggsave("results/manuscript/relapse/box_clusters_cmv_patient.pdf", width = 10, height = 7)

df %>%
  filter(cmv_status != "CMV unknown") %>%
  ggplot(aes(cmv_status, prop, fill = cmv_status)) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "") + theme_classic(base_size = 12) +
  geom_boxplot(outlier.shape = NA) + ggpubr::stat_compare_means(label = "p.format") + facet_wrap(~cluster, scales = "free_y") + ggpubr::rotate_x_text(45) + theme(legend.position = "none") + geom_jitter(size = 0.5)+ facets_nice
ggsave("results/manuscript/relapse/box_clusters_cmv__knownpatient.pdf", width = 10, height = 7)





df <- relapse_seurat@meta.data %>% group_by(orig.ident,relapse,cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n))
median_df     <- df %>% group_by(relapse, cluster) %>% summarise(median = median(prop)) %>% group_by(cluster) %>% top_n(n = 1, wt = median) %>% arrange(desc(median))
umap_means_df <- relapse_seurat@reductions$latent_umap@cell.embeddings %>% as.data.frame %>% bind_cols(relapse_seurat@meta.data) %>% group_by(cluster) %>% summarise(umap1 = median(latent_umap_1), umap2 = median(latent_umap_2)) %>%
  left_join(median_df)

relapse_seurat@reductions$latent_umap@cell.embeddings %>% as.data.frame %>% bind_cols(relapse_seurat@meta.data) %>%
  ggplot(aes(latent_umap_1, latent_umap_2, color = relapse)) + stat_density_2d(aes(fill = ..level..), geom = "polygon") + facet_wrap(~relapse, ncol = 4) + scale_fill_distiller(direction=1) + theme(legend.position = "none") + theme_bw(base_size = 12) + labs(x = "UMAP 1", y = "UMAP 2") + facets_nice +
  ggrepel::geom_text_repel(data = umap_means_df, aes(umap1,umap2,label=cluster), color = "black", size = 3.5) +
  scale_color_manual(values = getPalette3(4)) + guides(color=FALSE) + theme(legend.position = "none")
ggsave("results/manuscript/relapse//umap_dens.png", width = 10, height = 4)







## Find project specific clusters
q <- NULL

table(relapse_seurat$relapse)
min_cells = 5000
for(i in 1:100){
  set.seed(i)
  q[[i]] <- relapse_seurat@meta.data %>% group_by(relapse) %>% sample_n(min_cells) %>% group_by(project) %>%
    group_by(cluster, relapse) %>% summarise(n = n()) %>% mutate(prop = n / sum(n))
}

q <- q %>% rbindlist()

df <- q %>% group_by(cluster,relapse) %>% summarise(median = median(prop)) %>% mutate(median = median/sum(median))
order_df <- df %>% filter(relapse == "Slow") %>% arrange(desc(median)) %>% dplyr::rename(tot = median)
df <- df %>% left_join(order_df)
df$tot[is.na(df$tot)] <- 0
df %>% ggplot(aes(reorder(cluster, tot),median,fill=relapse)) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette3(6)) + ggpubr::rotate_x_text(angle = 45) + labs(x = "", y = "prop")
ggsave("results/manuscript/relapse/bar_cluster_freq.pdf", width = 6, height = 4)





relapse_seurat@meta.data %>% group_by(orig.ident) %>% sample_n(3e3) %>%
  group_by(orig.ident, cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = "", prop, fill = cluster),color="black") + geom_bar(stat = "identity") + facet_wrap(~orig.ident, ncol = 4) + coord_polar(theta="y") + scale_fill_manual(values = getPalette(20)) + labs(x = "") + theme_void() + theme(legend.position = "bottom")
ggsave("results/manuscript/relapse/pie_new_cluster_project_resampled2.pdf", width = 7, height = 4)


relapse_seurat@meta.data %>% group_by(orig.ident) %>% sample_n(3e3) %>%
  group_by(orig.ident, cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = "", prop, fill = cluster),color="black") + geom_bar(stat = "identity") + facet_wrap(~orig.ident, ncol = 4) + coord_polar(theta="y") + scale_fill_manual(values = getPalette(20)) + labs(x = "") + theme_void() + theme(legend.position = "bottom")
ggsave("results/manuscript/relapse/pie_new_cluster_project_resampled_legend.pdf", width = 14, height = 8)






## DEGs
getDEGbyCluster <- function(seurat_object, cluster){


  ## If under 50 cells to begin with
  if(table(Idents(seurat_object)) %>% as.data.frame() %>% filter(Var1 == cluster) %>% pull(Freq) <= 50) return(NULL)
  message(paste0("===== ", cluster, " ====="))

  ## Subet to only cluster
  seurat_cluster         <- subset(seurat_object, ident = cluster)
  Idents(seurat_cluster) <- seurat_cluster$timepoint

  ## Calculate DEG only if there's at least 5 cells per time point
  n_df <- table(Idents(seurat_cluster)) %>% as.data.frame()

  n1 <- n_df %>% filter(Var1 == "baseline") %>% pull(Freq) >= 10
  n2 <- n_df %>% filter(Var1 == "6m") %>% pull(Freq) >= 10
  n3 <- n_df %>% filter(Var1 == "12m") %>% pull(Freq) >= 10

  if(length(n1) == 0) n1 <- FALSE
  if(length(n2) == 0) n2 <- FALSE
  if(length(n3) == 0) n3 <- FALSE

  cluster_markers_2v1 <- NULL
  cluster_markers_3v1 <- NULL
  cluster_markers_3v2 <- NULL

  if(n1 & n2) cluster_markers_2v1 <- FindMarkers(object = seurat_cluster, ident.1 = "baseline", ident.2 = "6m", only.pos = F, min.pct = 0.05, logfc.threshold = 0.1, return.thresh = 1, do.print = T, test.use = "t") %>% add_rownames(var = "gene") %>% mutate(timepoint = "2v1") %>% mutate(p_val_adj = p.adjust(p_val, method = "BH"))
  if(n3 & n1) cluster_markers_3v1 <- FindMarkers(object = seurat_cluster, ident.1 = "baseline", ident.2 = "12m", only.pos = F, min.pct = 0.05, logfc.threshold = 0.1, return.thresh = 1, do.print = T, test.use = "t") %>% add_rownames(var = "gene") %>% mutate(timepoint = "3v1") %>% mutate(p_val_adj = p.adjust(p_val, method = "BH"))
  if(n3 & n2) cluster_markers_3v2 <- FindMarkers(object = seurat_cluster, ident.1 = "6m", ident.2 = "12m", only.pos = F, min.pct = 0.05, logfc.threshold = 0.1, return.thresh = 1, do.print = T, test.use = "t") %>% add_rownames(var = "gene") %>% mutate(timepoint = "3v2") %>% mutate(p_val_adj = p.adjust(p_val, method = "BH"))

  df <- rbind(cluster_markers_2v1, cluster_markers_3v1, cluster_markers_3v2)

  if(!is.null(df)) df <- df %>% filter(p_val_adj < 0.05) %>% mutate(cluster = cluster, direction = ifelse(avg_logFC > 0, "up", "down"))
  return(df)

}

## Get the DEGs
relapse_seurat$timepoint <- relapse_seurat$relapse
relapse_seurat$timepoint <- plyr::revalue(relapse_seurat$timepoint, replace = c("Slow" = "12m", "Fast" = "6m", "Diagnosis" = "baseline"))
DEG_rp_cluster           <- lapply(levels(Idents(relapse_seurat)), getDEGbyCluster, seurat_object = relapse_seurat)
DEG_rp_cluster_df        <- DEG_rp_cluster %>% rbindlist()
DEG_rp_cluster_df$timepoint <- plyr::revalue(as.factor(DEG_rp_cluster_df$timepoint), replace = c("2v1" = "FastvDiagnosis", "3v1" = "SlowvDiagnosis", "3v2" = "SlowvFast"))
write.table(DEG_rp_cluster_df, "results/manuscript/relapse/deg_tki.txt", sep = "\t", quote = F, row.names = F)


DEG_rp_cluster_df %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% summarise(n = n()) %>%
  mutate(celltype = extractCoarsePhenotype(cluster)) %>%
  ggplot(aes(reorder(cluster,n), n)) + geom_bar(fill = "lightgrey", stat = "identity") + coord_flip() + labs(x = "", y = "DEG") + scale_fill_manual(values = getPalette3(9))
ggsave("results/manuscript/relapse/bar_deg.pdf", width = 5, height = 4)


volcano_df <- DEG_rp_cluster_df %>% filter(cluster == "10 CD4 treg") %>% filter(p_val_adj < 0.05)

ggplot(volcano_df, aes(avg_logFC, -log10(p_val), color = ifelse(avg_logFC > 0, "yes", "no"))) + geom_point(alpha=0.3) + facet_wrap(~timepoint) + xlim(values = c(-5,5)) +
  ggrepel::geom_text_repel(data = subset(volcano_df, gene %in% inhibitory_long), aes(label=gene), nudge_y = 10, font = 3) + theme(legend.position = "none") + scale_color_manual(values = c("dodgerblue", "salmon")) + facets_nice
ggsave("results/manuscript/relapse/volcano_treg.pdf", width = 9, height = 4)





### NK subsets
diet_nk_seurat$timepoint == "relapse"
table(diet_nk_seurat$project)

cells.to.keep <- diet_nk_seurat@meta.data %>% filter(timepoint == "relapse" | project == "CML dg") %>% pull(barcode)
relapse_nk_seurat <- subset(diet_nk_seurat, cells = cells.to.keep)
relapse_nk_seurat <- relapse_nk_seurat %>% getLatentUMAP() %>% fixSeurat()

DimPlot(relapse_nk_seurat, label = T, repel = T)

relapse_nk_seurat$relapse[is.na(relapse_nk_seurat$relapse)] <- "Diagnosis"

table(relapse_nk_seurat$relapse)



## Find project specific clusters
q <- NULL

table(relapse_nk_seurat$relapse)
min_cells = 500
for(i in 1:100){
  set.seed(i)
  q[[i]] <- relapse_nk_seurat@meta.data %>% group_by(relapse) %>% sample_n(min_cells) %>% group_by(project) %>%
    group_by(cluster, relapse) %>% summarise(n = n()) %>% mutate(prop = n / sum(n))
}

q <- q %>% rbindlist()

df <- q %>% group_by(cluster,relapse) %>% summarise(median = median(prop)) %>% mutate(median = median/sum(median))
order_df <- df %>% filter(relapse == "Slow") %>% arrange(desc(median)) %>% dplyr::rename(tot = median)
df <- df %>% left_join(order_df)
df$tot[is.na(df$tot)] <- 0
df %>% ggplot(aes(reorder(cluster, tot),median,fill=relapse)) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette3(6)) + ggpubr::rotate_x_text(angle = 45) + labs(x = "", y = "prop")
ggsave("results/manuscript/relapse/bar_nk_cluster_freq.pdf", width = 6, height = 4)
