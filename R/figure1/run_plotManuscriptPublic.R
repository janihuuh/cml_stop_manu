
dir.create("results/figure1/", showWarnings = F)

Idents(public_seurat) <- public_seurat$RNA_snn_res.0.5 %>% extractClusterNumber() %>% as.numeric %>% as.factor() %>% getClusterPhenotypesPublic()
public_seurat$cluster <- Idents(public_seurat)

public_seurat@meta.data %>% ggplot(aes(cluster, nCount_RNA, fill = cluster)) + geom_violin() + scale_fill_manual(values = getPalette5(25)) + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) + labs(x = "")
ggsave("results/figure1/vln_nCount_rna.pdf", width = 7, height = 4)

public_seurat@meta.data %>% ggplot(aes(cluster, nFeature_RNA, fill = cluster)) + geom_violin() + scale_fill_manual(values = getPalette5(25)) + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) + labs(x = "")
ggsave("results/figure1/vln_nFeature_RNA.pdf", width = 7, height = 4)

public_seurat@meta.data %>% ggplot(aes(cluster, percent.mt, fill = cluster)) + geom_violin() + scale_fill_manual(values = getPalette5(25)) + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) + labs(x = "")
ggsave("results/figure1/vln_percent.mt.pdf", width = 7, height = 4)

public_seurat@meta.data %>% ggplot(aes(cluster, percent.ribo, fill = cluster)) + geom_violin() + scale_fill_manual(values = getPalette5(25)) + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) + labs(x = "")
ggsave("results/figure1/vln_percent.ribo.pdf", width = 7, height = 4)


## Plot most notable markers
big_markers           <- c("CD3E", "TRAC",            ## T cell
                           "SELL", "IL7R", "LEF1", "TCF7", "CD4", "IL2",
                           "CD8A", "CD8B", "GZMA", "PRF1", "GZMB", "GNLY", "CXCR3",     ## CD8+ cell
                           "IFNG",
                           "GZMK", "CD27", "CD28", "IL2RA", "CCL5",             ## memory
                           "FOXP3", "PDCD1", "LAG3", "HAVCR2", "CTLA4", ## Inhibitory genes
                           "NKG7", "FCGR3A","KLRG1", "KLRB1", "KLRD1", "NCAM1", ## NK-cell
                           "LYZ", "CD14", "CST3",              ## monocytes
                           "FCER1A", "CLEC10A",                         ## cDC
                           "CLEC4C", "PTPRS", "TCF4",          ## pDC
                           "MS4A1", "CD19", "IL4R",                   ## b cells
                           "TNFRSF17", "JCHAIN",       ## plasma cells
                           "MKI67" )


DimPlot(public_seurat, label.size = 5, cols = getPalette5(25), label = T, repel = T) + theme_bw(base_size = 15) + theme(legend.position = "none") + labs(x = "UMAP1", y = "UMAP2")
ggsave("results/figure1/latent_umap.png", width = 9, height = 8)

p <- DotPlot(public_seurat, features = rev(unique(big_markers)), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y = "") + theme(legend.position = "top")
ggsave(plot = p, "results/figure1/dotplot_big_markers.png", width = 11, height = 8)

p <- DotPlot(public_seurat, features = rev(unique(c("CD3E", "NCAM1", "GZMK", "SELL", "XCL1", "XCL2", "KLRC1", "IL7R", "LTB", "FCGR3A", "GZMA", "GZMB", "GZMH", "GZMM", "KLRC2", "ZEB2", "KLF2", "PRDM1", "GZMH", "LAG3"))), cols = "RdYlBu") +
  ggpubr::rotate_x_text(angle = 90) + labs(x = "", y  = "")
ggsave(plot = p, "results/figure1/dotplot_dufva_markers.png", width = 9, height = 6)

cells.to.keep   <- public_seurat@meta.data %>% filter(grepl("Monocyte", cluster)) %>% pull(barcode)
monocyte_seurat <- subset(public_seurat, cells = cells.to.keep) %>% getLatentUMAP() %>% fixSeurat()

p <- DotPlot(public_seurat, features = rev(unique(do.call(zhang19_markers, what = "c"))), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y = "") + theme(legend.position = "top")
ggsave(plot = p, "results/figure1/dotplot_monocyte_markers.png", width = 11, height = 8)

guo_genes <- guo_markers %>% do.call(what = "c") %>% as.character() %>% rev %>% unique()
p <- DotPlot(public_seurat, features = guo_genes, cols = "RdYlBu") +
  ggpubr::rotate_x_text(angle = 90) + labs(x = "", y  = "")
ggsave(plot = p, "results/figure1/dotplot_guo_markers.png", width = 12, height = 6)

p <- public_seurat@meta.data %>% group_by(cluster, singler_hpca_pred) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% filter(prop > 0.1) %>%
  ggplot(aes(singler_hpca_pred, prop, fill = singler_hpca_pred, label = singler_hpca_pred)) + geom_bar(stat = "identity") + add_guide + scale_fill_manual(values = getPalette(20)) + facet_wrap(~cluster, ncol = 6) + coord_flip() +
  theme_bw(base_size = 12) + labs(x = "") + theme(legend.position = "none") + geom_hline(yintercept = 0.5, linetype = "dotted") + ggrepel::geom_text_repel()
ggsave(plot = p, "results/figure1/bar_predictions_hpca.png", width = 12, height = 8)

p <- public_seurat@meta.data %>% group_by(cluster, singler_blueprint_pred) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% filter(prop > 0.1) %>%
  ggplot(aes(singler_blueprint_pred, prop, fill = singler_blueprint_pred, label = singler_blueprint_pred)) + geom_bar(stat = "identity") + add_guide + scale_fill_manual(values = getPalette(18)) + facet_wrap(~cluster, ncol = 6) + coord_flip() +
  theme_bw(base_size = 12) + labs(x = "") + theme(legend.position = "none") + geom_hline(yintercept = 0.5, linetype = "dotted") + ggrepel::geom_text_repel() + ylim(values = c(0,1))
ggsave(plot = p, "results/figure1/bar_predictions_blueprint.png", width = 12, height = 8)





clusters <- public_seurat@meta.data %>% mutate(clusters = Idents(public_seurat)) %>% group_by(clusters, !is.na(new_clonotypes_id)) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% filter(`!is.na(new_clonotypes_id)` == T) %>% filter(prop > 0.4) %>% pull(clusters) %>% as.character() %>% unique()

public_seurat@meta.data %>%
  mutate(cluster = Idents(public_seurat)) %>%
  # filter(cluster %in% clusters) %>%
  filter(!is.na(new_clonotypes_id)) %>%
  group_by(cluster, new_clonotypes_id) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  summarise(diversity = vegan::diversity(prop), gini = ineq::Gini(prop)) %>%

  ggplot(aes(gini,diversity,label=cluster)) + geom_point(shape = 21, fill = "lightgrey", size = 3) + ggrepel::geom_text_repel() + labs(x = "Gini index", y = "Shannon diveristy")
ggsave("results/figure1/scatter_gini_clonality.pdf", width = 5, height = 4)




public_seurat@meta.data %>% mutate(clusters = Idents(public_seurat)) %>%
  group_by(clusters, !is.na(new_clonotypes_id)) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  filter(`!is.na(new_clonotypes_id)` == T) %>%

  ggplot(aes(reorder(clusters, prop), prop)) + geom_bar(stat = "identity") + ggpubr::rotate_x_text(45) +
  geom_hline(yintercept = 0.5, linetype = "dotted") + labs(x = "", y = "proportion of cells with TCR") + ylim(c(0,1))
ggsave("results/figure1/bar_tcrab_cluster.pdf", width = 5, height = 4)

public_seurat@meta.data %>% mutate(clusters = Idents(public_seurat)) %>%
  group_by(clusters) %>% summarise(n = n()) %>%

  ggplot(aes(reorder(clusters, n), n, fill = clusters, label = n)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_manual(values = getPalette3(25))+ labs(x = "", y = "nCells") +
  geom_text() + theme_classic(base_size = 12)  + theme(legend.position = "none")
ggsave("results/figure1/bar_n_cluster.pdf", width = 5, height = 4)




## Find DE-genes between all the clusters
all_markers <- FindAllMarkers(public_seurat, test.use = "t", verbose = T, max.cells.per.ident = 1e3)
fwrite(all_markers, "results/figure1/all_markers_1e3.txt", sep = "\t", quote = F, row.names = F)
all_markers <- fread("results/figure1/all_markers_1e3.txt") %>% mutate(cluster = cluster %>% extractClusterNumber() %>% as.numeric %>% as.factor() %>% getClusterPhenotypesPublic())


set.seed(125)
top10  <- all_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

cells <- public_seurat@meta.data %>% mutate(cluster = Idents(public_seurat)) %>% group_by(cluster) %>% sample_n(29) %>% pull(barcode)
p <- DoHeatmap(public_seurat2, cells = cells, features = top10$gene, angle = 90, group.colors = getPalette(24)) + theme(legend.position = "none")
ggsave(plot = p, "results/figure1/heatmap_top10_markers.png", width = 16, height = 20)

p <- DotPlot(public_seurat, features = rev(unique(top10$gene)), cols = "RdYlBu") + theme(legend.position = "none") + coord_flip() + ggpubr::rotate_x_text(45) + labs(x = "", y = "")
ggsave(plot = p, "results/figure1/heatmap_top10_markers.png", width = 8, height = 27)




clonality_genes <- getClonalityGenes(public_seurat)
unwanted_genes  <- getUnwantedGenes(public_seurat)
public_seurat2  <- preprocessSeurat(public_seurat, cells.to.use = colnames(public_seurat))





### Find project specific clusters
patient_df <- public_seurat@meta.data %>% filter(cluster != "5 CLL cells") %>% group_by(orig.ident, project) %>% summarise(n = n()) %>% dplyr::select(-n)
df <- public_seurat@meta.data %>% filter(cluster != "5 CLL cells") %>% group_by(orig.ident, cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% left_join(patient_df)

p.df <- lapply(unique(df$cluster), FUN = function(x){
  message(x)
  y <- df %>% filter(cluster == x)
  if(length(unique(y$project)) > 2){
    kruskal.test(prop~project, data = y) %>% broom::tidy() %>% mutate(cluster = x)
  }
}) %>% rbindlist() %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% arrange(p.adj)

fwrite(p.df, "results/figure1/cml_kruskal_p_df.txt", sep = "\t", quote = F, row.names = F)

p.df %>% filter(p.adj < 0.05)


df %>% ggplot(aes(project, prop, fill = project)) +
  scale_fill_manual(values = getPalette2(7)) + labs(x = "") + theme_classic(base_size = 12) +
  geom_boxplot(outlier.shape = NA) + ggpubr::stat_compare_means(label = "p.format") + facet_wrap(~cluster, scales = "free_y") + ggpubr::rotate_x_text(45) + theme(legend.position = "none") + geom_jitter(size = 0.5)+ facets_nice
ggsave("results/figure1/box_clusters_patient.pdf", width = 10, height = 7)

patient_df <- public_seurat@meta.data %>% filter(cluster != "5 CLL cells") %>% mutate(project = ifelse(project %in% c("CML on TKI", "CML dg", "CML off TKI"), "CML", "other")) %>% group_by(orig.ident, project) %>% summarise(n = n()) %>% dplyr::select(-n)
df <- public_seurat@meta.data %>% filter(cluster != "5 CLL cells") %>% group_by(orig.ident, cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% left_join(patient_df)

df %>%
  ggplot(aes(project, prop, fill = project)) +
  scale_fill_manual(values = getPalette4(4)) + labs(x = "") + theme_classic(base_size = 12) +
  geom_boxplot(outlier.shape = NA) + ggpubr::stat_compare_means(label = "p.format") + facet_wrap(~cluster, scales = "free_y") + ggpubr::rotate_x_text(45) + theme(legend.position = "none") + geom_jitter(size = 0.5) + facets_nice
ggsave("results/figure1/box_cml_clusters_patient.pdf", width = 10, height = 7)

p.df <- lapply(unique(df$cluster), FUN = function(x){
  message(x)
  y <- df %>% filter(cluster == x)
  if(length(unique(y$project)) == 2){
    wilcox.test(prop~project, data = y) %>% broom::tidy() %>% mutate(cluster = x)
  }
}) %>% rbindlist() %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% arrange(p.adj)

fwrite(p.df, "results/figure1/cml_p_df.txt", sep = "\t", quote = F, row.names = F)

p.df %>% filter(p.adj < 0.01)

x@meta.data %>% group_by(cmv_status, patient) %>% summarise(n = n()) %>% dplyr::select(-n)

#### CMV specific clusters
public_seurat$cmv_status <- ifelse(public_seurat$orig.ident %in% c("706_dg", "716_dg", "720_dg",
                                                                   "Batch 2 baseline", "Batch 2 6m", "Batch 2 relapse",
                                                                   "Batch 3 baseline", "Batch 3 6m", "Batch 3 relapse",
                                                                   "Batch 4 baseline", "Batch 4 12m",
                                                                   "Batch 6 baseline", "Batch 6 relapse"
                                                                   ), "CMV pos", "CMV unknown")
public_seurat$cmv_status <- ifelse(public_seurat$orig.ident %in% c("Batch 5 baseline", "Batch 5 6m", "Batch 5 12m",
                                                                   "Batch 7 baseline", "Batch 7 relapse"), "CMV neg", public_seurat$cmv_status)
public_seurat$cmv_status <- factor(as.character(public_seurat$cmv_status), levels = c("CMV pos", "CMV neg", "CMV unknown"))

patient_df <- public_seurat@meta.data %>% filter(cluster != "5 CLL cells") %>% group_by(orig.ident, cmv_status) %>% summarise(n = n()) %>% dplyr::select(-n)
df <- public_seurat@meta.data %>% filter(cluster != "5 CLL cells") %>% group_by(orig.ident, cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% left_join(patient_df)

p.df <- lapply(unique(df$cluster), FUN = function(x){
  message(x)
  y <- df %>% filter(cluster == x)
  if(length(unique(y$cmv_status)) > 2){
    kruskal.test(prop~cmv_status, data = y) %>% broom::tidy() %>% mutate(cluster = x)
  }
}) %>% rbindlist() %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% arrange(p.adj)

fwrite(p.df, "results/figure1/cml_cmv_p_df.txt", sep = "\t", quote = F, row.names = F)



df %>% ggplot(aes(cmv_status, prop, fill = cmv_status)) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "") + theme_classic(base_size = 12) +
  geom_boxplot(outlier.shape = NA) + ggpubr::stat_compare_means(label = "p.format") + facet_wrap(~cluster, scales = "free_y") + ggpubr::rotate_x_text(45) + theme(legend.position = "none") + geom_jitter(size = 0.5)+ facets_nice
ggsave("results/figure1/box_clusters_cmv_patient.pdf", width = 10, height = 7)

df %>%
  filter(cmv_status != "CMV unknown") %>%
  ggplot(aes(cmv_status, prop, fill = cmv_status)) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "") + theme_classic(base_size = 12) +
  geom_boxplot(outlier.shape = NA) + ggpubr::stat_compare_means(label = "p.format") + facet_wrap(~cluster, scales = "free_y") + ggpubr::rotate_x_text(45) + theme(legend.position = "none") + geom_jitter(size = 0.5)+ facets_nice
ggsave("results/figure1/box_clusters_cmv_knownpatient.pdf", width = 10, height = 7)






median_df     <- df %>% group_by(project, cluster) %>% summarise(median = median(prop)) %>% group_by(cluster) %>% top_n(n = 1, wt = median) %>% arrange(desc(median))
umap_means_df <- public_seurat@reductions$latent_umap@cell.embeddings %>% as.data.frame %>% bind_cols(public_seurat@meta.data) %>% group_by(cluster) %>% summarise(umap1 = median(latent_umap_1), umap2 = median(latent_umap_2)) %>%
  left_join(median_df)

public_seurat@reductions$latent_umap@cell.embeddings %>% as.data.frame %>% bind_cols(public_seurat@meta.data) %>%
  ggplot(aes(latent_umap_1, latent_umap_2, color = project)) + stat_density_2d(aes(fill = ..level..), geom = "polygon") + facet_wrap(~project, ncol = 3) + scale_fill_distiller(direction=1) + theme(legend.position = "none") + theme_bw(base_size = 12) + labs(x = "UMAP 1", y = "UMAP 2") + facets_nice +
  ggrepel::geom_text_repel(data = umap_means_df, aes(umap1,umap2,label=cluster), color = "black", size = 3.5) +
  scale_color_manual(values = getPalette3(7)) + guides(color=FALSE) + theme(legend.position = "none")
ggsave("results/figure1//umap_dens.png", width = 7.5, height = 4)








### Study NK cells
diet_nk_seurat <- readRDS("results/diet_nk_seurat.rds")

diet_nk_seurat$project <- "none"
diet_nk_seurat$project <- ifelse(grepl(diet_nk_seurat$orig.ident, pattern = "^CLL"), "CLL", diet_nk_seurat$project)
diet_nk_seurat$project <- ifelse(grepl(diet_nk_seurat$orig.ident, pattern = "^Batch"), "CML off TKI", diet_nk_seurat$project)
diet_nk_seurat$project <- ifelse(grepl(diet_nk_seurat$orig.ident, pattern = "baseline"), "CML on TKI", diet_nk_seurat$project)
diet_nk_seurat$project <- ifelse(grepl(diet_nk_seurat$orig.ident, pattern = "^7"), "CML dg", diet_nk_seurat$project)
diet_nk_seurat$project <- ifelse(grepl(diet_nk_seurat$orig.ident, pattern = "^LB"), "NSCLC from blood", diet_nk_seurat$project)
diet_nk_seurat$project <- ifelse(grepl(diet_nk_seurat$orig.ident, pattern = "^RB"), "RCC from blood", diet_nk_seurat$project)
diet_nk_seurat$project <- ifelse(grepl(diet_nk_seurat$orig.ident, pattern = "^healthy"), "Healthy", diet_nk_seurat$project)

diet_nk_seurat$project <- factor(as.character(diet_nk_seurat$project), levels = c("CML dg", "CML on TKI", "CML off TKI", "Healthy", "CLL", "NSCLC from blood", "RCC from blood"))


p <- DimPlot(diet_nk_seurat, cols = getPalette5(8), label = T, repel = T) + theme_bw(base_size = 12) + theme(legend.position = "none") + labs(x = "UMAP1", y = "UMAP2")
ggsave(plot = p, "results/figure1/latent_nk_umap.png", width = 7, height = 6)



patient_df <- diet_nk_seurat@meta.data %>% group_by(orig.ident, project) %>% summarise(n = n()) %>% dplyr::select(-n)
df <- diet_nk_seurat@meta.data %>% group_by(orig.ident, cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% left_join(patient_df)
median_df     <- df %>% group_by(project, cluster) %>% summarise(median = median(prop)) %>% group_by(cluster) %>% top_n(n = 1, wt = median) %>% arrange(desc(median))
umap_means_df <- diet_nk_seurat@reductions$latent_umap@cell.embeddings %>% as.data.frame %>% bind_cols(diet_nk_seurat@meta.data) %>% group_by(cluster) %>% summarise(umap1 = median(latent_umap_1), umap2 = median(latent_umap_2)) %>%
  left_join(median_df)

diet_nk_seurat@reductions$latent_umap@cell.embeddings %>% as.data.frame %>% bind_cols(diet_nk_seurat@meta.data) %>%
  ggplot(aes(latent_umap_1, latent_umap_2, color = project)) + stat_density_2d(aes(fill = ..level..), geom = "polygon") + facet_wrap(~project, ncol = 4) + scale_fill_distiller(direction=1) + theme(legend.position = "none") + theme_bw(base_size = 12) + labs(x = "UMAP 1", y = "UMAP 2") + facets_nice +
  ggrepel::geom_text_repel(data = umap_means_df, aes(umap1,umap2,label=cluster), color = "black", size = 3.5) +
  scale_color_manual(values = getPalette3(7)) + guides(color=FALSE) + theme(legend.position = "none")
ggsave("results/figure1/umap_nk_dens.png", width = 7.5, height = 4)




patient_df <- diet_nk_seurat@meta.data %>% group_by(orig.ident, project) %>% summarise(n = n()) %>% dplyr::select(-n)
df <- diet_nk_seurat@meta.data %>% group_by(orig.ident, cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% left_join(patient_df)

df %>% ggplot(aes(project, prop, fill = project)) +
  scale_fill_manual(values = getPalette2(7)) + labs(x = "") + theme_classic(base_size = 12) +
  geom_boxplot(outlier.shape = NA) + ggpubr::stat_compare_means(label = "p.format") + facet_wrap(~cluster, scales = "free_y") + ggpubr::rotate_x_text(45) + theme(legend.position = "none") + geom_jitter(size = 0.5)+ facets_nice
ggsave("results/figure1/box_nk_clusters_patient.pdf", width = 10, height = 7)


patient_df <- diet_nk_seurat@meta.data %>% mutate(project = ifelse(project %in% c("CML on TKI", "CML dg", "CML off TKI"), "CML", "other")) %>% group_by(orig.ident, project) %>% summarise(n = n()) %>% dplyr::select(-n)
diet_nk_seurat@meta.data %>% group_by(orig.ident, cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% left_join(patient_df) %>%
  ggplot(aes(project, prop, fill = project)) +
  scale_fill_manual(values = getPalette4(4)) + labs(x = "") + theme_classic(base_size = 12) +
  geom_boxplot(outlier.shape = NA) + ggpubr::stat_compare_means(label = "p.format") + facet_wrap(~cluster, scales = "free_y") + ggpubr::rotate_x_text(45) + theme(legend.position = "none") + geom_jitter(size = 0.5) + facets_nice
ggsave("results/figure1/box_cml_nk_clusters_patient.pdf", width = 10, height = 7)



diet_nk_seurat@meta.data %>% mutate(cluster = Idents(diet_nk_seurat)) %>%
  group_by(project, cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  ggplot(aes(project, prop, fill = cluster)) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette5(8)) + labs(x = "") + ggpubr::rotate_x_text(45)
ggsave("results/figure1/bar_nk_clusters.pdf", width = 5, height = 5)



q <- NULL; min_cells <- 50

for(i in 1:100){
  set.seed(i)
  q[[i]] <- diet_nk_seurat@meta.data %>% mutate(cluster = Idents(diet_nk_seurat)) %>%
    group_by(project) %>% sample_n(min_cells) %>%
    group_by(cluster, project) %>% summarise(n = n()) %>% mutate(prop = n / sum(n))
}

q <- q %>% rbindlist()

tot_df <- q %>% group_by(cluster, project) %>% summarise(prop = median(prop)) %>% group_by(cluster) %>% summarise(tot = sum(prop))
q %>% group_by(cluster, project) %>% summarise(prop = median(prop)) %>% left_join(tot_df, by = "cluster") %>% mutate(prop = prop / tot) %>%
  ggplot(aes(cluster, prop, fill = project)) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette3(7)) + labs(x = "") + ggpubr::rotate_x_text(45)
ggsave("results/figure1/bar_nk_clusters2.pdf", width = 5, height = 4)





kir_genes <- grep("^KIR", rownames(diet_nk_seurat), value = T)
p <- FeaturePlot(diet_nk_seurat, features = kir_genes, ncol = 3, cols = c("lightgrey", "red"), order = T) + theme_bw(base_size = 12)
ggsave(plot = p, "results/figure1/umap_nk_kir.png", width = 9, height = 7)

diet_nk_seurat2 <- diet_nk_seurat
Idents(diet_nk_seurat2) <- diet_nk_seurat2$orig.ident
DotPlot(diet_nk_seurat2, features = kir_genes, cols = "RdYlBu") + labs(x = "", y = "") + coord_flip() + ggpubr::rotate_x_text(angle = 90) + theme(legend.position = "top")
ggsave("results/figure1/dotplot_nk_kir.pdf", width = 8, height = 4)












## Run slingshot
diet_nk_seurat <- readRDS("results/diet_nk_seurat.rds")

require(slingshot)
require(SingleCellExperiment)
require(SummarizedExperiment)

diet_nk_seurat$orig.clusters <- Idents(diet_nk_seurat)
diet_nk_sce <- as.SingleCellExperiment(diet_nk_seurat)

diet_nk_sling <- runSlingshot(diet_nk_sce, cluster_with_earliset_timepoint = "5 CD56 bright", reducedDim = "LATENT_UMAP")
diet_nk_sling <- readRDS("results/diet_nk_sling.rds")

slingshot_object=diet_nk_sling
reducedDim="LATENT_UMAP"

cololors                <- slingshot_object$orig.clusters %>% extractClusterNumber() %>% as.numeric()
cololors[cololors == 0] <- getPalette5(8)[1]
cololors[cololors == 1] <- getPalette5(8)[2]
cololors[cololors == 2] <- getPalette5(8)[3]
cololors[cololors == 3] <- getPalette5(8)[4]
cololors[cololors == 4] <- getPalette5(8)[5]
cololors[cololors == 5] <- getPalette5(8)[6]
cololors[cololors == 6] <- getPalette5(8)[7]
cololors[cololors == 7] <- getPalette5(8)[8]

curve_cols <- c("black", "darkred", "darkblue", "darkolivegreen4")
sling_curve <- SlingshotDataSet(slingshot_object)

pdf("results/public/umap_sling_nk.pdf", width = 5, height = 5)
plot(reducedDims(slingshot_object)[[reducedDim]], col = cololors, pch = 16, asp = 1)
for(i in 1:2){
  lines(sling_curve@curves[[i]], lwd = 4, col = curve_cols[i])
}
dev.off()


png("results/public/umap_sling_nk.png", width = 5, height = 5, units = "in", res = 1024)
plot(reducedDims(slingshot_object)[[reducedDim]], col = cololors, pch = 16, asp = 1)
for(i in 1:2){
  lines(sling_curve@curves[[i]], lwd = 4, col = curve_cols[i])
}
dev.off()





### Trajectories in different diseases
diet_nk_seurat$orig.clusters <- Idents(diet_nk_seurat)
diet_nk_sce <- as.SingleCellExperiment(diet_nk_seurat)

diet_nk_sling <- runSlingshot(diet_nk_sce, cluster_with_earliset_timepoint = "5 CD56 bright", reducedDim = "LATENT_UMAP")
diet_nk_sling <- readRDS("results/diet_nk_sling.rds")

pdf("results/public/umap_sling_nk.pdf", width = 5, height = 5)
plot(reducedDims(diet_nk_sling)[[reducedDim]], col = cololors, pch = 16, asp = 1)
for(i in 1:2){
  lines(sling_curve@curves[[i]], lwd = 4, col = curve_cols[i])
}
dev.off()




## In differet diseases
getSlingNk <- function(project_temp){

  cells.to.keep       <- diet_nk_seurat@meta.data %>% filter(project == project_temp) %>% pull(barcode)
  project_seurat      <- subset(diet_nk_seurat, cells = cells.to.keep) %>% as.SingleCellExperiment()
  proje_nk_sling      <- runSlingshot(project_seurat, cluster_with_earliset_timepoint = "5 CD56 bright", reducedDim = "LATENT_UMAP")
  return(proje_nk_sling)

}

cml_dg_sling      <- getSlingNk(project_temp = "CML dg")
cml_on_tki_sling  <- getSlingNk(project_temp = "CML on TKI")
cml_off_tki_sling <- getSlingNk(project_temp = "CML off TKI")
nsclc_sling       <- getSlingNk(project_temp = "NSCLC from blood")
rcc_sling         <- getSlingNk(project_temp = "RCC from blood")
cll_sling         <- getSlingNk(project_temp = "CLL")

cml_dg_curve      <- SlingshotDataSet(cml_dg_sling)
cml_on_tki_curve  <- SlingshotDataSet(cml_on_tki_sling)
cml_off_tki_curve <- SlingshotDataSet(cml_off_tki_sling)
nsclc_curve       <- SlingshotDataSet(nsclc_sling)
rcc_curve         <- SlingshotDataSet(rcc_sling)
cll_curve         <- SlingshotDataSet(cll_sling)

length(cml_dg_curve@curves)
length(cml_on_tki_curve@curves)
length(cml_off_tki_curve@curves)
length(nsclc_curve@curves)
length(rcc_curve@curves)
length(cll_curve@curves)

png("results/figure1/umap_sling_nk_disease.png", width = 1064*2/3, height = 1064*2/3*2/3, res = 80)

par(mfrow=c(2,3))

plot(reducedDims(diet_nk_sling)[[reducedDim]], col = cololors, pch = 16, asp = 1, main = "CML dg")
for(i in 1:2){lines(cml_dg_curve@curves[[i]], lwd = 4, col = curve_cols[1])}

plot(reducedDims(diet_nk_sling)[[reducedDim]], col = cololors, pch = 16, asp = 1, main = "CML on-TKI")
for(i in 1:2){lines(cml_on_tki_curve@curves[[i]], lwd = 4, col = curve_cols[2])}

plot(reducedDims(diet_nk_sling)[[reducedDim]], col = cololors, pch = 16, asp = 1, main = "CML off-TKI")
for(i in 1:2){lines(cml_off_tki_curve@curves[[i]], lwd = 4, col = curve_cols[3])}

plot(reducedDims(diet_nk_sling)[[reducedDim]], col = cololors, pch = 16, asp = 1, main = "NSCLC")
for(i in 1:1){lines(nsclc_curve@curves[[i]], lwd = 4, col = curve_cols[4])}

plot(reducedDims(diet_nk_sling)[[reducedDim]], col = cololors, pch = 16, asp = 1, main = "RCC")
for(i in c(1,3)){lines(rcc_curve@curves[[i]], lwd = 4, col = curve_cols[5])}

plot(reducedDims(diet_nk_sling)[[reducedDim]], col = cololors, pch = 16, asp = 1, main = "CLL")
for(i in 1:2){lines(cll_curve@curves[[i]], lwd = 4, col = curve_cols[6])}

dev.off()






curve_cols <- c("black", "darkred", "darkblue", "darkolivegreen4", "salmon", "dodgerblue")
pdf("results/public/umap_sling_nk.pdf", width = 5, height = 5)

plot(reducedDims(diet_nk_sling)[[reducedDim]], col = cololors, pch = 16, asp = 1)
for(i in 1){lines(cml_dg_curve@curves[[i]], lwd = 4, col = curve_cols[1])}
for(i in 1){lines(cml_on_tki_curve@curves[[i]], lwd = 4, col = curve_cols[2])}
for(i in 1){lines(cml_off_tki_curve@curves[[i]], lwd = 4, col = curve_cols[3])}
for(i in 1){lines(nsclc_curve@curves[[i]], lwd = 4, col = curve_cols[4])}
for(i in 1){lines(rcc_curve@curves[[i]], lwd = 4, col = curve_cols[5])}
for(i in 1){lines(cll_curve@curves[[i]], lwd = 4, col = curve_cols[6])}

dev.off()


