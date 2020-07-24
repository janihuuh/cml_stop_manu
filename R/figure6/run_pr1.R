

## Run PR1-specific
# pr1_tcrb <- fread("tcrb_data/pr1_total_unfiltered.txt")

pr1_tcrb <- lapply(list.files("tcrb_data/sorted/pr1_tolerable/", full.names = T), fread) %>% rbindlist() %>% filter(count > 1)
nlv_tcrb <- lapply(list.files("tcrb_data/sorted/", full.names = T, pattern = "NLV"), fread) %>% rbindlist() %>% filter(count > 1)

# pr1_tcrb2 <- fread("tcrb_data/pr1_total_filtered.txt")
cml_tcr_seurat$pr1_specific <- ifelse(cml_tcr_seurat$trb_cdr3s_aa %in% pr1_tcrb$cdr3aa, "pr1_specific", "other")
cml_tcr_seurat$nlv_specific <- ifelse(cml_tcr_seurat$trb_cdr3s_aa %in% nlv_tcrb$cdr3aa, "nlv_specific", "other")

cml_tcr_seurat@meta.data %>%
  group_by(orig.ident, pr1_specific) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  ggplot(aes(orig.ident,prop,fill=pr1_specific)) + geom_bar(stat = "identity") + scale_fill_manual(values = c("lightgrey", "salmon")) + labs(x = "") + theme_classic(base_size = 12) + ggpubr::rotate_x_text()
ggsave("results/tcr/bar_pr1.pdf", width = 5, height = 4)

cells.to.keep <- cml_tcr_seurat@meta.data %>% filter(pr1_specific == "pr1_specific" & patient %in% c("Batch 5", "Batch 7")) %>% pull(barcode)
pr1_seurat    <- subset(cml_tcr_seurat, cells = cells.to.keep)
pr1_seurat    <- pr1_seurat %>% getLatentUMAP() %>% fixSeurat()
pr1_seurat    <- pr1_seurat %>% getLatentClustering()

res                <- c(seq(0.1, 1, 0.1), seq(1.2, 2, 0.2), 2.5, 3)
clustering_columns <- grep("res", colnames(pr1_seurat@meta.data), value = T)
clustering_columns <- clustering_columns[order(substr(clustering_columns, 13, nchar(clustering_columns)) %>% as.numeric())]

q <- NULL
i <- 1

for(clustering_column in clustering_columns){

  message(clustering_column)
  q[[i]] <- pr1_seurat@meta.data[,clustering_column] %>% levels %>% length
  i <- i + 1

}

data.frame(resolution = res, nClusters = q) %>%
  ggplot(aes((resolution),nClusters), label=nClusters) + geom_point(shape = 21) + theme_bw()
ggsave("results/tcr/scatter_res_nCluster_pr1.png", width = 4, height = 3)

## Optim amount of clusters: res 0.3
DimPlot(pr1_seurat, group.by = "RNA_snn_res.0.3", label = T, repel = T) + theme_bw(base_size = 12) + theme(legend.position = "none") + labs(x = "UMAP 1", y = "UMAP 2") + facets_nice + scale_color_manual(values = getPalette3(4))

## Find markers
Idents(pr1_seurat) <- pr1_seurat$RNA_snn_res.0.3
pr1_markers <- FindAllMarkers(pr1_seurat, test.use = "t")
pr1_markers <- pr1_markers %>% mutate(dir = ifelse(avg_logFC > 0, "up", "down")) %>% filter(p_val_adj < 0.05)
write.table(pr1_markers, "results/tcr/deg_pr1.txt", sep = "\t", quote = F, row.names = F)

## Visualize known markers
pr1_markers %>% filter(dir == "up") %>% filter(gene %in% big_markers)
pr1_markers %>% filter(dir == "up") %>% filter(gene %in% inhibitory_long)
pr1_markers %>% filter(dir == "up") %>% filter(gene %in% big_markers)

Idents(pr1_seurat) <- pr1_seurat$RNA_snn_res.0.3 %>% getClusterPhenotypesPR1()

p <- DotPlot(pr1_seurat, features = rev(unique(big_markers)), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y = "")
ggsave(plot = p, "results/tcr/dotplot_big_markers_pr1.png", width = 14, height = 4)

guo_genes <- guo_markers %>% do.call(what = "c")
p <- DotPlot(pr1_seurat, features = rev(unique(guo_genes)), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y = "")
ggsave(plot = p, "results/tcr/dotplot_guo_markers_pr1.png", width = 14, height = 4)


## Visualize known scores
pr1_seurat <- AddModuleScore(pr1_seurat, features = list(unique(inhibitory_long)), nbin = 10, name = "inhibitory")
pr1_seurat <- AddModuleScore(pr1_seurat, features = list( c("LAG3", "CTLA4", "PDCD1", "HAVCR2", "TIGIT")), nbin = 10, name = "exhaustion")
pr1_seurat <- AddModuleScore(pr1_seurat, features = list(cytotoxic_markers), nbin = 10, name = "cytotoxicity")

FeaturePlot(pr1_seurat, label = F, features = "inhibitory1", cols = c("lightgrey", "salmon"), order = T, min.cutoff = 0)
ggsave("results/tcr/umap_pr1_inhibitory.png", width = 5, height = 4)

FeaturePlot(pr1_seurat, label = F, features = "exhaustion1", cols = c("lightgrey", "salmon"), order = T, min.cutoff = 0, max.cutoff = 1) #, max.cutoff = 0.2)
ggsave("results/tcr/umap_pr1_exhaustion.png", width = 5, height = 4)

FeaturePlot(pr1_seurat, features = "cytotoxicity1", cols = c("lightgrey", "salmon"), order = T, min.cutoff = 0) #, max.cutoff = 0.2)
ggsave("results/tcr/umap_pr1_cytotoxicity.png", width = 5, height = 4)

pr1_seurat@meta.data %>% mutate(cluster = Idents(pr1_seurat)) %>%
  ggplot(aes(cluster, cytotoxicity1, fill = cluster)) + geom_violin(draw_quantiles = 0.5) + scale_fill_manual(values = getPalette3(4)) + labs(x = "") + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) +
  ggpubr::stat_compare_means()
ggsave("results/tcr/vln_pr1_cytotoxicity.pdf", width = 5, height = 4)

pr1_seurat@meta.data %>% mutate(cluster = Idents(pr1_seurat)) %>%
  ggplot(aes(cluster, inhibitory1, fill = cluster)) + geom_violin(draw_quantiles = 0.5) + scale_fill_manual(values = getPalette3(4)) + labs(x = "") + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) +
  ggpubr::stat_compare_means()
ggsave("results/tcr/vln_pr1_inhibitory.pdf", width = 5, height = 4)

pr1_seurat@meta.data %>% mutate(cluster = Idents(pr1_seurat)) %>%
  ggplot(aes(cluster, exhaustion1, fill = cluster)) + geom_violin(draw_quantiles = 0.5) + scale_fill_manual(values = getPalette3(4)) + labs(x = "") + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) +
  ggpubr::stat_compare_means()
ggsave("results/tcr/vln_pr1_exhaustion.pdf", width = 5, height = 4)




getClusterPhenotypesPR1 <- function(clusters){

  plyr::revalue(clusters, replace = c("0"	 = "0 CD8+ effector/exhausted",
                                      "1"	 = "1 CD8+ EM",
                                      "2"	 = "2 CD8+ IFNg"))
}




Idents(pr1_seurat) %>% table()
DimPlot(pr1_seurat, label = F, repel = T)
DimPlot(pr1_seurat, label = F, repel = T, group.by = "patient") + labs(color = "Patient") + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(4))
ggsave("results/tcr/umap_pr1.png", width = 5, height = 4)

DimPlot(pr1_seurat, label = T, repel = T) + theme_bw(base_size = 12) + scale_color_manual(values = getPalette3(4)) + theme(legend.position = "none")
ggsave("results/tcr/umap_pr1_cluster_new.png", width = 5, height = 4)

patient_df <- pr1_seurat@meta.data %>% group_by(orig.ident, patient, timepoint) %>% summarise(n = n()) %>% dplyr::select(-n)

pr1_seurat@meta.data %>% mutate(cluster = Idents(pr1_seurat)) %>%
  group_by(patient, timepoint, cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  ggplot(aes(timepoint, prop, group = cluster, color = cluster)) + geom_point() + geom_path() + facet_wrap(~patient) + scale_color_manual(values = getPalette3(4)) + theme_classic(base_size = 12) + labs(x = "") + ggpubr::rotate_x_text(45)
ggsave("results/tcr/path_pr1_tcr.pdf", width = 7, height = 4)

library(treemapify)
df <- pr1_seurat@meta.data %>% dplyr::select(patient, new_clonotypes_id, orig.ident, timepoint) %>% group_by(orig.ident, new_clonotypes_id) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% left_join(patient_df)
df <- df %>% ungroup() %>% mutate(orig.ident = factor(as.character(orig.ident), levels = c("Batch 5 baseline", "Batch 5 6m", "Batch 5 12m", "Batch 7 baseline", "Batch 7 relapse")))
ggplot(df, aes(area = prop, fill = new_clonotypes_id)) + geom_treemap() + facet_wrap(~orig.ident) + theme_bw(base_size = 17) + theme(legend.position = "none") + scale_fill_manual(values = getPalette(408)) + facets_nice
ggsave("results/tcr/treemap_pr1.pdf", width = 7, height = 5)

total_df <- cml_tcr_seurat@meta.data %>% group_by(orig.ident) %>% summarise(n = n())

pr1_seurat@meta.data %>% mutate(cluster = Idents(pr1_seurat)) %>%
  group_by(orig.ident) %>% summarise(n = n()) %>% left_join(total_df, by = "orig.ident") %>% left_join(patient_df, by = "orig.ident") %>%
  mutate(prop = n.x / n.y) %>%
  ggplot(aes(timepoint, prop, group = patient, color = patient, fill = patient)) + geom_path() +  geom_point(shape = 21, size = 3, color = "black") + theme_classic(base_size = 12) + labs(x = "") +
  scale_fill_manual(values = getPalette3(4)) + scale_color_manual(values = getPalette3(4)) + theme_classic(base_size = 12) + ggpubr::rotate_x_text(45) + labs(x = "", y = "prop of TCR repertoire")
ggsave("results/tcr/path_pr1_amount.pdf", width = 5, height = 4)

cytotoxic_markers     <- c("GZMK", "IFNG", "GZMA", "GZMB", "GZMM", "GZMH", "PRF1", "GNLY", "NCAM1")
pr1_seurat <- AddModuleScore(pr1_seurat, features = c("LAG3", "CTLA4", "PDCD1", "HAVCR2"), name = "exhaustion", nbin = 10)
pr1_seurat <- AddModuleScore(pr1_seurat, features = cytotoxic_markers, name = "cytotoxicity", nbin = 10)

pr1_seurat@meta.data %>% mutate(cluster = Idents(pr1_seurat)) %>%
  ggplot(aes(timepoint, cytotoxicity1, fill = timepoint)) + geom_violin(draw_quantiles = c(0.5), adjust = 1.5) + facet_wrap(~patient, scales = "free_x") +
#  ggpubr::stat_compare_means() +
  ggsignif::geom_signif(comparisons = list(c("baseline", "relapse"))) +
  ggsignif::geom_signif(comparisons = list(c("baseline", "6m"), c("6m", "12m"), c("baseline", "12m")), step_increase = 0.05) +

  scale_fill_manual(values = getPalette4(4)) + theme(legend.position = "none") +
  theme_classic(base_size = 12) + labs(x = "", y = "cytotoxicity score") + ggpubr::rotate_x_text(45) + facets_nice + theme(legend.position = "none")
ggsave("results/tcr/vln_cytotoxicity.pdf", width = 5, height = 4)


pr1_seurat@meta.data %>% mutate(cluster = Idents(pr1_seurat)) %>%
  ggplot(aes(timepoint, exhaustion1, fill = timepoint)) + geom_violin(draw_quantiles = c(0.5), adjust = 5) + facet_wrap(~patient, scales = "free_x") +
  #  ggpubr::stat_compare_means() +
  ggsignif::geom_signif(comparisons = list(c("baseline", "relapse"))) +
  ggsignif::geom_signif(comparisons = list(c("baseline", "6m"), c("6m", "12m"), c("baseline", "12m")), step_increase = 0.05) +

  scale_fill_manual(values = getPalette4(4)) + theme(legend.position = "none") +
  theme_classic(base_size = 12) + labs(x = "", y = "exhaustion score") + ggpubr::rotate_x_text(45) + facets_nice + theme(legend.position = "none")
ggsave("results/tcr/vln_exhaustion.pdf", width = 5, height = 4)


pr1_seurat@meta.data %>% mutate(cluster = Idents(pr1_seurat)) %>%
group_by(orig.ident) %>% summarise(n = n()) %>% left_join(total_df, by = "orig.ident") %>% left_join(patient_df, by = "orig.ident") %>%
  mutate(prop = n.x / n.y) %>%
  ggplot(aes(timepoint, prop, group = patient, color = patient, fill = patient)) + geom_path() +  geom_point(shape = 21, size = 3, color = "black") + theme_classic(base_size = 12) + labs(x = "") +
  scale_fill_manual(values = getPalette3(4)) + scale_color_manual(values = getPalette3(4)) + theme_classic(base_size = 12) + ggpubr::rotate_x_text(45) + labs(x = "", y = "prop of TCR repertoire")
ggsave("results/tcr/path_pr1_amount.pdf", width = 5, height = 4)



## Get DEGs
pt5_cells <- pr1_seurat@meta.data %>% filter(patient == "5") %>% pull(barcode)
pt7_cells <- pr1_seurat@meta.data %>% filter(patient == "7") %>% pull(barcode)

batch5_pr1_seurat <- subset(pr1_seurat, cells = pt5_cells)
batch7_pr1_seurat <- subset(pr1_seurat, cells = pt7_cells)

DEG_batch5 <- lapply(levels(Idents(batch5_pr1_seurat)), getDEGbyCluster, seurat_object = batch5_pr1_seurat)
DEG_batch5_df <- DEG_batch5 %>% rbindlist()
write.table(DEG_batch5_df, "results/tcr/deg_batch5_pr1.txt", sep = "\t", quote = F, row.names = F)

DEG_batch7 <- lapply(levels(Idents(batch7_pr1_seurat)), getDEGbyCluster, seurat_object = batch7_pr1_seurat)
DEG_batch7_df <- DEG_batch7 %>% rbindlist()
write.table(DEG_batch7_df, "results/tcr/deg_batch7_pr1.txt", sep = "\t", quote = F, row.names = F)

DEG_batch5_df %>% filter(gene %in% inhibitory_long)
DEG_batch7_df %>% filter(gene %in% inhibitory_long)

saveRDS(pr1_seurat, "results/pr1_seurat.rds")
