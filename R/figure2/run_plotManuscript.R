
dir.create("results/figure2/", showWarnings = F)
cml_seurat <- readRDS("latest/data/seurat_objects/cml_fin_19feb.rds")

Idents(cml_seurat) <- Idents(cml_seurat) %>% extractClusterNumber() %>% as.numeric %>% as.factor() %>% getClusterPhenotypesCml()
cml_seurat$cluster <- Idents(cml_seurat)

##  QC
cml_seurat@meta.data %>% ggplot(aes(cluster, nCount_RNA, fill = cluster)) + geom_violin() + scale_fill_manual(values = getPalette2(25)) + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) + labs(x = "")
ggsave("results/figure2/vln_nCount_rna.pdf", width = 7, height = 4)

cml_seurat@meta.data %>% ggplot(aes(cluster, nFeature_RNA, fill = cluster)) + geom_violin() + scale_fill_manual(values = getPalette2(25)) + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) + labs(x = "")
ggsave("results/figure2/vln_nFeature_RNA.pdf", width = 7, height = 4)

cml_seurat@meta.data %>% ggplot(aes(cluster, percent.mt, fill = cluster)) + geom_violin() + scale_fill_manual(values = getPalette2(25)) + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) + labs(x = "")
ggsave("results/figure2/vln_percent.mt.pdf", width = 7, height = 4)

cml_seurat@meta.data %>% ggplot(aes(cluster, percent.rb, fill = cluster)) + geom_violin() + scale_fill_manual(values = getPalette2(25)) + theme_classic(base_size = 12) + theme(legend.position = "none") + ggpubr::rotate_x_text(45) + labs(x = "")
ggsave("results/figure2/vln_percent.ribo.pdf", width = 7, height = 4)



## Plot most notable markers
p <- DimPlot(cml_seurat, cols = getPalette2(24), label = T, repel = T) + theme_bw(base_size = 12) + theme(legend.position = "none") + labs(x = "UMAP1", y = "UMAP2")
ggsave(plot = p, "results/figure2/latent_umap.png", width = 6, height = 5)

p <- DotPlot(cml_seurat, features = rev(unique(big_markers)), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y = "")+ theme(legend.position = "top")
ggsave(plot = p, "results/figure2/dotplot_big_markers.png", width = 14, height = 6)

p <- DotPlot(cml_seurat, features = rev(unique(do.call(zhang_cd4_markers, what = "c"))), cols = "RdYlBu") + labs(x = "", y = "cluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(y = "") + theme(legend.position = "top")
ggsave(plot = p, "results/figure2/dotplot_zhang_cd4_markers.png", width = 14, height = 6)

p <- DotPlot(cml_seurat, features = rev(unique(c("CD3E", "NCAM1", "GZMK", "SELL", "XCL1", "XCL2", "KLRC1", "IL7R", "LTB", "FCGR3A", "GZMA", "GZMB", "GZMH", "GZMM", "KLRC2", "ZEB2", "KLF2", "PRDM1", "GZMH", "LAG3"))), cols = "RdYlBu") +
  ggpubr::rotate_x_text(angle = 90) + labs(x = "", y  = "") + theme(legend.position = "top")
ggsave(plot = p, "results/figure2/dotplot_dufva_markers.png", width = 9, height = 6)

guo_genes <- guo_markers %>% do.call(what = "c") %>% as.character() %>% rev %>% unique()
p <- DotPlot(cml_seurat, features = guo_genes, cols = "RdYlBu") +
  ggpubr::rotate_x_text(angle = 90) + labs(x = "", y  = "")
ggsave(plot = p, "results/figure2/dotplot_guo_markers.png", width = 12, height = 6)

p <- cml_seurat@meta.data %>% group_by(cluster, singler_hpca_pred) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% filter(prop > 0.1) %>%
  ggplot(aes(singler_hpca_pred, prop, fill = singler_hpca_pred, label = singler_hpca_pred)) + geom_bar(stat = "identity") + add_guide + scale_fill_manual(values = getPalette(12)) + facet_wrap(~cluster, ncol = 6) + coord_flip() +
  theme_bw(base_size = 12) + labs(x = "") + theme(legend.position = "none") + geom_hline(yintercept = 0.5, linetype = "dotted") + ggrepel::geom_text_repel()
ggsave(plot = p, "results/figure2/bar_predictions_hpca.png", width = 12, height = 8)

p <- cml_seurat@meta.data %>% group_by(cluster, singler_blueprint_pred) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% filter(prop > 0.1) %>%
  ggplot(aes(singler_blueprint_pred, prop, fill = singler_blueprint_pred, label = singler_blueprint_pred)) + geom_bar(stat = "identity") + add_guide + scale_fill_manual(values = getPalette(24)) + facet_wrap(~cluster, ncol = 6) + coord_flip() +
  theme_bw(base_size = 12) + labs(x = "") + theme(legend.position = "none") + geom_hline(yintercept = 0.5, linetype = "dotted") + ggrepel::geom_text_repel() + ylim(values = c(0,1))
ggsave(plot = p, "results/figure2/bar_predictions_blueprint.png", width = 12, height = 8)





clusters <- cml_seurat@meta.data %>% mutate(clusters = Idents(cml_seurat)) %>% group_by(clusters, !is.na(new_clonotypes_id)) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% filter(`!is.na(new_clonotypes_id)` == T) %>% filter(prop > 0.4) %>% pull(clusters) %>% as.character() %>% unique()

cml_seurat@meta.data %>% mutate(cluster = Idents(cml_seurat)) %>%
  filter(cluster %in% clusters) %>%
  filter(!is.na(new_clonotypes_id)) %>%
  group_by(cluster, new_clonotypes_id) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  summarise(diversity = vegan::diversity(prop), gini = ineq::Gini(prop)) %>%

  ggplot(aes(gini,diversity,label=cluster)) + geom_point(shape = 21, fill = "lightgrey", size = 3) + ggrepel::geom_text_repel() + labs(x = "Gini index", y = "Shannon diveristy")
ggsave("results/figure2/scatter_gini_clonality.pdf", width = 5, height = 4)




cml_seurat@meta.data %>% mutate(clusters = Idents(cml_seurat)) %>%
  group_by(clusters, !is.na(new_clonotypes_id)) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  filter(`!is.na(new_clonotypes_id)` == T) %>%

  ggplot(aes(reorder(clusters, prop), prop)) + geom_bar(stat = "identity") + ggpubr::rotate_x_text(45) +
  geom_hline(yintercept = 0.5, linetype = "dotted") + labs(x = "", y = "proportion of cells with TCR") + ylim(c(0,1))
ggsave("results/figure2/bar_tcrab_cluster.pdf", width = 5, height = 4)

cml_seurat@meta.data %>% mutate(clusters = Idents(cml_seurat)) %>%
  group_by(clusters) %>% summarise(n = n()) %>%

  ggplot(aes(reorder(clusters, n), n, fill = clusters, label = n)) + geom_bar(stat = "identity") + coord_flip() + scale_fill_manual(values = getPalette(24)) + theme(legend.position = "none") + labs(x = "", y = "nCells") +
  geom_text()
ggsave("results/figure2/bar_n_cluster.pdf", width = 5, height = 4)


## From public
df <- cml_seurat@meta.data %>% left_join(dplyr::select(public_seurat@meta.data, barcode, cluster) %>% dplyr::rename(public_cluster = cluster), by = "barcode")

p <- DimPlot(cml_seurat, cols = getPalette4(25), group.by = "public_cluster", label = T, repel = T) + theme_bw(base_size = 12) + theme(legend.position = "none") + labs(x = "UMAP1", y = "UMAP2")
ggsave(plot = p, "results/figure2/latent_umap_public.png", width = 7, height = 6)



## Find DE-genes between all the clusters
all_markers <- FindAllMarkers(cml_seurat, test.use = "t", verbose = T, max.cells.per.ident = 1e3)
fwrite(all_markers, "results/figure2/all_markers_1e3.txt", sep = "\t", quote = F, row.names = F)
all_markers <- fread("results/figure2/all_markers_1e3.txt") %>% mutate(cluster = cluster %>% extractClusterNumber() %>% as.numeric %>% as.factor() %>% getClusterPhenotypesCml())
fwrite(all_markers, "results/figure2/all_markers_1e3.txt", sep = "\t", quote = F, row.names = F)


set.seed(123)
top10  <- all_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
cells <- cml_seurat@meta.data %>% mutate(cluster = Idents(cml_seurat)) %>% group_by(cluster) %>% sample_n(50) %>% pull(barcode)
p <- DoHeatmap(cml_seurat, cells = cells, features = top10$gene, angle = 90, group.colors = getPalette(24)) + theme(legend.position = "none")
ggsave(plot = p, "results/figure2/heatmap_top10_markers.png", width = 16, height = 20)







cml_seurat$timepoint <- substr(cml_seurat$orig.ident, 9, nchar(cml_seurat$orig.ident))
cml_seurat$timepoint <- factor(as.character(cml_seurat$timepoint), levels = c("baseline", "6m", "12m", "relapse"))

## Lolliplot at base line
df <- cml_seurat@meta.data %>% mutate(cluster = as.character(Idents(cml_seurat))) %>%
  filter(timepoint %in% c("baseline")) %>%
  mutate(controller = ifelse(relapse == "Fast", "relapse", "control")) %>%
  group_by(controller,cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% ungroup()

df_controller <- df %>% filter(controller != "relapse")
df_relapse    <- df %>% filter(controller == "relapse")

df_relapse %>% left_join(df_controller, by = "cluster") %>% mutate(log2fc = log2(prop.y / prop.x)) %>%
  mutate(fill = ifelse(log2fc > 0, "up", "down")) %>%
  ggplot(aes(reorder(cluster, log2fc), log2fc,fill=fill)) + geom_bar(stat = "identity") + coord_flip() +
  geom_hline(yintercept = -1, linetype = "dotted") + geom_hline(yintercept = 1, linetype = "dotted")  + theme_classic(base_size = 12) + theme(legend.position = "none") + labs(x = "") + ylim(c(-2.5,2.5)) + scale_fill_manual(values = getPalette(6))
ggsave("results/figure2/bar_fc_baseline_contoller.pdf", width = 5, height = 4)

df <- cml_seurat@meta.data %>% mutate(cluster = as.character(Idents(cml_seurat))) %>%
  filter(timepoint %in% c("baseline", "6m") & relapse != "Fast") %>%

  mutate(timepoint = ifelse(timepoint == "baseline", "pre", "post")) %>%
  mutate(timepoint = factor(timepoint, levels = c("pre", "post"))) %>%
  group_by(timepoint,cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% ungroup()

df_pre  <- df %>% filter(timepoint == "pre")
df_post <- df %>% filter(timepoint == "post")

df_pre %>% left_join(df_post, by = "cluster") %>% mutate(log2fc = log2(prop.y / prop.x)) %>%
  mutate(fill = ifelse(log2fc > 0, "up", "down")) %>%
  ggplot(aes(reorder(cluster, log2fc), log2fc,fill=fill)) + geom_bar(stat = "identity") + coord_flip() +
  geom_hline(yintercept = -1, linetype = "dotted") + geom_hline(yintercept = 1, linetype = "dotted")  + theme_classic(base_size = 12) + theme(legend.position = "none") + labs(x = "") + ylim(c(-2.5,2.5)) + scale_fill_manual(values = getPalette(8))
ggsave("results/figure2/bar_fc_cluster.pdf", width = 5, height = 4)


## Line plot
p <- cml_seurat@meta.data %>% mutate(cluster = (Idents(cml_seurat))) %>%
  group_by(timepoint,relapse,cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% ungroup() %>%
  mutate(relapse = factor(relapse, levels = c("None", "Slow", "Fast"))) %>%

  ggplot(aes(as.factor(timepoint),prop,group=relapse,color=relapse)) + geom_path(lwd = 1.5) + facet_wrap(~cluster, scales = "free_y", ncol = 4) + theme_bw(base_size = 12) + facets_nice + scale_color_manual(values = getPalette3(4)[c(1,3,2)]) +
  labs(x = "", y = "prop") + ggpubr::rotate_x_text(45)
ggsave(plot = p, "results/figure2/path_cluster_relapse.pdf", width = 8, height = 12)


p <- cml_seurat@meta.data %>% mutate(cluster = (Idents(cml_seurat))) %>%
  mutate(controller = ifelse(relapse == "Fast", "relapse", "control")) %>%
  filter(controller == "control" & timepoint != "relapse") %>%
  group_by(timepoint,controller,cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% ungroup() %>%
  mutate(celltype = extractCoarsePhenotype(cluster)) %>%
  mutate(celltype = ifelse(celltype %in% c("B-cells", "innate", "low", "pDC", "T"), "other", celltype)) %>%

  ggplot(aes(as.factor(timepoint),prop,group=cluster,color=cluster)) +
  geom_path(lwd = 1.5) + facet_wrap(~celltype, scales = "free_y", ncol = 5) + theme_bw(base_size = 12) +
  facets_nice + scale_color_manual(values = getPalette3(24)) +
  labs(x = "", y = "prop") + ggpubr::rotate_x_text(45)
ggsave(plot = p, "results/figure2/path_cluster_celltype_controller.pdf", width = 14, height = 4)


p <- cml_seurat@meta.data %>% mutate(cluster = (Idents(cml_seurat))) %>%
  mutate(controller = ifelse(relapse == "Fast", "relapse", "control")) %>%
  filter(controller != "control") %>%
  group_by(timepoint,controller,cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% ungroup() %>%
  mutate(celltype = extractCoarsePhenotype(cluster)) %>%
  mutate(celltype = ifelse(celltype %in% c("B-cells", "innate", "low", "pDC", "T"), "other", celltype)) %>%

  ggplot(aes(as.factor(timepoint),prop,group=cluster,color=cluster)) +
  geom_path(lwd = 1.5) + facet_wrap(~celltype, scales = "free_y", ncol = 5) + theme_bw(base_size = 12) +
  facets_nice + scale_color_manual(values = getPalette3(24)) +
  labs(x = "", y = "prop") + ggpubr::rotate_x_text(45)
ggsave(plot = p, "results/figure2/path_cluster_celltype_relapse.pdf", width = 14, height = 4)


patient_df <- cml_seurat@meta.data %>% group_by(orig.ident, relapse) %>% summarise(n = n()) %>% ungroup %>% dplyr::select(-n)
diet_cml_nk_seurat$timepoint <- factor(as.character(diet_cml_nk_seurat$timepoint), levels = c("baseline", "6m", "12m", "relapse"))

diet_cml_nk_seurat@meta.data %>% mutate(cluster = (Idents(diet_cml_nk_seurat))) %>%
  filter(project %in% c("CML on TKI", "CML off TKI")) %>%
  left_join(patient_df) %>%
  group_by(timepoint,relapse,cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% ungroup() %>%
  mutate(relapse = factor(relapse, levels = c("None", "Slow", "Fast"))) %>%
  ggplot(aes(as.factor(timepoint),prop,group=relapse,color=relapse)) + geom_path(lwd = 1.5) + facet_wrap(~cluster, scales = "free_y", ncol = 4) + theme_bw(base_size = 12) + facets_nice + scale_color_manual(values = getPalette3(4)[c(1,3,2)]) +
  labs(x = "", y = "prop") + ggpubr::rotate_x_text(45)
ggsave("results/figure2/path_nk_cluster_celltype_relapse.pdf", width = 8, height = 4)


##### DEGs
DEG_cluster_no_df   <- read.delim("results/figure2/deg_no_rp_cluster.txt")    %>% mutate(cluster = cluster %>% extractClusterNumber() %>% as.numeric() %>% as.factor() %>% getClusterPhenotypesCml())
DEG_cluster_slow_df <- read.delim("results/figure2/deg_no_slow_cluster.txt")  %>% mutate(cluster = cluster %>% extractClusterNumber() %>% as.numeric() %>% as.factor() %>% getClusterPhenotypesCml())
DEG_cluster_fast_df <- read.delim("results/figure2/deg_no_fast_cluster.txt")  %>% mutate(cluster = cluster %>% extractClusterNumber() %>% as.numeric() %>% as.factor() %>% getClusterPhenotypesCml())



## Bar plot
no_relapse_2v1   <- DEG_cluster_no_df   %>% filter(timepoint == "2v1") %>% mutate(cluster = reorderClusters(cluster), type = "None")
slow_relapse_2v1 <- DEG_cluster_slow_df %>% filter(timepoint == "2v1") %>% mutate(cluster = reorderClusters(cluster), type = "Slow")
fast_relapse_2v1 <- DEG_cluster_fast_df %>% filter(timepoint == "4v1") %>% mutate(cluster = reorderClusters(cluster), type = "Fast")

df <- rbind(no_relapse_2v1, slow_relapse_2v1, fast_relapse_2v1) %>% filter(cluster == "3 CD8 effectory") %>% filter(gene %in% c(cytotoxic_markers, inhibitory_long, guo_markers, marker_genes))
heatmap <- dcast(gene~type, data = df, value.var = "avg_logFC")
rownames(heatmap) <- heatmap$gene
heatmap <- heatmap[,-1]
heatmap[is.na(heatmap)] <- 0
heatmap <- heatmap[,c(2,3,1)]

p <- pheatmap::pheatmap(heatmap, cluster_cols = F, color = rev(colorRampPalette(brewer.pal(11, "RdBu"))(100)), breaks = seq(-1, 1, 0.02)) #, annotation_col = "315")
ggsave(plot = print(p), "results/figure2/heatmap_cd8_deg.pdf", width = 3, height = 5)




rbind(no_relapse_2v1, slow_relapse_2v1, fast_relapse_2v1) %>%
  mutate(type = factor(type, levels = c("None", "Slow", "Fast"))) %>%
  group_by(cluster, type) %>% summarise(n = n()) %>%
  ggplot(aes(type,n,fill=type)) + geom_bar(stat = "identity") + facet_wrap(~cluster) + theme_classic(base_size = 12) + scale_fill_manual(values = getPalette3(4)[c(1,3,2)]) + labs(x = "", y = "nDEGs") + ggpubr::rotate_x_text(45) +
  facets_nice
ggsave("results/figure2/bar_deg_cluster_2v1.pdf", width = 8, height = 6)


rbind(no_relapse_2v1, slow_relapse_2v1, fast_relapse_2v1) %>%
  mutate(type = factor(type, levels = c("None", "Slow", "Fast"))) %>%
  group_by(cluster, type) %>% summarise(n = n()) %>%
  ggplot(aes(type,n,color=type,label=n)) +
  geom_segment(aes(x = type, xend = type, y = 0, yend = n), color = "black") + geom_point(size = 10) +
  geom_text(color = "white") + facet_wrap(~cluster) + theme_classic(base_size = 12) + scale_color_manual(values = getPalette3(4)[c(1,3,2)]) + labs(x = "", y = "nDEGs") + ggpubr::rotate_x_text(45) +
  facets_nice + theme(legend.position = "none")
ggsave("results/figure2/lolliplot_deg_cluster_2v1.pdf", width = 9, height = 7)




rbind(no_relapse_2v1, slow_relapse_2v1, fast_relapse_2v1) %>%
  mutate(type = factor(type, levels = c("None", "Slow", "Fast"))) %>%
  group_by(cluster, type) %>% summarise(n = n()) %>%
  ggplot(aes(reorder(cluster, n),n,fill=type)) + geom_bar(stat = "identity") + facet_wrap(~type) + theme_classic(base_size = 12) + scale_fill_manual(values = getPalette3(4)[c(1,3,2)]) + labs(x = "", y = "nDEGs") + ggpubr::rotate_x_text(45) + coord_flip() + facets_nice + theme(legend.position = "none")
ggsave("results/figure2/bar_deg_cluster_2v1_2.pdf", width = 8, height = 6)


rbind(no_relapse_2v1, slow_relapse_2v1, fast_relapse_2v1) %>%
  mutate(type = factor(type, levels = c("None", "Slow", "Fast"))) %>%
  group_by(cluster, type) %>% summarise(n = n()) %>%
  ggplot(aes(reorder(cluster, n),n,color=type,label=n)) +
  geom_segment(aes(x = reorder(cluster, n), xend = reorder(cluster, n), y = 0, yend = n), color = "black") + geom_point(size = 10) +
  geom_text(color = "white") +
  facet_wrap(~type) + theme_classic(base_size = 12) + scale_color_manual(values = getPalette3(4)[c(1,3,2)]) + labs(x = "", y = "nDEGs") + ggpubr::rotate_x_text(45) + coord_flip() + facets_nice + theme(legend.position = "none")
ggsave("results/figure2/lolliplot_deg_cluster_2v1_2.pdf", width = 8, height = 6)



df <- rbind(no_relapse_2v1, slow_relapse_2v1, fast_relapse_2v1) %>%
  group_by(cluster) %>% summarise(total = n())

rbind(no_relapse_2v1, slow_relapse_2v1, fast_relapse_2v1) %>%
  mutate(type = factor(type, levels = c("None", "Slow", "Fast"))) %>%
  group_by(cluster, type) %>% summarise(n = n()) %>%
  left_join(df) %>%

  ggplot(aes(reorder(cluster, total), n,color=type,label=n)) +
  geom_segment(aes(x = reorder(cluster, total), xend = reorder(cluster, total), y = 0, yend = n), color = "black") + geom_point(size = 10) +
  geom_text(color = "white") +
  facet_wrap(~type) + theme_classic(base_size = 12) + scale_color_manual(values = getPalette3(4)[c(1,3,2)]) + labs(x = "", y = "nDEGs") + ggpubr::rotate_x_text(45) + coord_flip() + facets_nice + theme(legend.position = "none")
ggsave("results/figure2/lolliplot_deg_cluster_2v1_2.pdf", width = 8, height = 6)




df <- rbind(no_relapse_2v1, slow_relapse_2v1, fast_relapse_2v1) %>%
  group_by(cluster) %>% summarise(total = n())

rbind(no_relapse_2v1, slow_relapse_2v1, fast_relapse_2v1) %>%
  mutate(type = factor(type, levels = c("None", "Slow", "Fast"))) %>%
  group_by(cluster, type) %>% summarise(n = n()) %>%
  left_join(df) %>%

  ggplot(aes(reorder(cluster, total), n,color=type,label=n,group=cluster)) +
  geom_point(size = 3) +
  geom_path(color="grey") +
  # geom_text(color = "white") +
  theme_classic(base_size = 12) + scale_color_manual(values = getPalette3(4)[c(1,3,2)]) + labs(x = "", y = "nDEGs") + ggpubr::rotate_x_text(45) + facets_nice #+ theme(legend.position = "none")
ggsave("results/figure2/line_deg_cluster_2v1.pdf", width = 6, height = 4)






## Enrichment; hallmark
output_dir = "results/figure2/enrichment/"
dir.create(output_dir, showWarnings = F)
universe_df <- rownames(cml_seurat)

hallmark_no_df_up <- lapply(levels(Idents(cml_seurat)), function(x){
  df <- no_relapse_2v1 %>% filter(cluster == x & direction == "up") %>% getHypergeometric(universe_df = universe_df, term_df = hallmark)
  if(!is.null(df)){df %>% mutate(cluster = x)}
}) %>% rbindlist() %>% filter(p.adjust < 0.05) %>% mutate(type = "None")

hallmark_slow_df_up <- lapply(levels(Idents(cml_seurat)), function(x){
  df <- slow_relapse_2v1 %>% filter(cluster == x & direction == "up") %>% getHypergeometric(universe_df = universe_df, term_df = hallmark)
  if(!is.null(df)){df %>% mutate(cluster = x)}
}) %>% rbindlist() %>% filter(p.adjust < 0.05) %>% mutate(type = "Slow")

hallmark_fast_df_up <- lapply(levels(Idents(cml_seurat)), function(x){
  df <- fast_relapse_2v1 %>% filter(cluster == x & direction == "up") %>% getHypergeometric(universe_df = universe_df, term_df = hallmark)
  if(!is.null(df)){df %>% mutate(cluster = x)}
}) %>% rbindlist() %>% filter(p.adjust < 0.05) %>% mutate(type = "Fast")

rbind(hallmark_no_df_up, hallmark_slow_df_up, hallmark_fast_df_up) %>%
  fwrite("results/figure2/pathways_up.txt", sep = "\t", quote = F, row.names = F)

hallmark_no_df_down <- lapply(levels(Idents(cml_seurat)), function(x){
  df <- no_relapse_2v1 %>% filter(cluster == x & direction == "down") %>% getHypergeometric(universe_df = universe_df, term_df = hallmark)
  if(!is.null(df)){df %>% mutate(cluster = x)}
}) %>% rbindlist() %>% filter(p.adjust < 0.05) %>% mutate(type = "None")

hallmark_slow_df_down <- lapply(levels(Idents(cml_seurat)), function(x){
  df <- slow_relapse_2v1 %>% filter(cluster == x & direction == "down") %>% getHypergeometric(universe_df = universe_df, term_df = hallmark)
  if(!is.null(df)){df %>% mutate(cluster = x)}
}) %>% rbindlist() %>% filter(p.adjust < 0.05) %>% mutate(type = "Slow")

hallmark_fast_df_down <- lapply(levels(Idents(cml_seurat)), function(x){
  df <- fast_relapse_2v1 %>% filter(cluster == x & direction == "down") %>% getHypergeometric(universe_df = universe_df, term_df = hallmark)
  if(!is.null(df)){df %>% mutate(cluster = x)}
}) %>% rbindlist() %>% filter(p.adjust < 0.05) %>% mutate(type = "Fast")

rbind(hallmark_no_df_down, hallmark_slow_df_down, hallmark_fast_df_down) %>%
  fwrite("results/figure2/pathways_down.txt", sep = "\t", quote = F, row.names = F)


df <- rbind(hallmark_no_df_up, hallmark_slow_df_up, hallmark_fast_df_up) %>%
  group_by(cluster, ID) %>% summarise(n = n()) # %>% filter(n > 1)
heatmap <- dcast(cluster~ID,data=df,values="n")
rownames(heatmap) <- heatmap$cluster
heatmap <- heatmap[,-1]
colnames(heatmap) <- gsub("HALLMARK\\_", "", colnames(heatmap))
heatmap[is.na(heatmap)] <- 0

p <- pheatmap::pheatmap(heatmap)
ggsave(plot = print(p), "results/figure2/heatmap_pathways_up.pdf", width = 8, height = 7)


df <- rbind(hallmark_no_df_up, hallmark_slow_df_up) %>%
  group_by(cluster, ID) %>% summarise(n = n()) # %>% filter(n > 1)
heatmap <- dcast(cluster~ID,data=df,values="n")
rownames(heatmap) <- heatmap$cluster
heatmap <- heatmap[,-1]
colnames(heatmap) <- gsub("HALLMARK\\_", "", colnames(heatmap))
heatmap[is.na(heatmap)] <- 0

p <- pheatmap::pheatmap(heatmap)
ggsave(plot = print(p), "results/figure2/heatmap_pathways_up_controller.pdf", width = 8, height = 7)



df <- rbind(hallmark_fast_df_up) %>%
  group_by(cluster, ID) %>% summarise(n = n()) # %>% filter(n > 1)
heatmap <- dcast(cluster~ID,data=df,values="n")
rownames(heatmap) <- heatmap$cluster
heatmap <- heatmap[,-1]
colnames(heatmap) <- gsub("HALLMARK\\_", "", colnames(heatmap))
heatmap[is.na(heatmap)] <- 0

p <- pheatmap::pheatmap(heatmap)
ggsave(plot = print(p), "results/figure2/heatmap_pathways_up_fast.pdf", width = 8, height = 7)







df <- rbind(hallmark_no_df_down, hallmark_slow_df_down, hallmark_fast_df_down) %>%
  group_by(cluster, ID) %>% summarise(n = n()) # %>% filter(n > 1)
heatmap <- dcast(cluster~ID,data=df,values="n")
rownames(heatmap) <- heatmap$cluster
heatmap <- heatmap[,-1]
colnames(heatmap) <- gsub("HALLMARK\\_", "", colnames(heatmap))
heatmap[is.na(heatmap)] <- 0
pheatmap::pheatmap(heatmap)

p <- pheatmap::pheatmap(heatmap)
ggsave(plot = print(p), "results/figure2/heatmap_pathways_down.pdf", width = 8, height = 7)

