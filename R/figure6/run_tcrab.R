
Idents(cml_seurat) %>% levels()
dir.create("results/tcr/", showWarnings = F)

tcr_clusters <- c("1 CD4 naive x", "2 CD4 CM x", "3 CD8 effectory x",
                  "4 CD4 Cyto/Th1", "5 CD8 effectory/exhausted x", "7 T cells activated x",
                  "9 CD4 naive x", "10 CD4 treg x", "11 CD8 CM/naive x", "14 CD8 EM x", "16 CD4", "18 CD4 Acute Act")
tcr_clusters <- cml_seurat@meta.data %>% group_by(cluster, tcr = !is.na(new_clonotypes_id)) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% filter(tcr == T) %>% filter(prop > 0.4) %>% pull(cluster)

tcr_cells <- cml_seurat@meta.data %>% mutate(cluster = Idents(cml_seurat)) %>% filter(!is.na(new_clonotypes_id)) %>%
  filter(cluster %in% tcr_clusters) %>% pull(barcode)

cml_tcr_seurat <- subset(cml_seurat, cells = tcr_cells)
cml_tcr_seurat <- cml_tcr_seurat %>% getLatentUMAP() %>% fixSeurat()

p <- DimPlot(cml_tcr_seurat, label = T, repel = T, cols = getPalette4(12)) + theme_bw(base_size = 12) + theme(legend.position = "none")
ggsave(plot = p, "results/tcr/latent_umap.png", width = 5, height = 4)

cml_tcr_seurat@meta.data %>% mutate(cluster = Idents(cml_tcr_seurat)) %>%
  mutate(relapse = factor(relapse, levels = c("None", "Slow", "Fast"))) %>%
  group_by(timepoint, relapse, cluster, new_clonotypes_id) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  summarise(diversity = vegan::diversity(prop), gini = ineq::Gini(prop)) %>%

  ggplot(aes(timepoint, gini, label = cluster, fill = relapse, group = relapse)) +
    geom_point(shape = 21, size = 3) + labs(y = "Gini index", x = "time point") + geom_path() + facet_wrap(~cluster, scales = "free_y", ncol = 5) +
    scale_fill_manual(values = getPalette3(4)[c(1,3,2)]) +
    facets_nice + ggpubr::rotate_x_text(45)
ggsave("results/tcr/path_cluster.png", width = 9, height = 4.5)


cml_tcr_seurat@meta.data %>% mutate(cluster = Idents(cml_tcr_seurat)) %>%
  filter(!is.na(new_clonotypes_id)) %>%
  group_by(cluster, new_clonotypes_id) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  summarise(diversity = vegan::diversity(prop), gini = ineq::Gini(prop)) %>%

  ggplot(aes(gini,diversity,label=cluster)) + geom_point(shape = 21, fill = "lightgrey", size = 3) + ggrepel::geom_text_repel() + labs(x = "Gini index", y = "Shannon diveristy")
ggsave("results/tcr/scatter_gini_clonality.pdf", width = 5, height = 4)

