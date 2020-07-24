

## No need to modify
#### =====================

## Init output options
output_dir       <- paste0("results/tcrgp/plots/")
dir.create(output_dir, showWarnings = F)


cml_seurat$patient   <- substr(cml_seurat$orig.ident, 1, 7)
cml_seurat$timepoint <- substr(cml_seurat$orig.ident, 9, nchar(cml_seurat$orig.ident))

## TCRGP results
cml_seurat@meta.data %>% filter(target != "multi") %>% group_by(pred_epitope) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>%
  mutate(pred_epitope = factor(pred_epitope, levels = pred_epitope[order((freq))])) %>%

  ggplot(aes(pred_epitope, freq, label = round(freq, 3), fill = pred_epitope)) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette3(10)) + coord_flip() +
  theme_bw() + theme(legend.position = "none") + geom_text() + labs(x = "TCRGP predicted epitopes") #ggrepel::geom_text_repel()
ggsave(paste0(output_dir, "bar_clonotype_tcrgp.pdf"), width = 6, height = 4)


cml_seurat@meta.data %>% filter(target != "multi") %>% group_by(patient, pred_epitope) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>%
  mutate(pred_epitope = factor(pred_epitope, levels = pred_epitope[order((freq))])) %>%

  ggplot(aes(pred_epitope, freq, label = round(freq, 3), fill = pred_epitope)) + geom_bar(stat = "identity") + scale_fill_manual(values = getPalette3(10)) + coord_flip() +
  theme_bw() + theme(legend.position = "none") + geom_text() + labs(x = "TCRGP predicted epitopes") + facet_wrap(~patient)
ggsave(paste0(output_dir, "bar_clonotype_tcrgp_per_pt.pdf"), width = 12, height = 8)




## Epitope specifities during time points?
cml_seurat@meta.data %>% filter(target != "multi") %>% group_by(timepoint, patient, pred_epitope) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>%
  ggplot(aes(timepoint, freq, fill = pred_epitope, group = pred_epitope)) + geom_path() + geom_point(pch = 21, size = 5) + scale_fill_manual(values = getPalette3(10)) + theme_bw() + facet_wrap(~patient)
ggsave(paste0(output_dir, "path_epitope1.pdf"), width = 12, height = 8)

cml_seurat@meta.data %>% filter(target != "multi") %>% group_by(timepoint, patient, pred_epitope) %>% summarise(n = n()) %>% mutate(freq = n/sum(n)) %>%
  filter(pred_epitope != "") %>%
  ggplot(aes(timepoint, freq, fill = pred_epitope, group = pred_epitope)) + geom_path() + geom_point(pch = 21, size = 5) + scale_fill_manual(values = getPalette3(10)) + theme_bw() + facet_wrap(~patient)
ggsave(paste0(output_dir, "path_epitope2.pdf"), width = 12, height = 8)




# Analyse the behaviour in UMAP from the clusters
viz_df    <- data.frame(cml_seurat[["umap"]]@cell.embeddings, "cluster" = Idents(cml_seurat), cml_seurat@meta.data) #%>%  dplyr::rename(umap_1 = 1, umap_2 = 2)
umap_mean <- data.frame(aggregate(umap_1 ~ cluster, viz_df, median), umap_2 = aggregate(umap_2 ~ cluster, viz_df, median)[,2])
nClusters <- viz_df$cluster %>% unique %>% length()

## Plot UMAPs with TCRGP predictions highlighted
p <- ggplot() +
  geom_point(data = viz_df, aes(x = umap_1, y = umap_2, color = cluster), size = 0.8, color = "lightgrey") +
  geom_point(data = subset(viz_df, target == "anti-viral"), aes(x = umap_1, y = umap_2, fill = pred_epitope), size = 2, shape = 21) +

  stat_ellipse(data = viz_df, geom = "polygon", aes(x = umap_1, y = umap_2, group = cluster), color = "lightgrey", fill = "lightgrey", alpha = 0.1, lty = "dotted") +
  ggrepel::geom_label_repel(data = umap_mean, aes(x = umap_1, y = umap_2, label = cluster), size = 3, color = "black") + theme_void() + add_guide + scale_fill_manual(values =getPalette3(10))
ggsave(plot = p, paste0(output_dir, "umap_tcrgp.png"), width = 10, height = 8)



## Every epitope
dir.create(paste0(output_dir, "perEpitope/"), showWarnings = F)

colors   <- getPalette(28)
epitopes <- viz_df %>% filter(target != "multi") %>% group_by(pred_epitope) %>% summarise(n = n()) %>% dplyr::select(pred_epitope) %>% as.vector()
i <- 1

for(epitope in as.character(epitopes$pred_epitope)[-1]){

  epitope_df <- fisherEpitope(seurat_object = cml_seurat, epitope = epitope) %>% filter(sigf == "Sigf")

  message(epitope)
  viz_df_temp <- viz_df %>% filter(target == "anti-viral") %>% filter(pred_epitope == epitope)

  p <- ggplot() +
    geom_point(data = viz_df,  aes(x = umap_1, y = umap_2), color = "lightgrey", size = 0.8) +
    geom_point(data = viz_df_temp, aes(x = umap_1, y = umap_2, color = cluster), size = 1.5) +

    ggrepel::geom_label_repel(data = umap_mean, aes(x = umap_1, y = umap_2, color = cluster, label = cluster), size = 10) +

    stat_ellipse(data = viz_df, aes(x = umap_1, y = umap_2, group = cluster, color = cluster), linetype = "dotted", size = 1) +
    stat_ellipse(data = subset(viz_df, cluster %in% epitope_df$cluster), aes(x = umap_1, y = umap_2, group = cluster), color = "darkred", size = 1.5, linetype = "dotted") +

    scale_color_manual(values = getPalette2(nClusters)) +
    theme_void() + theme(legend.position = "none") + labs(title = epitope)

  ggsave(plot = p, paste0(output_dir, "perEpitope/", epitope, ".png"), width = 14, height = 12)

  i <- i + 1

}


