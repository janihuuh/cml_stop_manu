
a <- DimPlot(cml2_seurat, group.by = "my.demux", cols = getPalette(10))
b <- DimPlot(cml2_seurat, group.by = "covariate", cols = getPalette(5))
c <- DimPlot(cml2_seurat, group.by = "expanded", cols = getPalette3(4))
d <- DimPlot(cml2_seurat, group.by = "lgals9_ko", cols = getPalette3(4))
a+b+c+d
FeaturePlot(cml2_seurat, features = "GZMB")


cml2_seurat <- cml_seurat_new_filt2



### Show that KO is working?
cml2_seurat_tumor <- subset(cml2_seurat, is_nk != "NK")
cml2_seurat_tumor <- subset(cml2_seurat_tumor, my.demux != "NE-NK-only")
cml2_seurat_tumor <- cml2_seurat_tumor %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml2_seurat_tumor), nPCs = 11)
cml2_seurat_tumor@reductions$pca@stdev %>% plot()

cml2_seurat_tumor <- AddModuleScore(cml2_seurat_tumor, features = list(c("LGALS9")), name = "LGALS9")

cml2_seurat_tumor@meta.data %>% 
  ggplot(aes(my.demux,LGALS91,fill=my.demux)) + geom_violin(draw_quantiles = 0.5, adjust = 3) + geom_jitter(size=0.001,alpha=0.3) +
  ggpubr::stat_compare_means(label="p", label.y.npc = 0.9) + ggpubr::rotate_x_text(angle = 45) + 
  #scale_fill_manual(values=c("blue3","green3","green4")) +
  NoLegend() + labs(x="",y="LGALS9") 
ggsave("results/functional/hto2/vln_lgals_ko_total.pdf", width = 4, height = 3)

cml2_seurat_tumor@meta.data %>% group_by(my.demux, lgals9_pos = LGALS91>0) %>% tally() %>% mutate(prop=n/sum(n)) %>% 
  filter(lgals9_pos) %>% 
  
  # ggplot(aes(reorder(my.demux,-prop),prop*100,fill=my.demux)) + geom_bar(stat="identity") +
  ggplot(aes(my.demux,prop*100,fill=my.demux)) + geom_bar(stat="identity") +
  
  # ggpubr::stat_compare_means(label="p") + 
  ggpubr::rotate_x_text(angle = 45) +
  #scale_fill_manual(values=c("blue3","green3","green4")) +
  labs(x="",y="% of LGALS9 pos cells") + NoLegend() 
ggsave("results/functional/hto2/bar_lgals_ko_total.pdf", width = 3, height = 4)





cml2_seurat_tumor_lgal <- subset(cml2_seurat, my.demux %in% c("LAMA84-LGAL1-T-","LAMA84-LGAL2-T-", "LAMA84-ctrl1-T-"))
cml2_seurat_tumor_lgal <- subset(cml2_seurat_tumor_lgal, is_nk != "NK")
cml2_seurat_tumor_lgal <- cml2_seurat_tumor_lgal %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml2_seurat_tumor_lgal), nPCs = 11)
cml2_seurat_tumor_lgal@reductions$pca@stdev %>% plot()

cml2_seurat_tumor_lgal <- AddModuleScore(cml2_seurat_tumor_lgal, features = list(c("LGALS9")), name = "LGALS9")

cml2_seurat_tumor_lgal@meta.data %>% 
  ggplot(aes(my.demux,LGALS91,fill=my.demux)) + geom_violin(draw_quantiles = 0.5, adjust = 3) + geom_jitter(size=0.001,alpha=0.3) +
  ggpubr::stat_compare_means(label="p", label.y.npc = 0.9) + ggpubr::rotate_x_text(angle = 45) + 
  #scale_fill_manual(values=c("blue3","green3","green4")) +
  NoLegend() + labs(x="",y="LGALS9") 
ggsave("results/functional/hto2/vln_lgals_ko_tumor.pdf", width = 4, height = 3)



cml2_seurat_tumor_lgal@meta.data %>% group_by(my.demux, lgals9_pos = LGALS91>0) %>% tally() %>% mutate(prop=n/sum(n)) %>% 
  filter(lgals9_pos) %>% 
  
  # ggplot(aes(reorder(my.demux,-prop),prop*100,fill=my.demux)) + geom_bar(stat="identity") +
  ggplot(aes(my.demux,prop*100,fill=my.demux)) + geom_bar(stat="identity") +
  
  # ggpubr::stat_compare_means(label="p") + 
  ggpubr::rotate_x_text(angle = 45) +
  #scale_fill_manual(values=c("blue3","green3","green4")) +
  labs(x="",y="% of LGALS9 pos cells") + NoLegend() 
ggsave("results/functional/hto2/bar_lgals_ko_tumor.pdf", width = 3, height = 4)



table(cml2_seurat$my.demux, cml2_seurat$my.demux.1)
table(cml2_seurat$my.demux, cml2_seurat$my.demux2.1)





## Focus on NK cells; select NK cells based on i) Essi's demultiplex and ii) RNA-profile of NK cells
# cml2_seurat$is_nk <- ifelse(cml2_seurat$RNA_snn_res.0.1 == 3 & grepl("NK", cml2_seurat$my.demux), "NK", "CML")

DimPlot(cml2_seurat, group.by = "is_nk")
cml2_seurat_nk <- subset(cml2_seurat, is_nk == "NK")
cells.to.keep  <- cml2_seurat_nk$barcode[grepl("NK", cml2_seurat_nk$my.demux)]
cml2_seurat_nk <- subset(cml2_seurat_nk, cells = cells.to.keep)

# cml2_seurat_nk <- cml2_seurat_nk %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml2_seurat_nk), nPCs = 20)
cml2_seurat_nk <- cml2_seurat_nk %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml2_seurat_nk), nPCs = 10)
# cml2_seurat_nk <- cml2_seurat_nk %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml2_seurat_nk), nPCs = 7)
cml2_seurat_nk@reductions$pca@stdev %>% plot()

cml2_seurat_nk$nk_presence <- ifelse(grepl(cml2_seurat_nk$my.demux, pattern = "NK"), "NK presence", "No NK")
cml2_seurat_nk$lgals9_ko   <- ifelse(grepl(cml2_seurat_nk$my.demux, pattern = "LGAL"), "LGALS9 KO", "LGALS9 WT")
cml2_seurat_nk$lgals9_ko   <- ifelse(grepl(cml2_seurat_nk$my.demux, pattern = "NK-only"), "No tumor", cml2_seurat_nk$lgals9_ko)
cml2_seurat_nk$expanded    <- ifelse(grepl(cml2_seurat_nk$my.demux, pattern = "NE"), "Non-expanded NK", "Expanded NK")
cml2_seurat_nk$covariate   <- paste0(cml2_seurat_nk$expanded , "\n", cml2_seurat_nk$lgals9_ko)

a <- DimPlot(cml2_seurat_nk, group.by = "my.demux", cols = getPalette(10))
b <- DimPlot(cml2_seurat_nk, group.by = "covariate", cols = getPalette(5))
c <- DimPlot(cml2_seurat_nk, group.by = "expanded", cols = getPalette3(4))
d <- DimPlot(cml2_seurat_nk, group.by = "lgals9_ko", cols = getPalette3(4))

a + b + c + d
ggsave("results/functional/hto2/umap_meta_nk.png", width = 10, height = 7)

saveRDS(cml2_seurat_nk, "results/functional/cml2_seurat_nk.rds")



## Only expanded based on Essi's demultiplexing
cml2_seurat_nk$expanded <- ifelse(grepl(cml2_seurat_nk$my.demux, pattern = "E\\-"), "Expanded NK", "Other")
cml2_seurat_nk$expanded <- ifelse(grepl(cml2_seurat_nk$my.demux, pattern = "NE\\-"), "Non-expanded NK", cml2_seurat_nk$expanded)

cml2_seurat_nk_expanded <- subset(cml2_seurat_nk, expanded == "Expanded NK")
# cml2_seurat_nk_expanded <- cml2_seurat_nk_expanded %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml2_seurat_nk_expanded), nPCs = 12)
# cml2_seurat_nk_expanded <- cml2_seurat_nk_expanded %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml2_seurat_nk_expanded), nPCs = 7)
cml2_seurat_nk_expanded <- cml2_seurat_nk_expanded %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml2_seurat_nk_expanded), nPCs = 18)
cml2_seurat_nk_expanded <- cml2_seurat_nk_expanded %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml2_seurat_nk_expanded), nPCs = 30)
cml2_seurat_nk_expanded <- cml2_seurat_nk_expanded %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml2_seurat_nk_expanded), nPCs = 4)
cml2_seurat_nk_expanded@reductions$pca@stdev %>% plot()

a <- DimPlot(cml2_seurat_nk_expanded, group.by = "my.demux", cols = getPalette(10))
b <- DimPlot(cml2_seurat_nk_expanded, group.by = "covariate", cols = getPalette(5))
c <- DimPlot(cml2_seurat_nk_expanded, group.by = "expanded", cols = getPalette3(4))
d <- DimPlot(cml2_seurat_nk_expanded, group.by = "lgals9_ko", cols = getPalette3(4))

a + b + c + d
ggsave("results/functional/hto2/umap_meta_exp_nk2.png", width = 8, height = 4)

DimPlot(cml2_seurat_nk_expanded, group.by = "Phase", cols = getPalette3(4))
DimPlot(cml2_seurat_nk_expanded, group.by = "RNA_snn_res.1", cols = getPalette3(6), label = T)

Idents(cml2_seurat_nk_expanded) <- cml2_seurat_nk_expanded$RNA_snn_res.1
cml2_seurat_nk_expanded_markers <- FindAllMarkers(cml2_seurat_nk_expanded, only.pos = T)

cml2_seurat_nk_expanded_markers %>% filter(cluster == 5) %>% View


Idents(cml2_seurat_nk_expanded) <- cml2_seurat_nk_expanded$my.demux
DimPlot(cml2_seurat_nk_expanded, cols = c("green4", "blue3", "green3")) 
ggsave("results/functional/hto2/umap_expanded_tumor_nk_cells_lgals9.png", width = 5, height = 3)

Idents(cml2_seurat_nk_expanded) <- cml2_seurat_nk_expanded$my.demux
DimPlot(cml2_seurat_nk_expanded, cols = c("blue3", "green4", "green3")) 
ggsave("results/functional/hto2/umap_expanded_tumor_nk_cells_lgals9.png", width = 5, height = 3)

saveRDS(cml2_seurat_nk_expanded, "results/functional/cml2_seurat_nk_expanded.rds")


Idents(cml2_seurat_nk_expanded) <- cml2_seurat_nk_expanded$lgals9_ko
cml2_seurat_nk_expanded_lgals9 <- FindMarkers(cml2_seurat_nk_expanded, logfc.threshold = 0.15, ident.1 = "LGALS9 KO", ident.2 = "LGALS9 WT", test.use = "t", max.cells.per.ident = 1e3) %>%
  filter(p_val_adj < 0.05) %>% add_rownames(var = "gene")

ggplot(data = cml2_seurat_nk_expanded_lgals9, aes(avg_log2FC, -log10(p_val_adj), label = gene)) + geom_point() + 
  ggrepel::geom_text_repel(data = subset(cml2_seurat_nk_expanded_lgals9, p_val_adj < 1e-5)) +
  ggpubr::theme_pubr() + xlim(values=c(-2,2))
ggsave("results/functional/hto2/volcano_exp_nk_lgals9_ko_vs_wt.png", width = 9, height = 9)

cml2_seurat_nk_expanded_lgals9_hallmark_up <- cml2_seurat_nk_expanded_lgals9 %>% filter(avg_log2FC > 0) %>% getHypergeometric(universe_df = rownames(cml2_seurat_nk_expanded), term_df = hallmark_pathways)
cml2_seurat_nk_expanded_lgals9_hallmark_dn <- cml2_seurat_nk_expanded_lgals9 %>% filter(avg_log2FC < 0) %>% getHypergeometric(universe_df = rownames(cml2_seurat_nk_expanded), term_df = hallmark_pathways)

cml2_seurat_nk_expanded_lgals9_reactome_up <- cml2_seurat_nk_expanded_lgals9 %>% filter(avg_log2FC > 0) %>% getHypergeometric(universe_df = rownames(cml2_seurat_nk_expanded), term_df = canonical_pathways)
cml2_seurat_nk_expanded_lgals9_reactome_dn <- cml2_seurat_nk_expanded_lgals9 %>% filter(avg_log2FC < 0) %>% getHypergeometric(universe_df = rownames(cml2_seurat_nk_expanded), term_df = canonical_pathways)

cml2_seurat_nk_expanded_lgals9_hallmark_up %>% filter(p.adjust < 0.05)
cml2_seurat_nk_expanded_lgals9_hallmark_dn %>% filter(p.adjust < 0.05)

cml2_seurat_nk_expanded_lgals9_reactome_up %>% filter(p.adjust < 0.05)
cml2_seurat_nk_expanded_lgals9_reactome_dn %>% filter(p.adjust < 0.05)



#### Only expanded, LGALS91, ctrl1
cml2_seurat_nk_expanded_lgals91 <- subset(cml2_seurat_nk_expanded, my.demux %in% c("LAMA84-LGAL1-T-E-NK-", "LAMA84-ctrl1-T-E-NK-"))
cml2_seurat_nk_expanded_lgals91 <- cml2_seurat_nk_expanded_lgals91 %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml2_seurat_nk_expanded_lgals91), nPCs = 12)
cml2_seurat_nk_expanded_lgals91@reductions$pca@stdev %>% plot()
Idents(cml2_seurat_nk_expanded_lgals91) <- cml2_seurat_nk_expanded_lgals91$my.demux

DimPlot(cml2_seurat_nk_expanded_lgals91, cols = c("green4", "blue3", "green3")) 
ggsave("results/functional/hto2/umap_expanded_tumor_nk_cells_lgals91.png", width = 5, height = 3)

Idents(cml2_seurat_nk_expanded_lgals91) <- cml2_seurat_nk_expanded_lgals91$lgals9_ko
cml2_seurat_nk_expanded_lgals91_markers <- FindMarkers(cml2_seurat_nk_expanded_lgals91, logfc.threshold = 0.15, ident.1 = "LGALS9 KO", ident.2 = "LGALS9 WT", test.use = "t", max.cells.per.ident = 1e3) %>%
  filter(p_val_adj < 0.05) %>% add_rownames(var = "gene")

ggplot(data = cml2_seurat_nk_expanded_lgals91_markers, aes(avg_log2FC, -log10(p_val_adj), label = gene)) + geom_point() + 
  ggrepel::geom_text_repel(data = subset(cml2_seurat_nk_expanded_lgals91_markers, p_val_adj < 1e-5)) +
  ggpubr::theme_pubr() + xlim(values=c(-2,2))
ggsave("results/functional/hto2/volcano_exp_nk_lgals9_ko_vs_wt_lgals91.png", width = 9, height = 9)

cml2_seurat_nk_expanded_lgals9_hallmark_up <- cml2_seurat_nk_expanded_lgals91_markers %>% filter(avg_log2FC > 0) %>% getHypergeometric(universe_df = rownames(cml2_seurat_nk_expanded), term_df = hallmark_pathways)
cml2_seurat_nk_expanded_lgals9_hallmark_dn <- cml2_seurat_nk_expanded_lgals91_markers %>% filter(avg_log2FC < 0) %>% getHypergeometric(universe_df = rownames(cml2_seurat_nk_expanded), term_df = hallmark_pathways)

cml2_seurat_nk_expanded_lgals9_reactome_up <- cml2_seurat_nk_expanded_lgals91_markers %>% filter(avg_log2FC > 0) %>% getHypergeometric(universe_df = rownames(cml2_seurat_nk_expanded), term_df = canonical_pathways)
cml2_seurat_nk_expanded_lgals9_reactome_dn <- cml2_seurat_nk_expanded_lgals91_markers %>% filter(avg_log2FC < 0) %>% getHypergeometric(universe_df = rownames(cml2_seurat_nk_expanded), term_df = canonical_pathways)

cml2_seurat_nk_expanded_lgals9_hallmark_up %>% filter(p.adjust < 0.05)
cml2_seurat_nk_expanded_lgals9_hallmark_dn %>% filter(p.adjust < 0.05)

cml2_seurat_nk_expanded_lgals9_reactome_up %>% filter(p.adjust < 0.05)
cml2_seurat_nk_expanded_lgals9_reactome_dn %>% filter(p.adjust < 0.05)




#### Only expanded, LGALS92, ctrl1
cml2_seurat_nk_expanded_lgals92 <- subset(cml2_seurat_nk_expanded, my.demux %in% c("LAMA84-LGAL2-T-E-NK-", "LAMA84-ctrl1-T-E-NK-"))
cml2_seurat_nk_expanded_lgals92 <- cml2_seurat_nk_expanded_lgals92 %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml2_seurat_nk_expanded_lgals92), nPCs = 12)
cml2_seurat_nk_expanded_lgals92@reductions$pca@stdev %>% plot()
Idents(cml2_seurat_nk_expanded_lgals92) <- cml2_seurat_nk_expanded_lgals92$my.demux

DimPlot(cml2_seurat_nk_expanded_lgals92, cols = c("green4", "blue3", "green3")) 
ggsave("results/functional/hto2/umap_expanded_tumor_nk_cells_lgals92.png", width = 5, height = 3)

Idents(cml2_seurat_nk_expanded_lgals92) <- cml2_seurat_nk_expanded_lgals92$lgals9_ko
cml2_seurat_nk_expanded_lgals92_markers <- FindMarkers(cml2_seurat_nk_expanded_lgals92, logfc.threshold = 0.15, ident.1 = "LGALS9 KO", ident.2 = "LGALS9 WT", test.use = "t", max.cells.per.ident = 1e3) %>%
  filter(p_val_adj < 0.05) %>% add_rownames(var = "gene")

ggplot(data = cml2_seurat_nk_expanded_lgals92_markers, aes(avg_log2FC, -log10(p_val_adj), label = gene)) + geom_point() + 
  ggrepel::geom_text_repel(data = subset(cml2_seurat_nk_expanded_lgals92_markers, p_val_adj < 1e-5)) +
  ggpubr::theme_pubr() + xlim(values=c(-2,2))
ggsave("results/functional/hto2/volcano_exp_nk_lgals9_ko_vs_wt_lgals92.png", width = 9, height = 9)





## Only non-expanded
cml2_seurat_nk_nonexpanded <- subset(cml2_seurat_nk, expanded != "Expanded NK")
cml2_seurat_nk_nonexpanded <- cml2_seurat_nk_nonexpanded %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml2_seurat_nk_nonexpanded), nPCs = 13)
cml2_seurat_nk_nonexpanded@reductions$pca@stdev %>% plot()
cml2_seurat_nk_nonexpanded <- cml2_seurat_nk_nonexpanded %>% getClustering()
Idents(cml2_seurat_nk_nonexpanded) <- cml2_seurat_nk_nonexpanded$RNA_snn_res.0.1


a <- DimPlot(cml2_seurat_nk_nonexpanded, group.by = "my.demux", cols = getPalette(10))
b <- DimPlot(cml2_seurat_nk_nonexpanded, group.by = "covariate", cols = getPalette(5))
c <- DimPlot(cml2_seurat_nk_nonexpanded, group.by = "expanded", cols = getPalette3(4))
d <- DimPlot(cml2_seurat_nk_nonexpanded, group.by = "lgals9_ko", cols = getPalette3(4))
e <- DimPlot(cml2_seurat_nk_nonexpanded, group.by = "Phase", cols = getPalette3(4))
f <- DimPlot(cml2_seurat_nk_nonexpanded, group.by = "RNA_snn_res.0.1", cols = getPalette3(4))

a + b + c + d + e + f
ggsave("results/functional/hto2/umap_meta_nonexp_nk.png", width = 12, height = 5)

Idents(cml2_seurat_nk_nonexpanded) <- ifelse(cml2_seurat_nk_nonexpanded$covariate == "Non-expanded NK\nNo tumor", "No tumor", "Tumor")
cml2_seurat_nk_nonexpanded_lgals9  <- FindMarkers(cml2_seurat_nk_nonexpanded, logfc.threshold = 0.15, ident.1 = "Tumor", ident.2 = "No tumor", test.use = "t", max.cells.per.ident = 1e3)
cml2_seurat_nk_nonexpanded_lgals9  <- cml2_seurat_nk_nonexpanded_lgals9 %>% filter(p_val_adj < 0.05) %>% add_rownames(var = "gene")

cml2_seurat_nk_nonexpanded_lgals9 %>% filter(avg_log2FC > 0) #%>% View

ggplot(data = cml2_seurat_nk_nonexpanded_lgals9, aes(avg_log2FC, -log10(p_val_adj), label = gene)) + geom_point() + ggrepel::geom_text_repel(data = subset(cml2_seurat_nk_nonexpanded_lgals9, p_val_adj < 1e-15)) +
  ggpubr::theme_pubr()
ggsave("results/functional/hto2/volcano_nonexp_nk_tumor_vs_non_tumor.png", width = 9, height = 9)









## Only non-expanded with knockouts
cml2_seurat_nk_nonexpanded$tumor <- ifelse(cml2_seurat_nk_nonexpanded$covariate == "Non-expanded NK\nNo tumor", "No tumor", "Tumor")
cml2_seurat_nk_nonexpanded2 <- subset(cml2_seurat_nk_nonexpanded, tumor == "Tumor")

cml2_seurat_nk_nonexpanded2 <- cml2_seurat_nk_nonexpanded2 %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml2_seurat_nk_nonexpanded2), nPCs = 13)
cml2_seurat_nk_nonexpanded2@reductions$pca@stdev %>% plot()
cml2_seurat_nk_nonexpanded2 <- cml2_seurat_nk_nonexpanded2 %>% getClustering()
Idents(cml2_seurat_nk_nonexpanded2) <- cml2_seurat_nk_nonexpanded2$RNA_snn_res.0.1

a <- DimPlot(cml2_seurat_nk_nonexpanded2, group.by = "my.demux", cols = getPalette(10))
b <- DimPlot(cml2_seurat_nk_nonexpanded2, group.by = "covariate", cols = getPalette(5))
c <- DimPlot(cml2_seurat_nk_nonexpanded2, group.by = "expanded", cols = getPalette3(4))
d <- DimPlot(cml2_seurat_nk_nonexpanded2, group.by = "lgals9_ko", cols = getPalette3(4))
e <- DimPlot(cml2_seurat_nk_nonexpanded2, group.by = "Phase", cols = getPalette3(4))
f <- DimPlot(cml2_seurat_nk_nonexpanded2, group.by = "RNA_snn_res.0.1", cols = getPalette3(4))

a + b + c + d + e + f
ggsave("results/functional/hto2/umap_meta_nonexp_nk_only_tumor.png", width = 12, height = 5)


DimPlot(cml2_seurat_nk_nonexpanded2, group.by = "my.demux", cols = c("blue3", "green3", "green4")) + NoLegend()
ggsave("results/functional/hto2/umap_nonexp_meta.png", width = 3, height = 3)



Idents(cml2_seurat_nk_nonexpanded) <- cml2_seurat_nk_nonexpanded$lgals9_ko
cml2_seurat_nk_nonexpanded_lgals9_markers  <- FindMarkers(cml2_seurat_nk_nonexpanded, logfc.threshold = 0.15, ident.1 = "LGALS9 KO", ident.2 = "LGALS9 WT", test.use = "t", max.cells.per.ident = 1e3, min.pct = 0.01)
cml2_seurat_nk_nonexpanded_lgals9_markers  <- cml2_seurat_nk_nonexpanded_lgals9_markers %>% filter(p_val_adj < 0.05) %>% add_rownames(var = "gene")
cml2_seurat_nk_nonexpanded_lgals9_markers %>% filter(avg_log2FC > 0) #%>% View

ggplot(data = cml2_seurat_nk_nonexpanded_lgals9_markers, aes(avg_log2FC, -log10(p_val_adj), label = gene)) + 
  geom_point() + ggrepel::geom_text_repel(data = subset(cml2_seurat_nk_nonexpanded_lgals9_markers, p_val_adj < 1e-15)) +
  ggpubr::theme_pubr() + xlim(values = c(-3,3))
ggsave("results/functional/hto2/volcano_nonexp_nk_tumor_vs_non_tumor.png", width = 9, height = 9)



#### Lgals9.1
cml2_seurat_nk_nonexpanded_lgals91 <- subset(cml2_seurat_nk_nonexpanded2, my.demux %in% c("LAMA84-LGAL1-T-NE-NK-", "LAMA84-ctrl1-T-NE-NK-"))
cml2_seurat_nk_nonexpanded_lgals91 <- cml2_seurat_nk_nonexpanded_lgals91 %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml2_seurat_nk_nonexpanded_lgals91), nPCs = 12)
cml2_seurat_nk_nonexpanded_lgals91@reductions$pca@stdev %>% plot()
Idents(cml2_seurat_nk_nonexpanded_lgals91) <- cml2_seurat_nk_nonexpanded_lgals91$my.demux

DimPlot(cml2_seurat_nk_nonexpanded_lgals91, cols = c("blue3", "green4")) + NoLegend()
ggsave("results/functional/hto2/umap_nonexpanded_tumor_nk_cells_lgals91.png", width = 3.5, height = 3)

Idents(cml2_seurat_nk_nonexpanded_lgals91) <- cml2_seurat_nk_nonexpanded_lgals91$lgals9_ko
cml2_seurat_nk_expanded_lgals91_markers <- FindMarkers(cml2_seurat_nk_nonexpanded_lgals91, logfc.threshold = 0.15, ident.1 = "LGALS9 KO", ident.2 = "LGALS9 WT", test.use = "t", max.cells.per.ident = 1e3) %>%
  filter(p_val_adj < 0.05) %>% add_rownames(var = "gene")

ggplot(data = cml2_seurat_nk_expanded_lgals91_markers, aes(avg_log2FC, -log10(p_val_adj), label = gene)) + geom_point() + 
  ggrepel::geom_text_repel(data = subset(cml2_seurat_nk_expanded_lgals91_markers, p_val_adj < 1e-5)) +
  ggpubr::theme_pubr() + xlim(values=c(-3,3))
ggsave("results/functional/hto2/volcano_nonexp_nk_lgals9_ko_vs_wt_lgals91.png", width = 9, height = 9)








#### Lgals9.2
cml2_seurat_nk_nonexpanded_lgals92 <- subset(cml2_seurat_nk_nonexpanded2, my.demux %in% c("LAMA84-LGAL2-T-NE-NK-", "LAMA84-ctrl1-T-NE-NK-"))
cml2_seurat_nk_nonexpanded_lgals92 <- cml2_seurat_nk_nonexpanded_lgals92 %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml2_seurat_nk_nonexpanded_lgals92), nPCs = 12)
cml2_seurat_nk_nonexpanded_lgals92@reductions$pca@stdev %>% plot()
Idents(cml2_seurat_nk_nonexpanded_lgals92) <- cml2_seurat_nk_nonexpanded_lgals92$my.demux

DimPlot(cml2_seurat_nk_nonexpanded_lgals92, cols = c("blue3", "green4")) + NoLegend()
ggsave("results/functional/hto2/umap_nonexpanded_tumor_nk_cells_lgals92.png", width = 3.5, height = 3)

Idents(cml2_seurat_nk_nonexpanded_lgals92) <- cml2_seurat_nk_nonexpanded_lgals92$lgals9_ko
cml2_seurat_nk_expanded_lgals92_markers <- FindMarkers(cml2_seurat_nk_nonexpanded_lgals92, logfc.threshold = 0.15, ident.1 = "LGALS9 KO", ident.2 = "LGALS9 WT", test.use = "t", max.cells.per.ident = 1e3) %>%
  filter(p_val_adj < 0.05) %>% add_rownames(var = "gene")

ggplot(data = cml2_seurat_nk_expanded_lgals92_markers, aes(avg_log2FC, -log10(p_val_adj), label = gene)) + geom_point() + 
  ggrepel::geom_text_repel(data = subset(cml2_seurat_nk_expanded_lgals92_markers, p_val_adj < 1e-5)) +
  ggpubr::theme_pubr() + xlim(values=c(-3,3))
ggsave("results/functional/hto2/volcano_nonexp_nk_lgals9_ko_vs_wt_lgals92.png", width = 9, height = 9)








#### ADD CML3 to CML2 
cml3_seurat_nk_expanded_only <- subset(cml3_seurat_nk_expanded, my.demux == "E-NK-only")
cml2_seurat_nk_expanded_meta <- merge(cml2_seurat_nk_expanded, cml3_seurat_nk_expanded_only)

cml2_seurat_nk_expanded_meta <- cml2_seurat_nk_expanded_meta %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml2_seurat_nk_expanded_meta), nPCs = 11)
cml2_seurat_nk_expanded_meta@reductions$pca@stdev %>% plot()
table(cml2_seurat_nk_expanded_meta$my.demux)


### Satu's edit
cml2_seurat_nk_expanded_meta$satu <- ifelse(cml2_seurat_nk_expanded_meta$my.demux %in% c("E-NK-only"), "No tumor", "other")
cml2_seurat_nk_expanded_meta$satu <- ifelse(cml2_seurat_nk_expanded_meta$my.demux %in% c("LAMA84-ctrl1-T-E-NK-", "LAMA84-ctrl2-T-E-NK-"), "Tumor LGAL wt", cml2_seurat_nk_expanded_meta$satu)
cml2_seurat_nk_expanded_meta$satu <- ifelse(cml2_seurat_nk_expanded_meta$my.demux %in% c("LAMA84-LGAL1-T-E-NK-", "LAMA84-LGAL2-T-E-NK-"), "Tumor LGAL ko", cml2_seurat_nk_expanded_meta$satu)
cml2_seurat_nk_expanded_meta$satu <- factor(cml2_seurat_nk_expanded_meta$satu, levels = c("Tumor LGAL ko", "Tumor LGAL wt", "No tumor"))

Idents(cml2_seurat_nk_expanded_meta) <- cml2_seurat_nk_expanded_meta$satu
DimPlot(cml2_seurat_nk_expanded_meta)

cml2_seurat_nk_expanded_meta <- AddModuleScore(cml2_seurat_nk_expanded_meta, features = list(activation$activated), name = "activation")

dufva_signatures_full <- fread("data/dufva_nk_signatures_full.txt", dec = ",") %>% filter(p_val_adj < 0.05) %>% dplyr::select(cluster, gene)
dufva_signatures_activated_full <- dufva_signatures_full %>% filter(cluster == "Activated (2)") %>% pull(gene)
dufva_signatures_ifna_full      <- dufva_signatures_full %>% filter(cluster == "Type I IFN (3)") %>% pull(gene)


cml2_seurat_nk_expanded_meta <- AddModuleScore(cml2_seurat_nk_expanded_meta, features = list(dufva_signatures_activated_full), name = "activation_full")
cml2_seurat_nk_expanded_meta <- AddModuleScore(cml2_seurat_nk_expanded_meta, features = list(head(dufva_signatures_activated_full,50)), name = "activation_top50")


markers2_1 <- FindMarkers(cml2_seurat_nk_expanded_meta, ident.1 = "Tumor LGAL wt", ident.2 = "No tumor", test.use = "t") %>% add_rownames(var = "gene")# %>% filter(p_val_adj < 0.05)
markers2_2 <- FindMarkers(cml2_seurat_nk_expanded_meta, ident.1 = "Tumor LGAL ko", ident.2 = "No tumor", test.use = "t") %>% add_rownames(var = "gene")# %>% filter(p_val_adj < 0.05)
markers2_3 <- FindMarkers(cml2_seurat_nk_expanded_meta, ident.1 = "Tumor LGAL ko", ident.2 = "Tumor LGAL wt", test.use = "t") %>% add_rownames(var = "gene")# %>% filter(p_val_adj < 0.05)

markers2_1_sigf <- subset(markers2_1, p_val_adj < 0.05)
markers2_2_sigf <- subset(markers2_2, p_val_adj < 0.05)

markers2_1_sigf_exclusive_genes <- markers2_1_sigf$gene[!markers2_1_sigf$gene %in% markers2_2_sigf$gene]
markers2_2_sigf_exclusive_genes <- markers2_2_sigf$gene[!markers2_2_sigf$gene %in% markers2_1_sigf$gene]

full_df2 <- markers2_2 %>% full_join(markers2_1, by = "gene") 

full_df2$avg_log2FC.x <- ifelse(is.na(full_df2$avg_log2FC.x), 0, full_df2$avg_log2FC.x)
full_df2$avg_log2FC.y <- ifelse(is.na(full_df2$avg_log2FC.y), 0, full_df2$avg_log2FC.y)

full_df2$gene_sigf <- ifelse(full_df2$gene %in% markers2_1_sigf_exclusive_genes, "Sigf only in LGAL WT", "Other")
full_df2$gene_sigf <- ifelse(full_df2$gene %in% markers2_2_sigf_exclusive_genes, "Sigf only in LGAL KO", full_df2$gene_sigf)
full_df2$activation <- ifelse(full_df2$gene %in% dufva_signatures_activated_full, "activation", "other")
full_df2$ifna <- ifelse(full_df2$gene %in% dufva_signatures_ifna_full, "ifna", "other")

full_df1$ifna <- ifelse(full_df1$gene %in% dufva_signatures_ifna_full, "ifna", "other")
full_df2$ifna <- ifelse(full_df2$gene %in% dufva_signatures_ifna_full, "ifna", "other")
full_df3$ifna <- ifelse(full_df3$gene %in% dufva_signatures_ifna_full, "ifna", "other")

ggplot() +
  geom_point(data = full_df2, aes(avg_log2FC.x, avg_log2FC.y, color = gene_sigf, label = gene, shape = activation)) + 
  # ggrepel::geom_text_repel(data = subset(full_df2, gene_sigf != "Other" & activation == "activation"), aes(avg_log2FC.x, avg_log2FC.y,  label = gene, shape = activation), font.face = "italic") + 
  ggrepel::geom_text_repel(data = subset(full_df2, gene_sigf == "Sigf only in LGAL KO" & activation == "activation"), aes(avg_log2FC.x, avg_log2FC.y,  label = gene, shape = activation), fontface = "italic") + 
  
  geom_abline() +
  xlim(values = c(0,4)) + ylim(values = c(0,4)) + 
  scale_shape_manual(values = c(17,16)) + scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = c("gray90", "red", "blue")) + add_guide + labs(x = "avg log2FC\nin NKs against LGAL KO\nand NK alone", y = "avg log2FC\nin NKs against LGAL WT\nand NK alone", color = "", shape = "associated\nwith NK activation")
ggsave("results/functional/scatter_cml2_LGAL_ko_LGAL_wt.png", width = 7, height = 5)

fwrite(full_df2, "results/functional/hto3/scatter_cml2_LGAL_ko_LGAL_wt.png", sep = "\t", quote = F, row.names = F)
