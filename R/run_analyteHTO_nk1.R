
cml1_seurat <- cml_seurat_new_filt1

### Show that KO is working?
cml1_seurat_tumor <- subset(cml1_seurat, is_nk != "NK")
cells.to.keep     <- cml1_seurat_tumor@meta.data %>% filter(!my.demux %in% c("E-NK-only", "NE-NK-only")) %>% pull(barcode)
cml1_seurat_tumor <- subset(cml1_seurat_tumor, cells = cells.to.keep)

cml1_seurat_tumor <- cml1_seurat_tumor %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml1_seurat_tumor), nPCs = 14)
cml1_seurat_tumor@reductions$pca@stdev %>% plot()
cml1_seurat_tumor <- AddModuleScore(cml1_seurat_tumor, features = list(c("PVR")), name = "PVR")
cml1_seurat_tumor <- AddModuleScore(cml1_seurat_tumor, features = list(c("PVR")), name = "PVR")
cml1_seurat_tumor <- AddModuleScore(cml1_seurat_tumor, features = list(c("PVR")), name = "PVR")

cml1_seurat_tumor@meta.data %>% 
  ggplot(aes(my.demux,PVR1,fill=my.demux)) + geom_violin(draw_quantiles = 0.5, adjust = 3) + geom_jitter(size=0.001,alpha=0.3) +
  ggpubr::stat_compare_means(label="p", label.y.npc = 0.9) + ggpubr::rotate_x_text(angle = 45) + 
  #scale_fill_manual(values=c("blue3","green3","green4")) +
  NoLegend() + labs(x="",y="PVR") 
ggsave("results/functional/hto1/vln_pvr_ko_total.pdf", width = 4, height = 3)

cml1_seurat_tumor@meta.data %>% group_by(my.demux, pvr_pos = PVR1>0) %>% tally() %>% mutate(prop=n/sum(n)) %>% 
  filter(pvr_pos) %>% 
  
  # ggplot(aes(reorder(my.demux,-prop),prop*100,fill=my.demux)) + geom_bar(stat="identity") +
  ggplot(aes(my.demux,prop*100,fill=my.demux)) + geom_bar(stat="identity") +
  
  # ggpubr::stat_compare_means(label="p") + 
  ggpubr::rotate_x_text(angle = 45) +
  #scale_fill_manual(values=c("blue3","green3","green4")) +
  labs(x="",y="% of PVR pos cells") + NoLegend() 
ggsave("results/functional/hto1/bar_pvr_ko_total.pdf", width = 3, height = 4)





cml1_seurat_tumor_pvr <- subset(cml1_seurat, my.demux %in% c("K562-PVR1-T-","K562-PVR2-T-", "K562-ctrl1-T-", "K562-ctrl2-T-"))
cml1_seurat_tumor_pvr <- subset(cml1_seurat_tumor_pvr, is_nk != "NK")
cml1_seurat_tumor_pvr <- cml1_seurat_tumor_pvr %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml1_seurat_tumor_pvr), nPCs = 11)
cml1_seurat_tumor_pvr@reductions$pca@stdev %>% plot()

cml1_seurat_tumor_pvr <- AddModuleScore(cml1_seurat_tumor_pvr, features = list(c("PVR")), name = "PVR")

cml1_seurat_tumor_pvr@meta.data %>% 
  ggplot(aes(my.demux,PVR1,fill=my.demux)) + geom_violin(draw_quantiles = 0.5, adjust = 3) + geom_jitter(size=0.001,alpha=0.3) +
  ggpubr::stat_compare_means(label="p", label.y.npc = 0.9) + ggpubr::rotate_x_text(angle = 45) + 
  scale_fill_manual(values=c("blue3", "dodgerblue", "green3","green4")) +
  NoLegend() + labs(x="",y="PVR") 
ggsave("results/functional/hto1/vln_pvr_ko.pdf", width = 2, height = 3)

cml1_seurat_tumor_pvr@meta.data %>% group_by(my.demux, pvr_pos = PVR1>0) %>% tally() %>% mutate(prop=n/sum(n)) %>% 
  filter(pvr_pos) %>% 
  
  ggplot(aes(reorder(my.demux,-prop),prop*100,fill=my.demux)) + geom_bar(stat="identity") +
  ggpubr::stat_compare_means(label="p") + ggpubr::rotate_x_text(angle = 45) +
  labs(x="",y="% of PVR pos cells") + NoLegend() + 
  scale_fill_manual(values=c("blue3", "dodgerblue", "green3","green4"))
ggsave("results/functional/hto1/bar_pvr_ko.pdf", width = 2, height = 3)







## Focus on NK cells; select NK cells based on i) Essi's demultiplex and ii) RNA-profile of NK cells
cml1_seurat$is_nk <- ifelse(cml1_seurat$RNA_snn_res.0.1 == 2 & grepl("NK", cml1_seurat$my.demux), "NK", "CML")

cml1_seurat_nk <- subset(cml1_seurat, is_nk == "NK")
cml1_seurat_nk <- cml1_seurat_nk %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml1_seurat_nk), nPCs = 12)
cml1_seurat_nk@reductions$pca@stdev %>% plot()
cml1_seurat_nk <- cml1_seurat_nk %>% getClustering()

cml1_seurat_nk$nk_presence <- ifelse(grepl(cml1_seurat_nk$my.demux, pattern = "NK"), "NK presence", "No NK")
cml1_seurat_nk$pvr_ko      <- ifelse(grepl(cml1_seurat_nk$my.demux, pattern = "PVR"), "PVR KO", "PVR WT")
cml1_seurat_nk$tumor       <- ifelse(grepl(cml1_seurat_nk$my.demux, pattern = "NK-only"), "No tumor", "Tumor")
cml1_seurat_nk$expanded    <- ifelse(grepl(cml1_seurat_nk$my.demux, pattern = "NE"), "Non-expanded NK", "Expanded NK")
cml1_seurat_nk$covariate   <- paste0(cml1_seurat_nk$expanded , "\n", cml1_seurat_nk$pvr_ko)

a <- DimPlot(cml1_seurat_nk, group.by = "my.demux", cols = getPalette(10))
b <- DimPlot(cml1_seurat_nk, group.by = "covariate", cols = getPalette(5))
c <- DimPlot(cml1_seurat_nk, group.by = "expanded", cols = getPalette3(4))
d <- DimPlot(cml1_seurat_nk, group.by = "pvr_ko", cols = getPalette3(4))
e <- DimPlot(cml1_seurat_nk, group.by = "tumor", cols = getPalette3(3))

a + b + c + d + e
ggsave("results/functional/hto1/umap_meta_nk.png", width = 10, height = 7)

saveRDS(cml1_seurat_nk, "results/functional/cml1_seurat_nk.rds")

## Only expanded based on Essi's demultiplexing
cml1_seurat_nk$expanded <- ifelse(grepl(cml1_seurat_nk$my.demux, pattern = "NE"), "Non-expanded NK", "Expanded NK")
cml1_seurat_nk_expanded <- subset(cml1_seurat_nk, expanded == "Expanded NK")
cml1_seurat_nk_expanded <- cml1_seurat_nk_expanded %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml1_seurat_nk_expanded), nPCs = 7)
cml1_seurat_nk_expanded@reductions$pca@stdev %>% plot()

a <- DimPlot(cml1_seurat_nk_expanded, group.by = "my.demux", cols = getPalette(10))
b <- DimPlot(cml1_seurat_nk_expanded, group.by = "covariate", cols = getPalette(5))
c <- DimPlot(cml1_seurat_nk_expanded, group.by = "expanded", cols = getPalette3(4))
d <- DimPlot(cml1_seurat_nk_expanded, group.by = "pvr_ko", cols = getPalette3(4))
e <- DimPlot(cml1_seurat_nk_expanded, group.by = "tumor", cols = getPalette3(3))

a + b + c + d + e
ggsave("results/functional/hto1/umap_meta_exp_nk.png", width = 8, height = 4)

saveRDS(cml1_seurat_nk_expanded, "results/functional/cml1_seurat_nk_expanded.rds")





## Tumor vs no tumor, expanded
cml1_seurat_nk_expanded_tumor_no_tumor <- subset(cml1_seurat_nk, my.demux %in% c("K562-ctrl1-T-E-NK-", "K562-ctrl2-T-E-NK-", "E-NK-only"))
cml1_seurat_nk_expanded_tumor_no_tumor <- cml1_seurat_nk_expanded_tumor_no_tumor %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml1_seurat_nk_expanded_tumor_no_tumor), nPCs = 12)
cml1_seurat_nk_expanded_tumor_no_tumor@reductions$pca@stdev %>% plot()
Idents(cml1_seurat_nk_expanded_tumor_no_tumor) <- cml1_seurat_nk_expanded_tumor_no_tumor$tumor

saveRDS(cml1_seurat_nk_expanded_tumor_no_tumor, "results/functional/cml1_seurat_nk_expanded_tumor_no_tumor.rds")

cml1_seurat_nk_expanded_tumor_no_tumor$tumor_presence <- ifelse(cml1_seurat_nk_expanded_tumor_no_tumor$my.demux == "E-NK-only", "Without tumor", "With tumor")

DimPlot(cml1_seurat_nk_expanded_tumor_no_tumor, group.by = "my.demux", cols = c("red3", "blue3", "dodgerblue"))
ggsave("results/functional/hto1/umap_meta_exp_nk_tumor_only.png", width = 5, height = 3)

DimPlot(cml1_seurat_nk_expanded_tumor_no_tumor, group.by = "tumor_presence", cols = c("darkgoldenrod4", "darkslategray4")) + theme_void(base_size = 17) + NoLegend() + labs(title = "")
ggsave("results/functional/umap_1_exp_nk_tumor_only.png", width = 2, height = 2)

DimPlot(cml1_seurat_nk_expanded_tumor_no_tumor, group.by = "tumor_presence", cols = c("darkgoldenrod4", "darkslategray4")) + theme_void(base_size = 17) + labs(title = "") + add_guide
ggsave("results/functional/umap_1_exp_nk_tumor_only2.png", width = 2, height = 2)

FeaturePlot(cml1_seurat_nk_expanded_tumor_no_tumor, features = c("TNFRSF18", "CRTAM", "GZMA", "NKG7"), cols = c("gray90", "red3"))
ggsave("results/functional/hto1/umap_meta_exp_nk_tumor_only_features.png", width = 6, height = 4)


cml1_seurat_nk_expanded_tumor_no_tumor_markers <- FindMarkers(cml1_seurat_nk_expanded, logfc.threshold = 0.15, ident.1 = "Tumor", ident.2 = "No tumor", 
                                                              test.use = "t", max.cells.per.ident = 1e3) %>%
  filter(p_val_adj < 0.05) %>% add_rownames(var = "gene")
fwrite(cml1_seurat_nk_expanded_tumor_no_tumor_markers, "results/functional/cml1_seurat_nk_expanded_tumor_no_tumor_markers.txt", sep = "\t", quote = F, row.names = F)

ggplot() +
  geom_point(data = cml1_seurat_nk_expanded_tumor_no_tumor_markers, aes(avg_log2FC, -log10(p_val), color = avg_log2FC > 0)) + 
  ggrepel::geom_text_repel(data = subset(cml1_seurat_nk_expanded_tumor_no_tumor_markers, -log10(cml1_seurat_nk_expanded_tumor_no_tumor_markers$p_val) > 20), aes(avg_log2FC, -log10(p_val), label = gene), fontface = "italic") + 
  theme_classic(base_size = 17) + scale_color_manual(values = c("darkslategray4", "darkgoldenrod4")) + NoLegend() + xlim(values = c(-3,3))
ggsave("results/functional/volcano_1_exp_nk_tumor_only.png", width = 5.5, height = 4.5)

FeaturePlot(cml1_seurat_nk_expanded_tumor_no_tumor, features = c("TNFRSF18", "TNFRSF4", "TNFRSF9", "CRTAM", "HAVCR2", "TIGIT"), cols = c("gray90", "red3"), ncol = 3) & theme_void() & NoLegend() & theme(title = element_text(face = "italic")) #& labs(title = "")
ggsave("results/functional/umap_1_exp_nk_tumor_only_features.png", width = 4, height = 3)



ggplot(data = cml1_seurat_nk_expanded_tumor_no_tumor_markers, aes(avg_log2FC, -log10(p_val_adj), label = gene)) + geom_point() + 
  ggrepel::geom_text_repel(data = subset(cml1_seurat_nk_expanded_tumor_no_tumor_markers, p_val_adj < 1e-10)) +
  ggpubr::theme_pubr() + xlim(values=c(-3.5,3.5))
ggsave("results/functional/hto1/volcano_exp_nk_tumor_vs_not_tumor.png", width = 9, height = 9)







##### Only expanded
cml1_seurat_nk_expanded <- subset(cml1_seurat_nk, expanded == "Expanded NK")
cml1_seurat_nk_expanded <- cml1_seurat_nk_expanded %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml1_seurat_nk_expanded), nPCs = 11)
cml1_seurat_nk_expanded@reductions$pca@stdev %>% plot()
Idents(cml1_seurat_nk_expanded) <- cml1_seurat_nk_expanded$my.demux

DimPlot(cml1_seurat_nk_expanded, group.by = "my.demux", cols = c("red3", "blue3", "dodgerblue", "green3", "green4")) + NoLegend()
ggsave("results/functional/hto1/umap_meta_exp_nk.png", width = 3, height = 3)

DimPlot(cml1_seurat_nk_expanded, group.by = "my.demux", cols = c("red3", "blue3", "dodgerblue", "green3", "green4")) #+ NoLegend()

FeaturePlot(cml1_seurat_nk_expanded, features = c("TNFRSF18", "CRTAM", "GZMA", "NKG7"), cols = c("gray90", "red3"))
ggsave("results/functional/hto1/umap_meta_exp_nk_features.png", width = 6, height = 4)

FeaturePlot(cml1_seurat_nk_expanded, features = c("TIGIT", "CD226", "CD96", "KIR2DL5A", "CD112R"), cols = c("gray90", "red3"))
cml1_seurat_nk_expanded_markers %>% filter(gene %in% c("TIGIT", "CD226", "CD96", "KIR2DL5A", "CD112R"))

Idents(cml1_seurat_nk_expanded) <- cml1_seurat_nk_expanded$tumor
cml1_seurat_nk_expanded_markers <- FindMarkers(cml1_seurat_nk_expanded, logfc.threshold = 0.15, ident.1 = "Tumor", ident.2 = "No tumor", 
                                                              test.use = "t", max.cells.per.ident = 1e3) %>%
  filter(p_val_adj < 0.05) %>% add_rownames(var = "gene")
fwrite(cml1_seurat_nk_expanded_tumor_no_tumor_markers, "results/functional/hto1/cml1_seurat_nk_expanded_tumor_no_tumor_markers.txt", sep = "\t", quote = F, row.names = F)

ggplot(data = cml1_seurat_nk_expanded_tumor_no_tumor_markers, aes(avg_log2FC, -log10(p_val_adj), label = gene)) + geom_point() + 
  ggrepel::geom_text_repel(data = subset(cml1_seurat_nk_expanded_tumor_no_tumor_markers, p_val_adj < 1e-10)) +
  ggpubr::theme_pubr() + xlim(values=c(-3.5,3.5))
ggsave("results/functional/hto1/volcano_exp_nk_tumor_vs_not_tumor.png", width = 9, height = 9)





### Satu's edit
cml1_seurat_nk_expanded$satu <- ifelse(cml1_seurat_nk_expanded$my.demux %in% c("E-NK-only"), "No tumor", "other")
cml1_seurat_nk_expanded$satu <- ifelse(cml1_seurat_nk_expanded$my.demux %in% c("K562-ctrl1-T-E-NK-", "K562-ctrl2-T-E-NK-"), "Tumor PVR wt", cml1_seurat_nk_expanded$satu)
cml1_seurat_nk_expanded$satu <- ifelse(cml1_seurat_nk_expanded$my.demux %in% c("K562-PVR1-T-E-NK-", "K562-PVR2-T-E-NK-"), "Tumor PVR ko", cml1_seurat_nk_expanded$satu)
cml1_seurat_nk_expanded$satu <- factor(cml1_seurat_nk_expanded$satu, levels = c("Tumor PVR ko", "Tumor PVR wt", "No tumor"))

cml1_seurat_nk_expanded <- AddModuleScore(cml1_seurat_nk_expanded, features = list(activation$activated), name = "activation")

dufva_signatures_full <- fread("data/dufva_nk_signatures_full.txt", dec = ",") %>% filter(p_val_adj < 0.05) %>% dplyr::select(cluster, gene)
dufva_signatures_activated_full <- dufva_signatures_full %>% filter(cluster == "Activated (2)") %>% pull(gene)

markers1_2 %>% filter(!gene %in% markers1_1$gene) %>% filter(!gene %in% dufva_signatures_full$gene) 
markers1_2 %>% filter(!gene %in% markers1_1$gene) %>% arrange(-(avg_log2FC)) %>% filter(avg_log2FC > 0) %>% filter(!gene %in% unwanted_genes) %>% filter(gene %in% dufva_signatures_activated_full) %>% pull(gene)

cml1_seurat_nk_expanded <- AddModuleScore(cml1_seurat_nk_expanded, features = list(dufva_signatures_activated_full), name = "activation_full")
cml1_seurat_nk_expanded <- AddModuleScore(cml1_seurat_nk_expanded, features = list(head(dufva_signatures_activated_full,50)), name = "activation_top50")

markers1_1 <- FindMarkers(cml1_seurat_nk_expanded, ident.1 = "Tumor PVR wt", ident.2 = "No tumor", test.use = "t") %>% add_rownames(var = "gene")# %>% filter(p_val_adj < 0.05)
markers1_2 <- FindMarkers(cml1_seurat_nk_expanded, ident.1 = "Tumor PVR ko", ident.2 = "No tumor", test.use = "t") %>% add_rownames(var = "gene")# %>% filter(p_val_adj < 0.05)
markers1_3 <- FindMarkers(cml1_seurat_nk_expanded, ident.1 = "Tumor PVR ko", ident.2 = "Tumor PVR wt", test.use = "t") %>% add_rownames(var = "gene")# %>% filter(p_val_adj < 0.05)

markers1_1_sigf <- subset(markers1_1, p_val_adj < 0.05)
markers1_2_sigf <- subset(markers1_2, p_val_adj < 0.05)

markers1_1_sigf_exclusive_genes <- markers1_1_sigf$gene[!markers1_1_sigf$gene %in% markers1_2_sigf$gene]
markers1_2_sigf_exclusive_genes <- markers1_2_sigf$gene[!markers1_2_sigf$gene %in% markers1_1_sigf$gene]

full_df1 <- markers1_2 %>% full_join(markers1_1, by = "gene") 

full_df1$avg_log2FC.x <- ifelse(is.na(full_df1$avg_log2FC.x), 0, full_df1$avg_log2FC.x)
full_df1$avg_log2FC.y <- ifelse(is.na(full_df1$avg_log2FC.y), 0, full_df1$avg_log2FC.y)

full_df1$gene_sigf <- ifelse(full_df1$gene %in% markers1_1_sigf_exclusive_genes, "Sigf only in PVR WT", "Other")
full_df1$gene_sigf <- ifelse(full_df1$gene %in% markers1_2_sigf_exclusive_genes, "Sigf only in PVR KO", full_df1$gene_sigf)
full_df1$activation <- ifelse(full_df1$gene %in% dufva_signatures_activated_full, "activation", "other")


ggplot() +
  geom_point(data = full_df1, aes(avg_log2FC.x, avg_log2FC.y, color = gene_sigf, label = gene, shape = activation)) + 
  # ggrepel::geom_text_repel(data = subset(full_df1, gene_sigf != "Other" & activation == "activation"), aes(avg_log2FC.x, avg_log2FC.y,  label = gene, shape = activation), font.face = "italic") + 
  ggrepel::geom_text_repel(data = subset(full_df1, gene_sigf == "Sigf only in PVR KO" & activation == "activation"), aes(avg_log2FC.x, avg_log2FC.y,  label = gene, shape = activation), fontface = "italic") + 
  
  geom_abline() +
  xlim(values = c(0,4)) + ylim(values = c(0,4)) + scale_shape_manual(values = c(17,16)) + scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = c("gray90", "red", "blue")) + add_guide + labs(x = "avg log2FC\nin NKs against PVR KO\nand NK alone", y = "avg log2FC\nin NKs against PVR WT\nand NK alone", color = "", shape = "associated\nwith NK activation")
ggsave("results/functional/scatter_cml1_pvr_ko_pvr_wt.png", width = 7, height = 5)

fwrite(full_df1, "results/functional/scatter_cml1_pvr_ko_pvr_wt.txt", sep = "\t", quote = F, row.names = F)
fwrite(full_df2, "results/functional/scatter_cml2_lgal_ko_pvr_wt.txt", sep = "\t", quote = F, row.names = F)
fwrite(full_df3, "results/functional/scatter_cml3_pvr_ko_pvr_wt.txt", sep = "\t", quote = F, row.names = F)

a <- full_df1 %>% filter(gene_sigf == "Sigf only in PVR KO") 
b <- full_df2 %>% filter(gene_sigf == "Sigf only in LGAL KO") 
c <- full_df3 %>% filter(gene_sigf == "Sigf only in PVR KO") 

View(a)
View(b)
View(c)

a %>% filter(gene %in% b$gene) %>% filter(gene %in% c$gene)

VlnPlot(cml1_seurat_nk_expanded, features = c("DUSP2", "NFKBIA", "LAG3"), cols = c("red3", "blue3", "gray90"))
ggsave("results/functional/vln_cml1_pvr_ko_pvr_wt.png", width = 4, height = 3)

VlnPlot(cml2_seurat_nk_expanded_meta, features = c("DUSP2", "NFKB1", "TRAF1"), cols = c("red3", "blue3", "gray90"))
ggsave("results/functional/vln_cml2_pvr_ko_pvr_wt.png", width = 4, height = 3)

VlnPlot(cml3_seurat_nk_expanded, features = c("MYC", "PIM3", "EGR2"), cols = c("red3", "blue3", "gray90"))
ggsave("results/functional/vln_cml3_pvr_ko_pvr_wt.png", width = 4, height = 3)


VlnPlot(cml1_seurat_nk_expanded, features = c("HAVCR2", "TIGIT"), cols = c("red3", "blue3", "gray90"))
ggsave("results/functional/vln_cml1_satu.png", width = 4, height = 3)
VlnPlot(cml2_seurat_nk_expanded, features = c("HAVCR2", "TIGIT"), cols = c("red3", "blue3", "gray90"))
ggsave("results/functional/vln_cml2_satu.png", width = 4, height = 3)
VlnPlot(cml3_seurat_nk_expanded, features = c("HAVCR2", "TIGIT"), cols = c("red3", "blue3", "gray90"))
ggsave("results/functional/vln_cml3_satu.png", width = 4, height = 3)


## Only expanded, with tumor
cml1_seurat_nk_expanded_tumor <- subset(cml1_seurat_nk, tumor == "Tumor")
cml1_seurat_nk_expanded_tumor <- subset(cml1_seurat_nk_expanded_tumor, expanded == "Expanded NK")
cml1_seurat_nk_expanded_tumor <- cml1_seurat_nk_expanded_tumor %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml1_seurat_nk_expanded_tumor), nPCs = 11)
cml1_seurat_nk_expanded_tumor@reductions$pca@stdev %>% plot()

Idents(cml1_seurat_nk_expanded_tumor) <- cml1_seurat_nk_expanded_tumor$my.demux
DimPlot(cml1_seurat_nk_expanded_tumor, group.by = "my.demux", cols = c("blue3", "dodgerblue", "green3", "green4")) + NoLegend()
ggsave("results/functional/hto1/umap_meta_exp_nk_tumor.png", width = 3, height = 3)





## Only expanded, with tumor; PVR1.1 guide and ctrl1 guide
cml1_seurat_nk_expanded_tumor_pvr1.1_ctrl1.1 <- subset(cml1_seurat_nk_expanded_tumor, my.demux %in% c("K562-ctrl1-T-E-NK-", "K562-PVR1-T-E-NK-"))
cml1_seurat_nk_expanded_tumor_pvr1.1_ctrl1.1 <- cml1_seurat_nk_expanded_tumor_pvr1.1_ctrl1.1 %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml1_seurat_nk_expanded_tumor_pvr1.1_ctrl1.1), nPCs = 22)
cml1_seurat_nk_expanded_tumor_pvr1.1_ctrl1.1@reductions$pca@stdev %>% plot()
Idents(cml1_seurat_nk_expanded_tumor_pvr1.1_ctrl1.1) <- cml1_seurat_nk_expanded_tumor_pvr1.1_ctrl1.1$pvr_ko

Idents(cml1_seurat_nk_expanded_tumor_pvr1.1_ctrl1.1) <- cml1_seurat_nk_expanded_tumor_pvr1.1_ctrl1.1$my.demux
DimPlot(cml1_seurat_nk_expanded_tumor_pvr1.1_ctrl1.1, group.by = "my.demux", cols = c("blue3", "green3"))+ NoLegend()
ggsave("results/functional/hto1/umap_meta_exp_pvr1_ctrl1.png", width = 3, height = 3)



Idents(cml1_seurat_nk_expanded_tumor_pvr1.1_ctrl1.1) <- cml1_seurat_nk_expanded_tumor_pvr1.1_ctrl1.1$PVR_ko

cml1_seurat_nk_expanded_tumor_pvr1.1_ctrl1.1_markers <- FindMarkers(cml1_seurat_nk_expanded_tumor_pvr1.1_ctrl1.1, logfc.threshold = 0.15, ident.1 = "PVR KO", ident.2 = "PVR WT", 
                                                            test.use = "t", max.cells.per.ident = 1e3) %>%
  filter(p_val_adj < 0.05) %>% add_rownames(var = "gene")

ggplot(data = cml1_seurat_nk_expanded_tumor_pvr1.1_ctrl1.1_markers, aes(avg_log2FC, -log10(p_val_adj), label = gene)) + geom_point() + 
  ggrepel::geom_text_repel(data = subset(cml1_seurat_nk_expanded_tumor_pvr1.1_ctrl1.1_markers, p_val_adj < 1e-10)) +
  ggpubr::theme_pubr()
ggsave("results/functional/hto1/volcano_exp_nk_pvr_ko_pvr_wt.png", width = 9, height = 9)








## Only expanded, with tumor; PVR1.1 guide and ctrl2 guide
cml1_seurat_nk_expanded_tumor_pvr1.1_ctrl1.2 <- subset(cml1_seurat_nk_expanded_tumor, my.demux %in% c("K562-ctrl2-T-E-NK-", "K562-PVR1-T-E-NK-"))
cml1_seurat_nk_expanded_tumor_pvr1.1_ctrl1.2 <- cml1_seurat_nk_expanded_tumor_pvr1.1_ctrl1.2 %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml1_seurat_nk_expanded_tumor_pvr1.1_ctrl1.2), nPCs = 22)
cml1_seurat_nk_expanded_tumor_pvr1.1_ctrl1.2@reductions$pca@stdev %>% plot()
Idents(cml1_seurat_nk_expanded_tumor_pvr1.1_ctrl1.2) <- cml1_seurat_nk_expanded_tumor_pvr1.1_ctrl1.2$pvr_ko

Idents(cml1_seurat_nk_expanded_tumor_pvr1.1_ctrl1.2) <- cml1_seurat_nk_expanded_tumor_pvr1.1_ctrl1.2$my.demux
DimPlot(cml1_seurat_nk_expanded_tumor_pvr1.1_ctrl1.2, group.by = "my.demux", cols = c("dodgerblue", "green3")) + NoLegend()
ggsave("results/functional/hto1/umap_meta_exp_pvr1_ctrl2.png", width = 3, height = 3)

Idents(cml1_seurat_nk_expanded_tumor_pvr1.1_ctrl1.2) <- cml1_seurat_nk_expanded_tumor_pvr1.1_ctrl1.2$pvr_ko
cml1_seurat_nk_expanded_tumor_pvr1.1_ctrl1.2_markers <- FindMarkers(cml1_seurat_nk_expanded_tumor_pvr1.1_ctrl1.2, logfc.threshold = 0.15, ident.1 = "PVR KO", ident.2 = "PVR WT", 
                                                                    test.use = "t", max.cells.per.ident = 1e3) %>%
  filter(p_val_adj < 0.05) %>% add_rownames(var = "gene")

ggplot(data = cml1_seurat_nk_expanded_tumor_pvr1.1_ctrl1.2_markers, aes(avg_log2FC, -log10(p_val_adj), label = gene)) + geom_point() + 
  ggrepel::geom_text_repel(data = subset(cml1_seurat_nk_expanded_tumor_pvr1.1_ctrl1.1_markers, p_val_adj < 1e-10)) +
  ggpubr::theme_pubr()
ggsave("results/functional/hto1/volcano_exp_nk_pvr1.1_ko_ctrl1.2.png", width = 9, height = 9)







## Only expanded, with tumor; PVR2 guide and ctrl1 guide
cml1_seurat_nk_expanded_tumor_pvr1.2_ctrl1.1 <- subset(cml1_seurat_nk_expanded_tumor, my.demux %in% c("K562-ctrl1-T-E-NK-", "K562-PVR2-T-E-NK-"))
cml1_seurat_nk_expanded_tumor_pvr1.2_ctrl1.1 <- cml1_seurat_nk_expanded_tumor_pvr1.2_ctrl1.1 %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml1_seurat_nk_expanded_tumor_pvr1.2_ctrl1.1), nPCs = 24)
cml1_seurat_nk_expanded_tumor_pvr1.2_ctrl1.1@reductions$pca@stdev %>% plot()
Idents(cml1_seurat_nk_expanded_tumor_pvr1.2_ctrl1.1) <- cml1_seurat_nk_expanded_tumor_pvr1.2_ctrl1.1$pvr_ko

Idents(cml1_seurat_nk_expanded_tumor_pvr1.2_ctrl1.1) <- cml1_seurat_nk_expanded_tumor_pvr1.2_ctrl1.1$my.demux
DimPlot(cml1_seurat_nk_expanded_tumor_pvr1.2_ctrl1.1, group.by = "my.demux", cols = c("blue3", "green4")) + NoLegend()
ggsave("results/functional/hto1/umap_meta_exp_pvr2_ctrl1.png", width = 3, height = 3)



Idents(cml1_seurat_nk_expanded_tumor_pvr1.2_ctrl1.1) <- cml1_seurat_nk_expanded_tumor_pvr1.2_ctrl1.1$pvr_ko
cml1_seurat_nk_expanded_tumor_pvr1.2_ctrl1.1_markers <- FindMarkers(cml1_seurat_nk_expanded_tumor_pvr1.2_ctrl1.1, logfc.threshold = 0.15, 
                                                                    ident.1 = "PVR KO", ident.2 = "PVR WT", 
                                                                    test.use = "t", max.cells.per.ident = 1e3) %>%
  filter(p_val_adj < 0.05) %>% add_rownames(var = "gene")

ggplot(data = cml1_seurat_nk_expanded_tumor_pvr1.2_ctrl1.1_markers, aes(avg_log2FC, -log10(p_val_adj), label = gene)) + geom_point() + 
  ggrepel::geom_text_repel(data = subset(cml1_seurat_nk_expanded_tumor_pvr1.2_ctrl1.1_markers, p_val_adj < 1e-10)) +
  ggpubr::theme_pubr()
ggsave("results/functional/hto1/volcano_exp_nk_pvr1.2_ko_ctrl1.1.png", width = 9, height = 9)







## Only expanded, with tumor; PVR1.2 guide and ctrl2 guide
cml1_seurat_nk_expanded_tumor_pvr1.2_ctrl1.2 <- subset(cml1_seurat_nk_expanded_tumor, my.demux %in% c("K562-ctrl2-T-E-NK-", "K562-PVR2-T-E-NK-"))
cml1_seurat_nk_expanded_tumor_pvr1.2_ctrl1.2 <- cml1_seurat_nk_expanded_tumor_pvr1.2_ctrl1.2 %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml1_seurat_nk_expanded_tumor_pvr1.2_ctrl1.2), nPCs = 24)
cml1_seurat_nk_expanded_tumor_pvr1.2_ctrl1.2@reductions$pca@stdev %>% plot()
Idents(cml1_seurat_nk_expanded_tumor_pvr1.2_ctrl1.2) <- cml1_seurat_nk_expanded_tumor_pvr1.2_ctrl1.2$pvr_ko

Idents(cml1_seurat_nk_expanded_tumor_pvr1.2_ctrl1.2) <- cml1_seurat_nk_expanded_tumor_pvr1.2_ctrl1.2$my.demux
DimPlot(cml1_seurat_nk_expanded_tumor_pvr1.2_ctrl1.2, group.by = "my.demux", cols = c("dodgerblue", "green3")) + NoLegend()
ggsave("results/functional/hto1/umap_meta_exp_pvr2_ctrl2.png", width = 3, height = 3)

cml1_seurat_nk_expanded_tumor_pvr1.2_ctrl1.2_markers <- FindMarkers(cml1_seurat_nk_expanded_tumor_pvr1.2_ctrl1.2, logfc.threshold = 0.15, 
                                                                    ident.1 = "PVR KO", ident.2 = "PVR WT", 
                                                                    test.use = "t", max.cells.per.ident = 1e3) %>%
  filter(p_val_adj < 0.05) %>% add_rownames(var = "gene")

ggplot(data = cml1_seurat_nk_expanded_tumor_pvr1.2_ctrl1.2_markers, aes(avg_log2FC, -log10(p_val_adj), label = gene)) + geom_point() + 
  ggrepel::geom_text_repel(data = subset(cml1_seurat_nk_expanded_tumor_pvr1.2_ctrl1.2_markers, p_val_adj < 1e-10)) +
  ggpubr::theme_pubr()
ggsave("results/functional/hto1/volcano_exp_nk_pvr1.2_ko_ctrl1.2.png", width = 9, height = 9)







#### =========


## Only non-expanded
cml1_seurat_nk_nonexpanded <- subset(cml1_seurat_nk, expanded != "Expanded NK")
cml1_seurat_nk_nonexpanded <- cml1_seurat_nk_nonexpanded %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml1_seurat_nk_nonexpanded), nPCs = 10)
cml1_seurat_nk_nonexpanded@reductions$pca@stdev %>% plot()
cml1_seurat_nk_nonexpanded <- cml1_seurat_nk_nonexpanded %>% getClustering()
Idents(cml1_seurat_nk_nonexpanded) <- cml1_seurat_nk_nonexpanded$RNA_snn_res.0.1

a <- DimPlot(cml1_seurat_nk_nonexpanded, group.by = "my.demux", cols = getPalette(10))
b <- DimPlot(cml1_seurat_nk_nonexpanded, group.by = "covariate", cols = getPalette(5))
c <- DimPlot(cml1_seurat_nk_nonexpanded, group.by = "expanded", cols = getPalette3(4))
d <- DimPlot(cml1_seurat_nk_nonexpanded, group.by = "pvr_ko", cols = getPalette3(4))
e <- DimPlot(cml1_seurat_nk_nonexpanded, group.by = "Phase", cols = getPalette3(4))
f <- DimPlot(cml1_seurat_nk_nonexpanded, group.by = "RNA_snn_res.0.1", cols = getPalette3(4))

a + b + c + d + e + f
ggsave("results/functional/hto1/umap_meta_nonexp_nk.png", width = 12, height = 5)







## Tumor vs no tumor, nonexpanded
cml1_seurat_nk_nonexpanded_tumor_no_tumor <- subset(cml1_seurat_nk, my.demux %in% c("K562-ctrl1-T-NE-NK-", "K562-ctrl2-T-NE-NK-", "NE-NK-only"))
cml1_seurat_nk_nonexpanded_tumor_no_tumor <- cml1_seurat_nk_nonexpanded_tumor_no_tumor %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml1_seurat_nk_nonexpanded_tumor_no_tumor), nPCs = 10)
cml1_seurat_nk_nonexpanded_tumor_no_tumor@reductions$pca@stdev %>% plot()
Idents(cml1_seurat_nk_nonexpanded_tumor_no_tumor) <- cml1_seurat_nk_nonexpanded_tumor_no_tumor$tumor

DimPlot(cml1_seurat_nk_nonexpanded_tumor_no_tumor, group.by = "my.demux", cols = c("blue3", "dodgerblue","red3"))
ggsave("results/functional/hto1/umap_meta_nonexp_nk_tumor_only.png", width = 5, height = 3)

FeaturePlot(cml1_seurat_nk_nonexpanded_tumor_no_tumor, features = c("TNFRSF18", "CRTAM", "GZMA", "NKG7"), cols = c("gray90", "red3"))
ggsave("results/functional/hto1/umap_meta_nonexp_nk_tumor_only_features.png", width = 6, height = 4)

cml1_seurat_nk_nonexpanded_tumor_no_tumor_markers <- FindMarkers(cml1_seurat_nk_nonexpanded_tumor_no_tumor, logfc.threshold = 0.15, ident.1 = "Tumor", ident.2 = "No tumor", 
                                                              test.use = "t", max.cells.per.ident = 1e3) %>%
  filter(p_val_adj < 0.05) %>% add_rownames(var = "gene")
fwrite(cml1_seurat_nk_nonexpanded_tumor_no_tumor_markers, "results/functional/hto1/volcano_nonexp_nk_tumor_vs_not_tumor.txt", sep = "\t", quote = F, row.names = F)


ggplot(data = cml1_seurat_nk_nonexpanded_tumor_no_tumor_markers, aes(avg_log2FC, -log10(p_val_adj), label = gene)) + geom_point() + 
  ggrepel::geom_text_repel(data = subset(cml1_seurat_nk_nonexpanded_tumor_no_tumor_markers, p_val_adj < 1e-10)) +
  ggpubr::theme_pubr() + xlim(values=c(-3.5,3.5))
ggsave("results/functional/hto1/volcano_nonexp_nk_tumor_vs_not_tumor.png", width = 9, height = 9)


## Only nonexpanded, with tumor
cml1_seurat_nk_nonexpanded_tumor <- subset(cml1_seurat_nk, tumor == "Tumor")
cells.to.keep <- cml1_seurat_nk_nonexpanded_tumor@meta.data %>% filter(grepl("NE", cml1_seurat_nk_nonexpanded_tumor$my.demux)) %>% pull(barcode)
cml1_seurat_nk_nonexpanded_tumor <- subset(cml1_seurat_nk_nonexpanded_tumor, cells = cells.to.keep)
cml1_seurat_nk_nonexpanded_tumor <- cml1_seurat_nk_nonexpanded_tumor %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml1_seurat_nk_nonexpanded_tumor), nPCs = 20)
cml1_seurat_nk_nonexpanded_tumor@reductions$pca@stdev %>% plot()

Idents(cml1_seurat_nk_nonexpanded_tumor) <- cml1_seurat_nk_nonexpanded_tumor$my.demux
DimPlot(cml1_seurat_nk_nonexpanded_tumor)

Idents(cml1_seurat_nk_nonexpanded_tumor) <- cml1_seurat_nk_nonexpanded_tumor$pvr_ko
DimPlot(cml1_seurat_nk_nonexpanded_tumor)



## Only nonexpanded
cml1_seurat_nk_nonexpanded <- subset(cml1_seurat_nk, expanded != "Expanded NK")
cml1_seurat_nk_nonexpanded <- cml1_seurat_nk_nonexpanded %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml1_seurat_nk_nonexpanded), nPCs = 12)
cml1_seurat_nk_nonexpanded@reductions$pca@stdev %>% plot()

Idents(cml1_seurat_nk_nonexpanded) <- cml1_seurat_nk_nonexpanded$my.demux
DimPlot(cml1_seurat_nk_nonexpanded, group.by = "my.demux", cols = c("red3", "blue3", "dodgerblue", "green3", "green4")) + NoLegend()
ggsave("results/functional/hto1/umap_meta_nonexp_nk.png", width = 3, height = 3)



## Only nonexpanded, with tumor
cml1_seurat_nk_nonexpanded_tumor <- subset(cml1_seurat_nk_nonexpanded, tumor == "Tumor")
cml1_seurat_nk_nonexpanded_tumor <- cml1_seurat_nk_nonexpanded_tumor %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml1_seurat_nk_nonexpanded_tumor), nPCs = 16)
cml1_seurat_nk_nonexpanded_tumor@reductions$pca@stdev %>% plot()

# cml1_seurat_nk_nonexpanded_tumor <- RunUMAP(cml1_seurat_nk_nonexpanded_tumor, dims = 2:12)
Idents(cml1_seurat_nk_nonexpanded_tumor) <- cml1_seurat_nk_nonexpanded_tumor$my.demux
DimPlot(cml1_seurat_nk_nonexpanded_tumor, group.by = "my.demux", cols = c("blue3", "dodgerblue", "green3", "green4")) + NoLegend()
ggsave("results/functional/hto1/umap_meta_nonexp_nk_tumor.png", width = 3, height = 3)

DimPlot(cml1_seurat_nk_nonexpanded_tumor, group.by = "Phase") #, cols = c("blue3", "dodgerblue", "green3", "green4")) + NoLegend()

## Only nonexpanded, with tumor
cml1_seurat_nk_nonexpanded_tumor <- subset(cml1_seurat_nk, expanded != "Expanded NK")
cml1_seurat_nk_nonexpanded_tumor <- subset(cml1_seurat_nk_nonexpanded_tumor, tumor != "Tumor")

cml1_seurat_nk_nonexpanded_tumor <- cml1_seurat_nk_nonexpanded_tumor %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml1_seurat_nk_nonexpanded_tumor), nPCs = 6)
cml1_seurat_nk_nonexpanded_tumor@reductions$pca@stdev %>% plot()

Idents(cml1_seurat_nk_nonexpanded) <- cml1_seurat_nk_nonexpanded$my.demux
DimPlot(cml1_seurat_nk_nonexpanded, group.by = "my.demux", cols = c("red3", "blue3", "dodgerblue", "green3", "green4")) + NoLegend()
ggsave("results/functional/hto1/umap_meta_nonexp_nk.png", width = 3, height = 3)



## Only nonexpanded, with tumor; PVR1.1 guide and ctrl1 guide
cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.1 <- subset(cml1_seurat_nk_nonexpanded_tumor, my.demux %in% c("K562-ctrl1-T-NE-NK-", "K562-PVR1-T-NE-NK-"))
cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.1 <- cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.1 %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.1), nPCs = 22)
cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.1@reductions$pca@stdev %>% plot()
Idents(cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.1) <- cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.1$pvr_ko

Idents(cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.1) <- cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.1$my.demux
DimPlot(cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.1, group.by = "my.demux", cols = c("blue3", "green4")) + NoLegend()
ggsave("results/functional/hto1/umap_pvr1_1_ctrl1_1.png", width = 3, height = 3)

Idents(cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.1) <- cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.1$pvr_ko
DimPlot(cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.1)
ggsave("results/functional/hto1/umap_pvr1_1_ctrl1_1_2.png", width = 4, height = 3)

cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.1_markers <- FindMarkers(cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.1, logfc.threshold = 0.15, ident.1 = "PVR KO", ident.2 = "PVR WT", 
                                                                    test.use = "t", max.cells.per.ident = 1e3) %>%
  filter(p_val_adj < 0.05) %>% add_rownames(var = "gene")

ggplot(data = cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.1_markers, aes(avg_log2FC, -log10(p_val_adj), label = gene)) + geom_point() + 
  ggrepel::geom_text_repel(data = subset(cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.1_markers, p_val_adj < 1e-10)) +
  ggpubr::theme_pubr()
ggsave("results/functional/hto1/volcano_nonexp_nk_pvr_ko_pvr_wt.png", width = 9, height = 9)








## Only nonexpanded, with tumor; PVR1.1 guide and ctrl2 guide
cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.2 <- subset(cml1_seurat_nk_nonexpanded_tumor, my.demux %in% c("K562-ctrl2-T-NE-NK-", "K562-PVR1-T-NE-NK-"))
cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.2 <- cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.2 %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.2), nPCs = 10)
cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.2@reductions$pca@stdev %>% plot()
Idents(cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.2) <- cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.2$pvr_ko

Idents(cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.2) <- cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.2$my.demux
DimPlot(cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.2, group.by = "my.demux", cols = c("blue3", "green3", "green4")) + NoLegend()
ggsave("results/functional/hto1/umap_pvr1_1_ctrl1_2.png", width = 3, height = 3)

Idents(cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.2) <- cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.2$pvr_ko
DimPlot(cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.2)
ggsave("results/functional/hto1/umap_pvr1_1_ctrl1_2_2.png", width = 4, height = 3)

cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.2_markers <- FindMarkers(cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.2, logfc.threshold = 0.15, ident.1 = "PVR KO", ident.2 = "PVR WT", 
                                                                    test.use = "t", max.cells.per.ident = 1e3) %>%
  filter(p_val_adj < 0.05) %>% add_rownames(var = "gene")

ggplot(data = cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.2_markers, aes(avg_log2FC, -log10(p_val_adj), label = gene)) + geom_point() + 
  ggrepel::geom_text_repel(data = subset(cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.2_markers, p_val_adj < 1e-1)) +
  ggpubr::theme_pubr()
ggsave("results/functional/hto1/volcano_enonxp_nk_pvr1.1_ko_ctrl1.2.png", width = 9, height = 9)







## Only nonexpanded, with tumor; PVR2 guide and ctrl1 guide
cml1_seurat_nk_nonexpanded_tumor_pvr1.2_ctrl1.1 <- subset(cml1_seurat_nk_nonexpanded_tumor, my.demux %in% c("K562-ctrl1-T-NE-NK-", "K562-PVR2-T-NE-NK-"))
cml1_seurat_nk_nonexpanded_tumor_pvr1.2_ctrl1.1 <- cml1_seurat_nk_nonexpanded_tumor_pvr1.2_ctrl1.1 %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml1_seurat_nk_nonexpanded_tumor_pvr1.2_ctrl1.1), nPCs = 10)
cml1_seurat_nk_nonexpanded_tumor_pvr1.2_ctrl1.1@reductions$pca@stdev %>% plot()
Idents(cml1_seurat_nk_nonexpanded_tumor_pvr1.2_ctrl1.1) <- cml1_seurat_nk_nonexpanded_tumor_pvr1.2_ctrl1.1$pvr_ko

Idents(cml1_seurat_nk_nonexpanded_tumor_pvr1.2_ctrl1.1) <- cml1_seurat_nk_nonexpanded_tumor_pvr1.2_ctrl1.1$my.demux
DimPlot(cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.1, group.by = "my.demux", cols = c("blue3", "green4")) + NoLegend()
ggsave("results/functional/hto1/umap_pvr1_2_ctrl1_1.png", width = 3, height = 3)

Idents(cml1_seurat_nk_nonexpanded_tumor_pvr1.2_ctrl1.1) <- cml1_seurat_nk_nonexpanded_tumor_pvr1.2_ctrl1.1$pvr_ko
DimPlot(cml1_seurat_nk_nonexpanded_tumor_pvr1.2_ctrl1.1)
ggsave("results/functional/hto1/umap_pvr1_2_ctrl1_2_1.png", width = 4, height = 3)

cml1_seurat_nk_nonexpanded_tumor_pvr1.2_ctrl1.1_markers <- FindMarkers(cml1_seurat_nk_nonexpanded_tumor_pvr1.2_ctrl1.1, logfc.threshold = 0.15, 
                                                                    ident.1 = "PVR KO", ident.2 = "PVR WT", 
                                                                    test.use = "t", max.cells.per.ident = 1e3) %>%
  filter(p_val_adj < 0.05) %>% add_rownames(var = "gene")

ggplot(data = cml1_seurat_nk_nonexpanded_tumor_pvr1.2_ctrl1.1_markers, aes(avg_log2FC, -log10(p_val_adj), label = gene)) + geom_point() + 
  ggrepel::geom_text_repel(data = subset(cml1_seurat_nk_nonexpanded_tumor_pvr1.2_ctrl1.1_markers, p_val_adj < 1e-10)) +
  ggpubr::theme_pubr()
ggsave("results/functional/hto1/volcano_exp_nk_pvr1.2_ko_ctrl1.1.png", width = 9, height = 9)






## Only nonexpanded, with tumor; PVR1.2 guide and ctrl2 guide
cml1_seurat_nk_nonexpanded_tumor_pvr1.2_ctrl1.2 <- subset(cml1_seurat_nk_nonexpanded_tumor, my.demux %in% c("K562-ctrl2-T-NE-NK-", "K562-PVR2-T-NE-NK-"))
cml1_seurat_nk_nonexpanded_tumor_pvr1.2_ctrl1.2 <- cml1_seurat_nk_nonexpanded_tumor_pvr1.2_ctrl1.2 %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml1_seurat_nk_nonexpanded_tumor_pvr1.2_ctrl1.2), nPCs = 10)
cml1_seurat_nk_nonexpanded_tumor_pvr1.2_ctrl1.2@reductions$pca@stdev %>% plot()
Idents(cml1_seurat_nk_nonexpanded_tumor_pvr1.2_ctrl1.2) <- cml1_seurat_nk_nonexpanded_tumor_pvr1.2_ctrl1.2$pvr_ko

Idents(cml1_seurat_nk_nonexpanded_tumor_pvr1.2_ctrl1.2) <- cml1_seurat_nk_nonexpanded_tumor_pvr1.2_ctrl1.2$my.demux
DimPlot(cml1_seurat_nk_nonexpanded_tumor_pvr1.1_ctrl1.1, group.by = "my.demux", cols = c("dodgerblue", "green4")) + NoLegend()
ggsave("results/functional/hto1/umap_pvr1_2_ctrl1_2.png", width = 3, height = 3)

Idents(cml1_seurat_nk_nonexpanded_tumor_pvr1.2_ctrl1.2) <- cml1_seurat_nk_nonexpanded_tumor_pvr1.2_ctrl1.2$pvr_ko
DimPlot(cml1_seurat_nk_nonexpanded_tumor_pvr1.2_ctrl1.2)
ggsave("results/functional/hto1/umap_pvr1_1_ctrl1_2_2.png", width = 4, height = 3)

cml1_seurat_nk_nonexpanded_tumor_pvr1.2_ctrl1.2_markers <- FindMarkers(cml1_seurat_nk_nonexpanded_tumor_pvr1.2_ctrl1.2, logfc.threshold = 0.15, 
                                                                    ident.1 = "PVR KO", ident.2 = "PVR WT", 
                                                                    test.use = "t", max.cells.per.ident = 1e3) %>%
  filter(p_val_adj < 0.05) %>% add_rownames(var = "gene")

ggplot(data = cml1_seurat_nk_nonexpanded_tumor_pvr1.2_ctrl1.2_markers, aes(avg_log2FC, -log10(p_val_adj), label = gene)) + geom_point() + 
  ggrepel::geom_text_repel(data = subset(cml1_seurat_nk_nonexpanded_tumor_pvr1.2_ctrl1.2_markers, p_val_adj < 1e-10)) +
  ggpubr::theme_pubr()
ggsave("results/functional/hto1/volcano_nonexp_nk_pvr1.2_ko_ctrl1.2.png", width = 9, height = 9)

