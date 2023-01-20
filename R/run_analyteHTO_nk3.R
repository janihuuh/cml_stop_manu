
cml3_seurat <- cml_seurat_new_filt3

### Show that KO is working?
cml3_seurat_tumor <- subset(cml3_seurat, is_nk != "NK")
cml3_seurat_tumor <- subset(cml3_seurat_tumor, my.demux != "E-NK-only")
cml3_seurat_tumor <- cml3_seurat_tumor %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml3_seurat_tumor), nPCs = 14)
cml3_seurat_tumor@reductions$pca@stdev %>% plot()

cml3_seurat_tumor <- AddModuleScore(cml3_seurat_tumor, features = list(c("PVR")), name = "PVR")


cml3_seurat_tumor@meta.data %>% 
  ggplot(aes(my.demux,PVR1,fill=my.demux)) + geom_violin(draw_quantiles = 0.5, adjust = 3) + geom_jitter(size=0.001,alpha=0.3) +
  ggpubr::stat_compare_means(label="p", label.y.npc = 0.9) + ggpubr::rotate_x_text(angle = 45) + 
  #scale_fill_manual(values=c("blue3","green3","green4")) +
  NoLegend() + labs(x="",y="PVR") 
ggsave("results/functional/hto3/vln_pvr_ko_total.pdf", width = 4, height = 3)

cml3_seurat_tumor@meta.data %>% group_by(my.demux, pvr_pos = PVR1>0) %>% tally() %>% mutate(prop=n/sum(n)) %>% 
  filter(pvr_pos) %>% 
  
  # ggplot(aes(reorder(my.demux,-prop),prop*100,fill=my.demux)) + geom_bar(stat="identity") +
  ggplot(aes(my.demux,prop*100,fill=my.demux)) + geom_bar(stat="identity") +

  # ggpubr::stat_compare_means(label="p") + 
  ggpubr::rotate_x_text(angle = 45) +
  #scale_fill_manual(values=c("blue3","green3","green4")) +
  labs(x="",y="% of PVR pos cells") + NoLegend() 
ggsave("results/functional/hto3/bar_pvr_ko_total.pdf", width = 3, height = 4)





cml3_seurat_tumor_pvr <- subset(cml3_seurat, my.demux %in% c("LAMA84-PVR1-T-","LAMA84-PVR2-T-", "LAMA84-ctrl2-T-"))
cml3_seurat_tumor_pvr <- subset(cml3_seurat_tumor_pvr, is_nk != "NK")
cml3_seurat_tumor_pvr <- cml3_seurat_tumor_pvr %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml3_seurat_tumor_pvr), nPCs = 11)
cml3_seurat_tumor_pvr@reductions$pca@stdev %>% plot()

cml3_seurat_tumor_pvr <- AddModuleScore(cml3_seurat_tumor_pvr, features = list(c("PVR")), name = "PVR")

cml3_seurat_tumor_pvr@meta.data %>% 
  ggplot(aes(my.demux,PVR1,fill=my.demux)) + geom_violin(draw_quantiles = 0.5, adjust = 3) + geom_jitter(size=0.001,alpha=0.3) +
  ggpubr::stat_compare_means(label="p", label.y.npc = 0.9) + ggpubr::rotate_x_text(angle = 45) + scale_fill_manual(values=c("blue3","green3","green4")) +
  NoLegend() + labs(x="",y="PVR") 
ggsave("results/functional/hto3/vln_pvr_ko.pdf", width = 2, height = 3)

cml3_seurat_tumor_pvr@meta.data %>% group_by(my.demux, pvr_pos = PVR1>0) %>% tally() %>% mutate(prop=n/sum(n)) %>% 
  filter(pvr_pos) %>% 
  
  ggplot(aes(reorder(my.demux,-prop),prop*100,fill=my.demux)) + geom_bar(stat="identity") +
  ggpubr::stat_compare_means(label="p") + ggpubr::rotate_x_text(angle = 45) +
  labs(x="",y="% of PVR pos cells") + NoLegend() + scale_fill_manual(values=c("blue3","green3","green4"))
ggsave("results/functional/hto3/bar_pvr_ko.pdf", width = 2, height = 3)


## Focus on NK cells; select NK cells based on i) Essi's demultiplex and ii) RNA-profile of NK cells
cml3_seurat$is_nk <- ifelse(cml3_seurat$RNA_snn_res.0.1 == 1 & grepl("NK", cml3_seurat$my.demux), "NK", "CML")

cml3_seurat_nk <- subset(cml3_seurat, is_nk == "NK")
cml3_seurat_nk <- cml3_seurat_nk %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml3_seurat_nk), nPCs = 10)
cml3_seurat_nk@reductions$pca@stdev %>% plot()
cml3_seurat_nk <- cml3_seurat_nk %>% getClustering()

cml3_seurat_nk$nk_presence <- ifelse(grepl(cml3_seurat_nk$my.demux, pattern = "NK"), "NK presence", "No NK")
cml3_seurat_nk$pvr_ko      <- ifelse(grepl(cml3_seurat_nk$my.demux, pattern = "PVR"), "PVR KO", "PVR WT")
cml3_seurat_nk$tumor       <- ifelse(grepl(cml3_seurat_nk$my.demux, pattern = "NK-only"), "No tumor", "Tumor")

cml3_seurat_nk$expanded    <- ifelse(grepl(cml3_seurat_nk$my.demux, pattern = "NE"), "Non-expanded NK", "Expanded NK")
cml3_seurat_nk$covariate   <- paste0(cml3_seurat_nk$expanded , "\n", cml3_seurat_nk$pvr_ko)

a <- DimPlot(cml3_seurat_nk, group.by = "my.demux", cols = getPalette(10))
b <- DimPlot(cml3_seurat_nk, group.by = "covariate", cols = getPalette(5))
c <- DimPlot(cml3_seurat_nk, group.by = "expanded", cols = getPalette3(4))
d <- DimPlot(cml3_seurat_nk, group.by = "pvr_ko", cols = getPalette3(4))
e <- DimPlot(cml3_seurat_nk, group.by = "tumor", cols = getPalette3(3))

a + b + c + d + e
ggsave("results/functional/hto3/umap_meta_nk.png", width = 15, height = 7)




## Only expanded based on Essi's demultiplexing
cml3_seurat_nk$expanded <- ifelse(grepl(cml3_seurat_nk$my.demux, pattern = "NE"), "Non-expanded NK", "Expanded NK")
cml3_seurat_nk_expanded <- subset(cml3_seurat_nk, expanded == "Expanded NK")
cml3_seurat_nk_expanded <- cml3_seurat_nk_expanded %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml3_seurat_nk_expanded), nPCs = 11)
cml3_seurat_nk_expanded@reductions$pca@stdev %>% plot()

a <- DimPlot(cml3_seurat_nk_expanded, group.by = "my.demux", cols = getPalette(10))
b <- DimPlot(cml3_seurat_nk_expanded, group.by = "covariate", cols = getPalette(5))
c <- DimPlot(cml3_seurat_nk_expanded, group.by = "expanded", cols = getPalette3(4))
d <- DimPlot(cml3_seurat_nk_expanded, group.by = "pvr_ko", cols = getPalette3(4))
e <- DimPlot(cml3_seurat_nk_expanded, group.by = "tumor", cols = getPalette3(3))

a + b + c + d + e
ggsave("results/functional/hto3/umap_meta_exp_nk.png", width = 12, height = 4)

### Satu's edit
cml3_seurat_nk_expanded$satu <- ifelse(cml3_seurat_nk_expanded$my.demux %in% c("E-NK-only"), "No tumor", "other")
cml3_seurat_nk_expanded$satu <- ifelse(cml3_seurat_nk_expanded$my.demux %in% c("LAMA84-ctrl1-T-E-NK-", "LAMA84-ctrl2-T-E-NK-"), "Tumor PVR wt", cml3_seurat_nk_expanded$satu)
cml3_seurat_nk_expanded$satu <- ifelse(cml3_seurat_nk_expanded$my.demux %in% c("LAMA84-PVR1-T-E-NK-", "LAMA84-PVR2-T-E-NK-"), "Tumor PVR ko", cml3_seurat_nk_expanded$satu)
table(cml3_seurat_nk_expanded$my.demux)
table(cml3_seurat_nk_expanded$satu)
cml3_seurat_nk_expanded$satu <- factor(cml3_seurat_nk_expanded$satu, levels = c("Tumor PVR ko", "Tumor PVR wt", "No tumor"))

cml3_seurat_nk_expanded <- AddModuleScore(cml3_seurat_nk_expanded, features = list(activation$activated), name = "activation")

dufva_signatures_full <- fread("data/dufva_nk_signatures_full.txt", dec = ",") %>% filter(p_val_adj < 0.05) %>% dplyr::select(cluster, gene)
dufva_signatures_activated_full <- dufva_signatures_full %>% filter(cluster == "Activated (2)") %>% pull(gene)

markers3_2 %>% filter(!gene %in% markers3_1$gene) %>% filter(!gene %in% dufva_signatures_full$gene) 
markers3_2 %>% filter(!gene %in% markers3_1$gene) %>% arrange(-(avg_log2FC)) %>% filter(avg_log2FC > 0) %>% filter(!gene %in% unwanted_genes) %>% filter(gene %in% dufva_signatures_activated_full) %>% pull(gene)


cml3_seurat_nk_expanded <- AddModuleScore(cml3_seurat_nk_expanded, features = list(dufva_signatures_activated_full), name = "activation_full")
cml3_seurat_nk_expanded <- AddModuleScore(cml3_seurat_nk_expanded, features = list(head(dufva_signatures_activated_full,50)), name = "activation_top50")


dufva_signatures_activated_full

Idents(cml3_seurat_nk_expanded) <- cml3_seurat_nk_expanded$satu
DimPlot(cml3_seurat_nk_expanded)

markers3_1 <- FindMarkers(cml3_seurat_nk_expanded, ident.1 = "Tumor PVR wt", ident.2 = "No tumor", test.use = "t") %>% add_rownames(var = "gene")# %>% filter(p_val_adj < 0.05)
markers3_2 <- FindMarkers(cml3_seurat_nk_expanded, ident.1 = "Tumor PVR ko", ident.2 = "No tumor", test.use = "t") %>% add_rownames(var = "gene")# %>% filter(p_val_adj < 0.05)
markers3_3 <- FindMarkers(cml3_seurat_nk_expanded, ident.1 = "Tumor PVR ko", ident.2 = "Tumor PVR wt", test.use = "t") %>% add_rownames(var = "gene")# %>% filter(p_val_adj < 0.05)

markers3_1_sigf <- subset(markers3_1, p_val_adj < 0.05)
markers3_2_sigf <- subset(markers3_2, p_val_adj < 0.05)

markers3_1_sigf_exclusive_genes <- markers3_1_sigf$gene[!markers3_1_sigf$gene %in% markers3_2_sigf$gene]
markers3_2_sigf_exclusive_genes <- markers3_2_sigf$gene[!markers3_2_sigf$gene %in% markers3_1_sigf$gene]

full_df3 <- markers3_2 %>% full_join(markers3_1, by = "gene") 

full_df3$avg_log2FC.x <- ifelse(is.na(full_df3$avg_log2FC.x), 0, full_df3$avg_log2FC.x)
full_df3$avg_log2FC.y <- ifelse(is.na(full_df3$avg_log2FC.y), 0, full_df3$avg_log2FC.y)

full_df3$gene_sigf <- ifelse(full_df3$gene %in% markers3_1_sigf_exclusive_genes, "Sigf only in PVR WT", "Other")
full_df3$gene_sigf <- ifelse(full_df3$gene %in% markers3_2_sigf_exclusive_genes, "Sigf only in PVR KO", full_df3$gene_sigf)
full_df3$activation <- ifelse(full_df3$gene %in% dufva_signatures_activated_full, "activation", "other")

ggplot() +
  geom_point(data = full_df3, aes(avg_log2FC.x, avg_log2FC.y, color = gene_sigf, label = gene, shape = activation)) + 
  ggrepel::geom_text_repel(data = subset(full_df3, gene_sigf == "Sigf only in PVR KO" & activation == "activation"), aes(avg_log2FC.x, avg_log2FC.y,  label = gene, shape = activation), fontface = "italic") + 
  
  geom_abline() +
  xlim(values = c(0,4)) + ylim(values = c(0,4)) + scale_shape_manual(values = c(17,16)) + scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = c("gray90", "red", "blue")) + add_guide + labs(x = "avg log2FC\nin NKs against PVR KO\nand NK alone", y = "avg log2FC\nin NKs against PVR WT\nand NK alone", color = "", shape = "associated\nwith NK activation")
ggsave("results/functional/scatter_cml3_pvr_ko_pvr_wt.png", width = 7, height = 5)
  
fwrite(full_df3, "results/functional/hto3/volcano_nonexp_nk_pvr_ko_pvr1.2_wt.txt", sep = "\t", quote = F, row.names = F)


  
  # + ggrepel::geom_text_repel(color="black")



markers3_2 %>% filter(!gene %in% markers3_1$gene) %>% arrange(-(avg_log2FC)) %>% filter(avg_log2FC > 0) %>% filter(!gene %in% unwanted_genes) %>% pull(gene)
VlnPlot(cml3_seurat_nk_expanded, features = c("ODC1", "BATF3", "MYC", "SOCS1", "TIGIT", "HAVCR2", "activation1"))

cml3_seurat_nk_expanded$activation2 <- ifelse(cml3_seurat_nk_expanded$activation1 < 0, 0, cml3_seurat_nk_expanded$activation1)
cml3_seurat_nk_expanded$activation2 <- ifelse(cml3_seurat_nk_expanded$activation1 < 0, 0, cml3_seurat_nk_expanded$activation1)

cml3_seurat_nk_expanded@meta.data %>% 
  ggplot(aes(satu, activation1)) + geom_violin(draw_quantiles = c(0.5), adjust = 1) + geom_jitter(size = 0.1) + ggpubr::stat_compare_means(label = "p") + 
  ggsignif::geom_signif(comparisons = list(c("Tumor PVR ko", "Tumor PVR wt")))

cml3_seurat_nk_expanded@meta.data %>% 
  ggplot(aes(satu, activation2)) + geom_violin(draw_quantiles = c(0.5), adjust = 1) + geom_jitter(size = 0.1) + ggpubr::stat_compare_means(label = "p") + 
  ggsignif::geom_signif(comparisons = list(c("Tumor PVR ko", "Tumor PVR wt")))

cml3_seurat_nk_expanded@meta.data %>% 
  ggplot(aes(satu, activation_full1)) + geom_violin(draw_quantiles = c(0.5), adjust = 1) + geom_jitter(size = 0.1) + ggpubr::stat_compare_means(label = "p") + 
  ggsignif::geom_signif(comparisons = list(c("Tumor PVR ko", "Tumor PVR wt")))

cml3_seurat_nk_expanded@meta.data %>% 
  ggplot(aes(satu, activation_top501)) + geom_violin(draw_quantiles = c(0.5), adjust = 1) + geom_jitter(size = 0.1) + ggpubr::stat_compare_means(label = "p") + 
  ggsignif::geom_signif(comparisons = list(c("Tumor PVR ko", "Tumor PVR wt")))

DimPlot(cml3_seurat_nk_expanded, group.by = "RNA_snn_res.0.1", cols = getPalette3(4))



a <- markers3_2 %>% filter(!gene %in% markers3_1$gene) %>% arrange(-(avg_log2FC)) %>% filter(avg_log2FC > 0) %>% filter(!gene %in% unwanted_genes) #%>% pull(gene) 
getHypergeometric(genes_df = a, universe_df = rownames(cml_seurat), term_df = hallmark) %>% filter(qvalue < 0.05)
getHypergeometric(genes_df = a, universe_df = rownames(cml_seurat), term_df = kegg) %>% filter(qvalue < 0.05)
getHypergeometric(genes_df = a, universe_df = rownames(cml_seurat), term_df = reactome) %>% filter(qvalue < 0.05)
getHypergeometric(genes_df = a, universe_df = rownames(cml_seurat), term_df = go_pathways) %>% filter(qvalue < 0.05)

a <- markers3_2 %>% filter(!gene %in% markers3_1$gene) %>% arrange(-(avg_log2FC)) %>% filter(avg_log2FC < 0) %>% filter(!gene %in% unwanted_genes) #%>% pull(gene) 
getHypergeometric(genes_df = a, universe_df = rownames(cml_seurat), term_df = hallmark) %>% filter(qvalue < 0.05)


# getHypergeometric(genes_df = a, universe_df = rownames(cml_seurat), term_df = imm) %>% filter(qvalue < 0.05)


### Take only the activated NK cells
cml3_seurat_nk_expanded_activated <- subset(cml3_seurat_nk_expanded, RNA_snn_res.0.1 == 1)
cml3_seurat_nk_expanded_activated <- subset(cml3_seurat_nk_expanded_activated, my.demux != "E-NK-only")
cml3_seurat_nk_expanded_activated <- cml3_seurat_nk_expanded_activated %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml3_seurat_nk_expanded_activated), nPCs = 6)
cml3_seurat_nk_expanded_activated@reductions$pca@stdev %>% plot()

DimPlot(cml3_seurat_nk_expanded_activated, group.by = "my.demux", cols = getPalette(10))
DimPlot(cml3_seurat_nk_expanded_activated, group.by = "RNA_snn_res.0.5", cols = getPalette(10))

cml3_seurat_nk_expanded_activated@meta.data %>% 
  group_by(my.demux,cl=RNA_snn_res.0.5) %>% tally() %>% mutate(prop=n/sum(n)) %>% 
  
  ggplot(aes(my.demux,prop,fill=cl)) + geom_bar(stat="identity") + ggpubr::rotate_x_text(angle = 45)

Idents(cml3_seurat_nk_expanded_activated) <- cml3_seurat_nk_expanded_activated$RNA_snn_res.0.5
cml3_seurat_nk_expanded_activated_markers <- FindAllMarkers(cml3_seurat_nk_expanded_activated, logfc.threshold = 0.15, test.use = "t", only.pos = T) %>% filter(p_val_adj < 0.05)
cml3_seurat_nk_expanded_activated_markers

## Tumor vs no tumor, expanded
cml3_seurat_nk_expanded_tumor_no_tumor <- subset(cml3_seurat_nk, my.demux %in% c("LAMA84-ctrl2-T-E-NK-", "E-NK-only"))
cml3_seurat_nk_expanded_tumor_no_tumor <- cml3_seurat_nk_expanded_tumor_no_tumor %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml3_seurat_nk_expanded_tumor_no_tumor), nPCs = 11)
cml3_seurat_nk_expanded_tumor_no_tumor@reductions$pca@stdev %>% plot()
Idents(cml3_seurat_nk_expanded_tumor_no_tumor) <- cml3_seurat_nk_expanded_tumor_no_tumor$tumor

a <- DimPlot(cml3_seurat_nk_expanded_tumor_no_tumor, group.by = "my.demux", cols = getPalette(10))
b <- DimPlot(cml3_seurat_nk_expanded_tumor_no_tumor, group.by = "covariate", cols = getPalette(5))
c <- DimPlot(cml3_seurat_nk_expanded_tumor_no_tumor, group.by = "expanded", cols = getPalette3(4))
d <- DimPlot(cml3_seurat_nk_expanded_tumor_no_tumor, group.by = "pvr_ko", cols = getPalette3(4))
e <- DimPlot(cml3_seurat_nk_expanded_tumor_no_tumor, group.by = "tumor", cols = getPalette3(3))

a + b+ c+ d+ e 
ggsave("results/functional/hto3/umap_meta_exp_nk_tumor_no_tumor.png", width = 5, height = 3)


cml3_seurat_nk_expanded_tumor_no_tumor@reductions$umap@cell.embeddings

cml3_seurat_nk_expanded_tumor_no_tumor@reductions$umap@cell.embeddings[,1] <- -cml3_seurat_nk_expanded_tumor_no_tumor@reductions$umap@cell.embeddings[,1]


DimPlot(cml3_seurat_nk_expanded_tumor_no_tumor, group.by = "tumor", cols = rev(c("darkgoldenrod4", "darkslategray4"))) + theme_void(base_size = 17) + NoLegend() + labs(title = "")
ggsave("results/functional/umap_3_exp_nk_tumor_only.png", width = 2, height = 2)

DimPlot(cml3_seurat_nk_expanded_tumor_no_tumor, group.by = "tumor", cols = c("darkgoldenrod4", "darkslategray4")) + theme_void(base_size = 17) + labs(title = "") + add_guide
ggsave("results/functional/umap_3_exp_nk_tumor_only2.png", width = 2, height = 2)


cml3_seurat_nk_expanded_tumor_no_tumor_markers <- FindMarkers(cml3_seurat_nk_expanded_tumor_no_tumor, logfc.threshold = 0.15, ident.1 = "Tumor", ident.2 = "No tumor", 
                                                              test.use = "t", max.cells.per.ident = 1e3) %>%
  filter(p_val_adj < 0.05) %>% add_rownames(var = "gene")
fwrite(cml3_seurat_nk_expanded_tumor_no_tumor_markers, "results/functional/cml3_seurat_nk_expanded_tumor_no_tumor_markers.txt", sep = "\t", quote = F, row.names = F)

ggplot() +
  geom_point(data = cml3_seurat_nk_expanded_tumor_no_tumor_markers, aes(avg_log2FC, -log10(p_val), color = avg_log2FC > 0)) + 
  ggrepel::geom_text_repel(data = subset(cml3_seurat_nk_expanded_tumor_no_tumor_markers, -log10(cml3_seurat_nk_expanded_tumor_no_tumor_markers$p_val) > 15), aes(avg_log2FC, -log10(p_val), label = gene), fontface = "italic") + 
  theme_classic(base_size = 17) + scale_color_manual(values = c("darkslategray4", "darkgoldenrod4")) + NoLegend() + xlim(values = c(-3,3))
ggsave("results/functional/volcano_3_exp_nk_tumor_only.png", width = 5.5, height = 4.5)

FeaturePlot(cml3_seurat_nk_expanded_tumor_no_tumor, features = c("TNFRSF18", "TNFRSF4", "TNFRSF9", "CRTAM", "HAVCR2", "TIGIT"), cols = c("gray90", "red3"), ncol = 3) & theme_void() & NoLegend() & theme(title = element_text(face = "italic")) #& labs(title = "")
ggsave("results/functional/umap_3_exp_nk_tumor_only_features.png", width = 4, height = 3)







## Only expanded
cml3_seurat_nk_expanded <- subset(cml3_seurat_nk, expanded == "Expanded NK")
cml3_seurat_nk_expanded <- cml3_seurat_nk_expanded %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml3_seurat_nk_expanded), nPCs = 11)
cml3_seurat_nk_expanded@reductions$pca@stdev %>% plot()
Idents(cml3_seurat_nk_expanded) <- cml3_seurat_nk_expanded$pvr_ko

Idents(cml3_seurat_nk_expanded) <- cml3_seurat_nk_expanded$my.demux

cml3_seurat_nk_expanded$my.demux

DimPlot(cml3_seurat_nk_expanded, cols = c("green4", "blue3", "red3", "green3")) 
ggsave("results/functional/hto3/umap_expanded_nk_cells.png", width = 5, height = 3)

Idents(cml3_seurat_nk_expanded) <- cml3_seurat_nk_expanded$pvr_ko
DimPlot(cml3_seurat_nk_expanded)

DimPlot(cml3_seurat_nk_expanded, group.by = "RNA_snn_res.1")


## Only expanded, with tumor
cml3_seurat_nk_expanded_tumor <- subset(cml3_seurat_nk, tumor == "Tumor")
cml3_seurat_nk_expanded_tumor <- subset(cml3_seurat_nk_expanded_tumor, expanded == "Expanded NK")
cml3_seurat_nk_expanded_tumor <- cml3_seurat_nk_expanded_tumor %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml3_seurat_nk_expanded_tumor), nPCs = 19)
cml3_seurat_nk_expanded_tumor@reductions$pca@stdev %>% plot()

Idents(cml3_seurat_nk_expanded_tumor) <- cml3_seurat_nk_expanded_tumor$my.demux
DimPlot(cml3_seurat_nk_expanded_tumor, cols = c("green4", "blue3", "green3")) 
ggsave("results/functional/hto3/umap_expanded_tumor_nk_cells.png", width = 5, height = 3)

Idents(cml3_seurat_nk_expanded_tumor) <- cml3_seurat_nk_expanded_tumor$my.demux
DimPlot(cml3_seurat_nk_expanded_tumor)

Idents(cml3_seurat_nk_expanded_tumor) <- cml3_seurat_nk_expanded_tumor$pvr_ko
DimPlot(cml3_seurat_nk_expanded_tumor)



## Only expanded, with tumor; PVR1.1 guide
cml3_seurat_nk_expanded_tumor_pvr1.1 <- subset(cml3_seurat_nk_expanded_tumor, my.demux %in% c("LAMA84-PVR1-T-E-NK-", "LAMA84-ctrl2-T-E-NK-"))
cml3_seurat_nk_expanded_tumor_pvr1.1 <- cml3_seurat_nk_expanded_tumor_pvr1.1 %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml3_seurat_nk_expanded_tumor_pvr1.1), nPCs = 16)
cml3_seurat_nk_expanded_tumor_pvr1.1@reductions$pca@stdev %>% plot()
Idents(cml3_seurat_nk_expanded_tumor_pvr1.1) <- cml3_seurat_nk_expanded_tumor_pvr1.1$pvr_ko

Idents(cml3_seurat_nk_expanded_tumor_pvr1.1) <- cml3_seurat_nk_expanded_tumor_pvr1.1$my.demux
DimPlot(cml3_seurat_nk_expanded_tumor_pvr1.1, cols = c("green4", "blue3", "green3")) 
ggsave("results/functional/hto3/umap_expanded_tumor_nk_cells_pvr1.png", width = 5, height = 3)

DimPlot(cml3_seurat_nk_expanded_tumor_pvr1.1)
ggsave("results/functional/hto3/umap_pvr1.png", width = 4, height = 3)

FeaturePlot(cml3_seurat_nk_expanded_tumor_pvr1.1, features = "TIGIT")

DimPlot(cml3_seurat_nk_expanded_tumor_pvr1.1, group.by = "RNA_snn_res.0.3")


cml3_seurat_nk_expanded_tumor_pvr1.1@meta.data %>% 
  group_by(cl=RNA_snn_res.0.3,pvr_ko) %>% tally() %>% mutate(prop=n/sum(n))

DimPlot(cml3_seurat_nk_expanded_tumor_pvr1.1, group.by = "Phase")
ggsave("results/functional/hto3/umap_pvr1.png", width = 4, height = 3)


Idents(cml3_seurat_nk_expanded_tumor_pvr1.1) <- cml3_seurat_nk_expanded_tumor_pvr1.1$pvr_ko
DimPlot(cml3_seurat_nk_expanded_tumor_pvr1.1)
ggsave("results/functional/hto3/umap_pvr1_2.png", width = 4, height = 3)

cml3_seurat_nk_expanded_tumor_pvr1.1_markers <- FindMarkers(cml3_seurat_nk_expanded_tumor_pvr1.1, logfc.threshold = 0.15, ident.1 = "PVR KO", ident.2 = "PVR WT", 
                                                              test.use = "t", max.cells.per.ident = 1e3) %>%
  filter(p_val_adj < 0.05) %>% add_rownames(var = "gene")

ggplot(data = cml3_seurat_nk_expanded_tumor_pvr1.1_markers, aes(avg_log2FC, -log10(p_val_adj), label = gene)) + geom_point() + 
  ggrepel::geom_text_repel(data = subset(cml3_seurat_nk_expanded_tumor_pvr1.1_markers, p_val_adj < 1e-10)) +
  ggpubr::theme_pubr()
ggsave("results/functional/hto3/volcano_exp_nk_pvr_ko_pvr_wt.png", width = 9, height = 9)



## Only expanded, with tumor; PVR1.2 guide
cml3_seurat_nk_expanded_tumor_pvr1.2 <- subset(cml3_seurat_nk_expanded_tumor, my.demux %in% c("LAMA84-PVR2-T-E-NK-", "LAMA84-ctrl2-T-E-NK-"))
cml3_seurat_nk_expanded_tumor_pvr1.2 <- cml3_seurat_nk_expanded_tumor_pvr1.2 %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml3_seurat_nk_expanded_tumor_pvr1.2), nPCs = 20)
cml3_seurat_nk_expanded_tumor_pvr1.2@reductions$pca@stdev %>% plot()
Idents(cml3_seurat_nk_expanded_tumor_pvr1.2) <- cml3_seurat_nk_expanded_tumor_pvr1.2$pvr_ko

Idents(cml3_seurat_nk_expanded_tumor_pvr1.2) <- cml3_seurat_nk_expanded_tumor_pvr1.2$my.demux
DimPlot(cml3_seurat_nk_expanded_tumor_pvr1.2, cols = c("blue3", "green3")) 
ggsave("results/functional/hto3/umap_expanded_tumor_nk_cells_pvr2.png", width = 5, height = 3)


Idents(cml3_seurat_nk_expanded_tumor_pvr1.2) <- cml3_seurat_nk_expanded_tumor_pvr1.2$my.demux
DimPlot(cml3_seurat_nk_expanded_tumor_pvr1.2)

Idents(cml3_seurat_nk_expanded_tumor_pvr1.2) <- cml3_seurat_nk_expanded_tumor_pvr1.2$pvr_ko
DimPlot(cml3_seurat_nk_expanded_tumor_pvr1.2)

cml3_seurat_nk_expanded_tumor_pvr1.2_markers <- FindMarkers(cml3_seurat_nk_expanded_tumor_pvr1.2, logfc.threshold = 0.15, ident.1 = "PVR KO", ident.2 = "PVR WT", 
                                                            test.use = "t", max.cells.per.ident = 1e3) %>%
  filter(p_val_adj < 0.05) %>% add_rownames(var = "gene")

ggplot(data = cml3_seurat_nk_expanded_tumor_pvr1.2_markers, aes(avg_log2FC, -log10(p_val_adj), label = gene)) + geom_point() + 
  ggrepel::geom_text_repel(data = subset(cml3_seurat_nk_expanded_tumor_pvr1.2_markers, p_val_adj < 1e-5)) +
  ggpubr::theme_pubr()
ggsave("results/functional/hto3/volcano_exp_nk_pvr2_ko_pvr_wt.png", width = 9, height = 9)


## Only non-expanded
cml3_seurat_nk_nonexpanded <- subset(cml3_seurat_nk, expanded != "Expanded NK")
cml3_seurat_nk_nonexpanded <- cml3_seurat_nk_nonexpanded %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml3_seurat_nk_nonexpanded), nPCs = 10)
cml3_seurat_nk_nonexpanded@reductions$pca@stdev %>% plot()
cml3_seurat_nk_nonexpanded <- cml3_seurat_nk_nonexpanded %>% getClustering()
Idents(cml3_seurat_nk_nonexpanded) <- cml3_seurat_nk_nonexpanded$RNA_snn_res.0.1

a <- DimPlot(cml3_seurat_nk_nonexpanded, group.by = "my.demux", cols = getPalette(10))
b <- DimPlot(cml3_seurat_nk_nonexpanded, group.by = "covariate", cols = getPalette(5))
c <- DimPlot(cml3_seurat_nk_nonexpanded, group.by = "expanded", cols = getPalette3(4))
d <- DimPlot(cml3_seurat_nk_nonexpanded, group.by = "pvr_ko", cols = getPalette3(4))
e <- DimPlot(cml3_seurat_nk_nonexpanded, group.by = "Phase", cols = getPalette3(4))
f <- DimPlot(cml3_seurat_nk_nonexpanded, group.by = "RNA_snn_res.0.1", cols = getPalette3(4))

a + b + c + d + e + f
ggsave("results/functional/hto3/umap_meta_nonexp_nk.png", width = 12, height = 5)

DimPlot(cml3_seurat_nk_nonexpanded, group.by = "my.demux", cols = c("blue3","green3","green4"))
ggsave("results/functional/hto3/umap_meta_nonexp_nk_tot.png", width = 5, height = 3)





## Only non-nonexpanded with PVR1 guide
cml3_seurat_nk_nonexpanded_tumor_pvr1.1 <- subset(cml3_seurat_nk_nonexpanded, my.demux %in% c("LAMA84-PVR1-T-NE-NK-", "LAMA84-ctrl2-T-NE-NK-"))
cml3_seurat_nk_nonexpanded_tumor_pvr1.1 <- cml3_seurat_nk_nonexpanded_tumor_pvr1.1 %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml3_seurat_nk_nonexpanded_tumor_pvr1.1), nPCs = 10)
cml3_seurat_nk_nonexpanded_tumor_pvr1.1@reductions$pca@stdev %>% plot()
Idents(cml3_seurat_nk_nonexpanded_tumor_pvr1.1) <- cml3_seurat_nk_nonexpanded_tumor_pvr1.1$pvr_ko

Idents(cml3_seurat_nk_nonexpanded_tumor_pvr1.1) <- cml3_seurat_nk_nonexpanded_tumor_pvr1.1$my.demux
DimPlot(cml3_seurat_nk_nonexpanded_tumor_pvr1.1)
ggsave("results/functional/hto3/umap_pvr1_nonexp.png", width = 4, height = 3)

Idents(cml3_seurat_nk_nonexpanded_tumor_pvr1.1) <- cml3_seurat_nk_nonexpanded_tumor_pvr1.1$pvr_ko
DimPlot(cml3_seurat_nk_nonexpanded_tumor_pvr1.1)
ggsave("results/functional/hto3/umap_pvr1_n2_onexp.png", width = 4, height = 3)

cml3_seurat_nk_nonexpanded_tumor_pvr1.1_markers <- FindMarkers(cml3_seurat_nk_nonexpanded_tumor_pvr1.1, logfc.threshold = 0.15, ident.1 = "PVR KO", ident.2 = "PVR WT", 
                                                            test.use = "t", max.cells.per.ident = 1e3) %>%
  filter(p_val_adj < 0.05) %>% add_rownames(var = "gene")

ggplot(data = cml3_seurat_nk_nonexpanded_tumor_pvr1.1_markers, aes(avg_log2FC, -log10(p_val_adj), label = gene)) + geom_point() + 
  ggrepel::geom_text_repel(data = subset(cml3_seurat_nk_nonexpanded_tumor_pvr1.1_markers, p_val_adj < 1e-10)) +
  ggpubr::theme_pubr()
ggsave("results/functional/hto3/volcano_nonexp_nk_pvr_ko_pvr_wt.png", width = 9, height = 9)



DimPlot(cml3_seurat_nk_nonexpanded_tumor_pvr1.1, group.by = "my.demux", cols = c("blue3","green3"))
ggsave("results/functional/hto3/umap_meta_nonexp_nk_tot_pvr.png", width = 5, height = 3)





## Only non-nonexpanded with PVR1.2 guide
cml3_seurat_nk_nonexpanded_tumor_pvr1.2 <- subset(cml3_seurat_nk_nonexpanded, my.demux %in% c("LAMA84-PVR2-T-NE-NK-", "LAMA84-ctrl2-T-NE-NK-"))
cml3_seurat_nk_nonexpanded_tumor_pvr1.2 <- cml3_seurat_nk_nonexpanded_tumor_pvr1.2 %>% preprocessSeuratCellCycle(cells.to.use = colnames(cml3_seurat_nk_nonexpanded_tumor_pvr1.2), nPCs = 10)
cml3_seurat_nk_nonexpanded_tumor_pvr1.2@reductions$pca@stdev %>% plot()
Idents(cml3_seurat_nk_nonexpanded_tumor_pvr1.2) <- cml3_seurat_nk_nonexpanded_tumor_pvr1.2$pvr_ko

Idents(cml3_seurat_nk_nonexpanded_tumor_pvr1.2) <- cml3_seurat_nk_nonexpanded_tumor_pvr1.2$my.demux
DimPlot(cml3_seurat_nk_nonexpanded_tumor_pvr1.2)
ggsave("results/functional/hto3/umap_pvr1.2_nonexp.png", width = 4, height = 3)

Idents(cml3_seurat_nk_nonexpanded_tumor_pvr1.2) <- cml3_seurat_nk_nonexpanded_tumor_pvr1.2$pvr_ko
DimPlot(cml3_seurat_nk_nonexpanded_tumor_pvr1.2)
ggsave("results/functional/hto3/umap_pvr.21_n2_onexp.png", width = 4, height = 3)

DimPlot(cml3_seurat_nk_nonexpanded_tumor_pvr1.2, group.by = "my.demux", cols = c("blue3","green4"))
ggsave("results/functional/hto3/umap_meta_nonexp_nk_tot_pvr2.png", width = 5, height = 3)

cml3_seurat_nk_nonexpanded_tumor_pvr1.2_markers <- FindMarkers(cml3_seurat_nk_nonexpanded_tumor_pvr1.2, logfc.threshold = 0.15, ident.1 = "PVR KO", ident.2 = "PVR WT", 
                                                               test.use = "t", max.cells.per.ident = 1e3) %>%
  filter(p_val_adj < 0.05) %>% add_rownames(var = "gene")

ggplot(data = cml3_seurat_nk_nonexpanded_tumor_pvr1.2_markers, aes(avg_log2FC, -log10(p_val_adj), label = gene)) + geom_point() + 
  ggrepel::geom_text_repel(data = subset(cml3_seurat_nk_nonexpanded_tumor_pvr1.2_markers, p_val_adj < 1e-10)) +
  ggpubr::theme_pubr()
ggsave("results/functional/hto3/volcano_nonexp_nk_pvr_ko_pvr1.2_wt.png", width = 9, height = 9)

