
## Run cellphone db. Make subselection of 1000 cells from responders at time 0 and time 1
dir.create("results/cellphonedb/", showWarnings = F)
dir.create("results/cellphonedb/input_files", showWarnings = F)
dir.create("results/cellphonedb/analysis", showWarnings = F)

sample_n = 100

################################ Init cellphone input files



### no relapse
data.frame(cml_seurat@meta.data, cluster = Idents(cml_seurat)) %>% filter(relapse == "None", timepoint == "baseline") %>%
  initCellphonedb(seurat_object = cml_seurat, name = "no_0m", sample_n = sample_n)

data.frame(cml_seurat@meta.data, cluster = Idents(cml_seurat)) %>% filter(relapse == "None", timepoint == "6m") %>%
  initCellphonedb(seurat_object = cml_seurat, name = "no_6m", sample_n = sample_n)

data.frame(cml_seurat@meta.data, cluster = Idents(cml_seurat)) %>% filter(relapse == "None", timepoint == "12m") %>%
  initCellphonedb(seurat_object = cml_seurat, name = "no_12m", sample_n = sample_n)



### slow relapse
data.frame(cml_seurat@meta.data, cluster = Idents(cml_seurat)) %>% filter(relapse == "Slow", timepoint == "baseline") %>%
  initCellphonedb(seurat_object = cml_seurat, name = "slow_0m", sample_n = sample_n)

data.frame(cml_seurat@meta.data, cluster = Idents(cml_seurat)) %>% filter(relapse == "Slow", timepoint == "6m") %>%
  initCellphonedb(seurat_object = cml_seurat, name = "slow_6m", sample_n = sample_n)

data.frame(cml_seurat@meta.data, cluster = Idents(cml_seurat)) %>% filter(relapse == "Slow", timepoint == "relapse") %>%
  initCellphonedb(seurat_object = cml_seurat, name = "slow_relapse", sample_n = sample_n)



### fast relapse
data.frame(cml_seurat@meta.data, cluster = Idents(cml_seurat)) %>% filter(relapse == "Fast", timepoint == "baseline") %>%
  initCellphonedb(seurat_object = cml_seurat, name = "fast_0m", sample_n = sample_n)

data.frame(cml_seurat@meta.data, cluster = Idents(cml_seurat)) %>% filter(relapse == "Fast", timepoint == "relapse") %>%
  initCellphonedb(seurat_object = cml_seurat, name = "fast_relapse", sample_n = sample_n)


### disease controllers
data.frame(cml_seurat@meta.data, cluster = Idents(cml_seurat)) %>% filter(relapse != "Fast", timepoint == "baseline") %>%
  initCellphonedb(seurat_object = cml_seurat, name = "disease_ctrl_0m", sample_n = sample_n)

### disease controllers
data.frame(cml_seurat@meta.data, cluster = Idents(cml_seurat)) %>% filter(relapse != "Fast", timepoint == "baseline") %>%
  initCellphonedb(seurat_object = cml_seurat, name = "disease_ctrl_0m", sample_n = sample_n)


### tki controllers
df <- data.frame(cml_tki_seurat@meta.data, cluster = as.character(cml_tki_seurat$public_clusters)) %>% filter(project == "CML dg")
df$cluster <- as.character(df$cluster)
df$cluster <- df$public_clusters
df %>% initCellphonedb(seurat_object = cml_tki_seurat, name = "cml_dg_0m", sample_n = sample_n)

data.frame(cml_tki_seurat@meta.data, cluster = Idents(cml_tki_seurat)) %>% filter(project == "CML on TKI") %>%
  initCellphonedb(seurat_object = cml_tki_seurat, name = "cml_tki_0m", sample_n = sample_n)

data.frame(cml_tki_seurat@meta.data, cluster = Idents(cml_tki_seurat)) %>% filter(project == "CML dg") %>%
  initCellphonedb(seurat_object = cml_tki_seurat, name = "cml_dg_0m", sample_n = sample_n)

################################ Analyze

list.dirs("results/cellphonedb/out/")

sigf_means_no_relapse_0m  <- fread("results/cellphonedb/out/no_0m/significant_means.txt")
sigf_means_no_relapse_6m  <- fread("results/cellphonedb/out/no_6m/significant_means.txt")
sigf_means_no_relapse_12m <- fread("results/cellphonedb/out/no_12m/significant_means.txt")

sigf_means_slow_relapse_0m  <- fread("results/cellphonedb/out/slow_0m/significant_means.txt")
sigf_means_slow_relapse_6m  <- fread("results/cellphonedb/out/slow_6m/significant_means.txt")
sigf_means_slow_relapse_rp  <- fread("results/cellphonedb/out/slow_relapse//significant_means.txt")

sigf_means_fast_relapse_0m  <- fread("results/cellphonedb/out/fast_0m/significant_means.txt")
sigf_means_fast_relapse_rp  <- fread("results/cellphonedb/out/fast_relapse/significant_means.txt")

sigf_means_disease_ctrl_0m  <- fread("results/cellphonedb/out/disease_ctrl_0m/significant_means.txt")

sigf_means_disease_cml_dg_0m  <- fread("results/cellphonedb/out/cml_dg_0m/significant_means.txt")
sigf_means_disease_cml_tki_0m <- fread("results/cellphonedb/out/on_tki_0m//significant_means.txt")



getNewClusters <- function(clusters){
  clusters %>% extractClusterNumber() %>% as.numeric() %>% as.factor() %>% getClusterPhenotypesCml()
}


## How many interactions are there?
sigf_means_no_relapse_0m_pairs    <- getPairAmounts(sigf_means_no_relapse_0m)  %>% filterRedundantPairs() %>% mutate(cluster1 = getNewClusters(cluster1), cluster2 = getNewClusters(cluster2))
sigf_means_no_relapse_6m_pairs    <- getPairAmounts(sigf_means_no_relapse_6m)  %>% filterRedundantPairs() %>% mutate(cluster1 = getNewClusters(cluster1), cluster2 = getNewClusters(cluster2))
sigf_means_no_relapse_12m_pairs   <- getPairAmounts(sigf_means_no_relapse_12m) %>% filterRedundantPairs() %>% mutate(cluster1 = getNewClusters(cluster1), cluster2 = getNewClusters(cluster2))

sigf_means_slow_relapse_0m_pairs  <- getPairAmounts(sigf_means_slow_relapse_0m)  %>% filterRedundantPairs() %>% mutate(cluster1 = getNewClusters(cluster1), cluster2 = getNewClusters(cluster2))
sigf_means_slow_relapse_6m_pairs  <- getPairAmounts(sigf_means_slow_relapse_6m)  %>% filterRedundantPairs() %>% mutate(cluster1 = getNewClusters(cluster1), cluster2 = getNewClusters(cluster2))
sigf_means_slow_relapse_12m_pairs <- getPairAmounts(sigf_means_slow_relapse_rp)  %>% filterRedundantPairs() %>% mutate(cluster1 = getNewClusters(cluster1), cluster2 = getNewClusters(cluster2))

sigf_means_fast_relapse_0m_pairs  <- getPairAmounts(sigf_means_fast_relapse_0m) %>% filterRedundantPairs() %>% mutate(cluster1 = getNewClusters(cluster1), cluster2 = getNewClusters(cluster2))
sigf_means_fast_relapse_12m_pairs <- getPairAmounts(sigf_means_fast_relapse_rp) %>% filterRedundantPairs() %>% mutate(cluster1 = getNewClusters(cluster1), cluster2 = getNewClusters(cluster2))

sigf_means_disease_ctrl_0m_pairs  <- getPairAmounts(sigf_means_disease_ctrl_0m) %>% filterRedundantPairs() %>% mutate(cluster1 = getNewClusters(cluster1), cluster2 = getNewClusters(cluster2))

sigf_means_disease_cml_dg_0m_pairs  <- getPairAmounts(sigf_means_disease_cml_dg_0m) %>% filterRedundantPairs() %>% mutate(cluster1 = getNewClusters(cluster1), cluster2 = getNewClusters(cluster2))
sigf_means_disease_cml_tki_0m_pairs <- getPairAmounts(sigf_means_disease_cml_tki_0m) %>% filterRedundantPairs() %>% mutate(cluster1 = getNewClusters(cluster1), cluster2 = getNewClusters(cluster2))



## Heatmaps
p <- plotDiffHeatmap(df1 = sigf_means_fast_relapse_0m_pairs, df2 = sigf_means_disease_ctrl_0m_pairs)
ggsave(plot = print(p), "results/cellphonedb/analysis/heatmap_relapse_vs_ctrl_0m.png", width = 8, height = 6)


p <- plotDiffHeatmap(df1 = sigf_means_fast_relapse_0m_pairs, df2 = sigf_means_no_relapse_0m_pairs)
ggsave(plot = print(p), "results/cellphonedb/analysis/heatmap_relapse_vs_no_0m.png", width = 8, height = 6)

p <- plotDiffHeatmap(df1 = sigf_means_slow_relapse_0m_pairs, df2 = sigf_means_no_relapse_0m_pairs)
ggsave(plot = print(p), "results/cellphonedb/analysis/heatmap_slow_relapse_vs_no_0m.png", width = 8, height = 6)


p <- plotDiffHeatmap(df1 = sigf_means_no_relapse_0m_pairs, df2 = sigf_means_no_relapse_6m_pairs)
ggsave(plot = print(p), "results/cellphonedb/analysis/heatmap_no_relapse_0m_vs_no_6m.png", width = 8, height = 6)

p <- plotDiffHeatmap(df1 = sigf_means_slow_relapse_0m_pairs, df2 = sigf_means_slow_relapse_6m_pairs)
ggsave(plot = print(p), "results/cellphonedb/analysis/heatmap_slow_relapse_0m_vs_no_6m.png", width = 8, height = 6)

p <- plotDiffHeatmap(df1 = sigf_means_fast_relapse_0m_pairs, df2 = sigf_means_fast_relapse_12m_pairs)
ggsave(plot = print(p), "results/cellphonedb/analysis/heatmap_fast_relapse_0m_vs_no_6m.png", width = 8, height = 6)


p <- plotDiffHeatmap(df1 = sigf_means_no_relapse_0m_pairs, df2 = sigf_means_no_relapse_6m_pairs)
ggsave(plot = print(p), "results/cellphonedb/analysis/heatmap_slow_relapse_vs_no_0m.png", width = 8, height = 6)


p <- plotDiffHeatmap(df1 = sigf_means_disease_cml_dg_0m_pairs, df2 = sigf_means_disease_cml_tki_0m_pairs)
ggsave(plot = print(p), "results/cellphonedb/analysis/heatmap_dg_vs_tki_0m.png", width = 8, height = 6)





## Heatmaps
p <- plotDiffHeatmap(df1 = sigf_means_fast_relapse_0m_pairs, df2 = sigf_means_disease_ctrl_0m_pairs)
ggsave(plot = print(p), "results/cellphonedb/analysis/heatmap_relapse_vs_ctrl_0m.png", width = 8, height = 6)

p <- plotDiffHeatmap(df1 = sigf_means_disease_cml_dg_0m_pairs, df2 = sigf_means_disease_cml_tki_0m_pairs)
ggsave(plot = print(p), "results/cellphonedb/analysis/heatmap_dg_vs_tki_0m.png", width = 8, height = 6)


interactions_baseline <- sigf_means_no_relapse_0m_pairs %>% left_join(dplyr::select(sigf_means_slow_relapse_0m_pairs, pair, value), by = "pair") %>% left_join(dplyr::select(sigf_means_fast_relapse_0m_pairs, pair, value), by = "pair")

interactions_baseline_var <- interactions_baseline %>% dplyr::select(value.x,value.y,value)
interactions_baseline_var[is.na(interactions_baseline_var)] <- 0
RowVar(interactions_baseline_var)


interactions_baseline_tot <- interactions_baseline %>% group_by(pair) %>% summarise(tot = sum(value))
interactions_baseline_tot <- interactions_baseline_tot[order(interactions_baseline_var %>% RowVar()), ]

interactions_baseline <- interactions_baseline %>% dplyr::select(pair,value.x,value.y,value) %>% melt(id = "pair")

interactions_baseline %>% left_join(interactions_baseline_tot) %>%
  ggplot(aes(reorder(pair, tot), value, color = variable)) + geom_point() + ggpubr::rotate_x_text(angle = 45) + geom_path(aes(group=pair))
ggsave("results/manuscript/tki_stop/path_pairs.pdf", width = 24, height = 5)


RowVar <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}

interactions_baseline




## At baseline, different fates
p <- cast(sigf_means_no_relapse_0m_pairs, cluster1 ~ cluster2, sum)
rownames(p) <- colnames(p)[-1]

col_mat = rand_color(length(p[,-1]), transparency = 0.5)
dim(col_mat) = dim(p[,-1])  # to make sure it is a matrix

pdf("results/manuscript/tki_stop/circos_cellphone_no_0m.pdf", width = 8, height = 6)
chordDiagram(as.data.frame(p[,-1]), transparency = 0.5, col = col_mat) # + scale_fill_manual(values = getPalette5(23))
circos.clear()
dev.off()

p <- p %>% pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(seq(0, 60, 0.6))), breaks = seq(0, 60, 0.6))
ggsave(plot = p, "results/manuscript/tki_stop//heatmap_cellphone_no_0m.png", width = 8, height = 6)



p <- cast(sigf_means_slow_relapse_0m_pairs, cluster1 ~ cluster2, sum)
rownames(p) <- colnames(p)[-1]

pdf("results/manuscript/tki_stop/circos_cellphone_slow_0m.pdf", width = 8, height = 6)
chordDiagram(as.data.frame(p[,-1]), transparency = 0.5, col = col_mat) # + scale_fill_manual(values = getPalette5(23))
circos.clear()
dev.off()

p <- p %>% pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(seq(0, 60, 0.6))), breaks = seq(0, 60, 0.6))
ggsave(plot = p, "results/manuscript/tki_stop//heatmap_cellphone_slow_0m.png", width = 8, height = 6)




p <- cast(sigf_means_fast_relapse_0m_pairs, cluster1 ~ cluster2, sum)
rownames(p) <- colnames(p)[-1]

pdf("results/manuscript/tki_stop/circos_cellphone_fast_0m.pdf", width = 8, height = 6)
chordDiagram(as.data.frame(p[,-1]), transparency = 0.5, col = col_mat) # + scale_fill_manual(values = getPalette5(23))
circos.clear()
dev.off()


p <- p %>% pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(seq(0, 60, 0.6))), breaks = seq(0, 60, 0.6))
ggsave(plot = p, "results/manuscript/tki_stop//heatmap_cellphone_fast_0m.png", width = 8, height = 6)


## At follow-up, different fates
p <- cast(sigf_means_no_relapse_6m_pairs, cluster1 ~ cluster2, sum)
rownames(p) <- colnames(p)[-1]
p <- p %>% pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(seq(0, 60, 0.6))), breaks = seq(0, 60, 0.6))
ggsave(plot = p, "results/manuscript/tki_stop//heatmap_cellphone_no_6m.png", width = 8, height = 6)

p <- cast(sigf_means_slow_relapse_6m_pairs, cluster1 ~ cluster2, sum)
rownames(p) <- colnames(p)[-1]
p <- p %>% pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(seq(0, 60, 0.6))), breaks = seq(0, 60, 0.6))
ggsave(plot = p, "results/manuscript/tki_stop//heatmap_cellphone_slow_6m.png", width = 8, height = 6)

p <- cast(sigf_means_fast_relapse_12m_pairs, cluster1 ~ cluster2, sum)
rownames(p) <- colnames(p)[-1]
p <- p %>% pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(seq(0, 60, 0.6))), breaks = seq(0, 60, 0.6))
ggsave(plot = p, "results/manuscript/tki_stop//heatmap_cellphone_fast_6m.png", width = 8, height = 6)







## Dg vs on-TKI
p <- cast(sigf_means_disease_cml_dg_0m_pairs, cluster1 ~ cluster2, sum)
rownames(p) <- colnames(p)[-1]
p <- p %>% pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(seq(0, 60, 0.6))), breaks = seq(0, 60, 0.6))
ggsave(plot = p, "results/manuscript/cml_tki/heatmap_cellphone_dg.png", width = 8, height = 6)

p <- cast(sigf_means_disease_cml_tki_0m_pairs, cluster1 ~ cluster2, sum)
rownames(p) <- colnames(p)[-1]
p <- p %>% pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(seq(0, 60, 0.6))), breaks = seq(0, 60, 0.6))
ggsave(plot = p, "results/manuscript/cml_tki/heatmap_cellphone_on_tki.png", width = 8, height = 6)




## Dg 0m
p <- cast(sigf_means_disease_cml_dg_0m_pairs, cluster1 ~ cluster2, sum)
rownames(p) <- colnames(p)[-1]
p <- p %>% pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(seq(0, 60, 0.6))), breaks = seq(0, 60, 0.6))
ggsave(plot = print(p), "results/manuscript/relapse/heatmap_cml_dg_0m.pdf", width = 8, height = 6)

## fast relapse
p <- cast(sigf_means_fast_relapse_12m_pairs, cluster1 ~ cluster2, sum)
rownames(p) <- colnames(p)[-1]
p <- p %>% pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(seq(0, 60, 0.6))), breaks = seq(0, 60, 0.6))
ggsave(plot = p, "results/manuscript/relapse/heatmap_fast_relapse.pdf", width = 8, height = 6)

## late relapse
p <- cast(sigf_means_slow_relapse_12m_pairs, cluster1 ~ cluster2, sum)
rownames(p) <- colnames(p)[-1]
p <- p %>% pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(seq(0, 60, 0.6))), breaks = seq(0, 60, 0.6))
ggsave(plot = p, "results/manuscript/relapse/heatmap_slow_relapse.pdf", width = 8, height = 6)


## late relapse, treg interactions
slow_treg_interactome <- sigf_means_slow_relapse_rp %>% dplyr::select(c(c("id_cp_interaction","interacting_pair","partner_a","partner_b","gene_a","gene_b","secreted","receptor_a","receptor_b","annotation_strategy","is_integrin","rank"), grep("Treg", colnames(sigf_means_slow_relapse_rp), value = T)))
slow_treg_interactome

slow_treg_interactome_mtx <- slow_treg_interactome %>% dplyr::select(-c("id_cp_interaction","interacting_pair","partner_a","partner_b","gene_a","gene_b","secreted","receptor_a","receptor_b","annotation_strategy","is_integrin","rank"))
rows.to.keep <- rowSums(slow_treg_interactome_mtx, na.rm = T) > 0
slow_treg_interactome <- slow_treg_interactome[rows.to.keep, ]
slow_treg_interactome_mtx <- slow_treg_interactome_mtx[rows.to.keep, ]
rownames(slow_treg_interactome_mtx) <- slow_treg_interactome$interacting_pair

write.table(slow_treg_interactome, "results/cellphonedb/analysis/slow_treg_interactome.txt", sep = "\t", quote = F, row.names = F)
slow_treg_interactome_inh <- slow_treg_interactome %>% filter(gene_a %in% inhibitory_long | gene_b %in% inhibitory_long)
write.table(slow_treg_interactome_inh, "results/cellphonedb/analysis/slow_treg_interactome_inhibitory.txt", sep = "\t", quote = F, row.names = F)

slow_treg_interactome_inh


p <- cast(sigf_means_disease_ctrl_0m_pairs, cluster1 ~ cluster2, sum) %>% pheatmap::pheatmap()
ggsave(plot = p, "results/cellphonedb/analysis/heatmap_sigf_tot_0m.png", width = 8, height = 6)


## Which phenotypes gets/loses sigf amount of interactions?
sigf_r    <- testSigfPairs(sigf_means_r_0m_pairs, sigf_means_r_1m_pairs) %>% filter(p.value < 0.05)
write.table(sigf_r, "results/cellphonedb/analysis/de_clusters_r_2v1.txt", sep = "\t", quote = F, row.names = F)

sigf_r    <- testSigfPairs(sigf_means_r_1m_pairs, sigf_means_r_3m_pairs) %>% filter(p.value < 0.05)
write.table(sigf_r, "results/cellphonedb/analysis/de_clusters_r_3v2.txt", sep = "\t", quote = F, row.names = F)


sigf_n    <- testSigfPairs(sigf_means_n_0m_pairs, sigf_means_n_1m_pairs) %>% filter(p.value < 0.05)
write.table(sigf_n, "results/cellphonedb/analysis/de_clusters_n_2v1.txt", sep = "\t", quote = F, row.names = F)

sigf_n    <- testSigfPairs(sigf_means_n_1m_pairs, sigf_means_n_3m_pairs) %>% filter(p.value < 0.05)
write.table(sigf_n, "results/cellphonedb/analysis/de_clusters_n_3v2.txt", sep = "\t", quote = F, row.names = F)


sigf_diff_0m <- testSigfPairs(sigf_means_r_0m_pairs, sigf_means_n_0m_pairs, paired = F) %>% filter(p.value < 0.05)
write.table(sigf_diff_0m, "results/cellphonedb/analysis/de_clusters_r_vs_n_0m.txt", sep = "\t", quote = F, row.names = F)

sigf_diff_1m <- testSigfPairs(sigf_means_r_1m_pairs, sigf_means_n_1m_pairs, paired = F) %>% filter(p.value < 0.05)
write.table(sigf_diff_1m, "results/cellphonedb/analysis/de_clusters_r_vs_n_1m.txt", sep = "\t", quote = F, row.names = F)

sigf_diff_3m <- testSigfPairs(sigf_means_r_3m_pairs, sigf_means_n_3m_pairs, paired = F) %>% filter(p.value < 0.05)
write.table(sigf_diff_3m, "results/cellphonedb/analysis/de_clusters_r_vs_n_3m.txt", sep = "\t", quote = F, row.names = F)


## Visualize Which phenotypes gets/loses sigf amount of interactions?
p <- plotPairs(sigf_means_tot_0m_pairs, sigf_means_tot_1m_pairs)
ggsave(plot = p, "results/cellphonedb/analysis/lineplot_sigf_tot_2v1.png", width = 12, height = 10)

p <- plotPairs(sigf_means_r_0m_pairs, sigf_means_r_1m_pairs)
ggsave(plot = p, "results/cellphonedb/analysis/lineplot_sigf_r_2v1.png", width = 12, height = 10)

p <- plotPairs(sigf_means_n_0m_pairs, sigf_means_n_1m_pairs)
ggsave(plot = p, "results/cellphonedb/analysis/lineplot_sigf_n_2v1.png", width = 12, height = 10)



p <- plotNewPairs(sigf_means_tot_0m_pairs, sigf_means_tot_1m_pairs)
ggsave(plot = p, "results/cellphonedb/analysis/bar_new_sigf_tot_2v1.png", width = 12, height = 10)

p <- plotNewPairs(sigf_means_r_0m_pairs, sigf_means_r_1m_pairs)
ggsave(plot = p, "results/cellphonedb/analysis/bar_new_sigf_r_2v1.png", width = 12, height = 10)

p <- plotNewPairs(sigf_means_n_0m_pairs, sigf_means_n_1m_pairs)
ggsave(plot = p, "results/cellphonedb/analysis/bar_new_sigf_n_2v1.png", width = 12, height = 10)

r_new_pairs_2v1 <- getNewPairs(df1 = sigf_means_r_0m_pairs, df2 = sigf_means_r_1m_pairs) %>% mutate(overall = "R")
n_new_pairs_2v1 <- getNewPairs(df1 = sigf_means_n_0m_pairs, df2 = sigf_means_n_1m_pairs) %>% mutate(overall = "N")


p <- rbind(r_new_pairs_2v1, n_new_pairs_2v1) %>%
  ggplot(aes(overall,value, fill = overall)) + geom_boxplot(outlier.shape = 21) +
  facet_wrap(~cluster2.x, scales = "free") +
  ggsignif::geom_signif(comparisons = list(c("N", "R"))) + scale_fill_manual(values = c("lightgrey", "salmon"))
ggsave(plot = p, "results/cellphonedb/analysis/boxplot_sigf_r_v_n_2v1.png", width = 15, height = 13)




p <- plotPairs(sigf_means_tot_1m_pairs, sigf_means_tot_3m_pairs)
ggsave(plot = p, "results/cellphonedb/analysis/lineplot_sigf_tot_3v2.png", width = 12, height = 10)

p <- plotPairs(sigf_means_r_1m_pairs, sigf_means_r_3m_pairs)
ggsave(plot = p, "results/cellphonedb/analysis/lineplot_sigf_r_3v2.png", width = 12, height = 10)

p <- plotPairs(sigf_means_n_1m_pairs, sigf_means_n_3m_pairs)
ggsave(plot = p, "results/cellphonedb/analysis/lineplot_sigf_n_3v2.png", width = 12, height = 10)



p <- plotNewPairs(sigf_means_tot_1m_pairs, sigf_means_tot_3m_pairs)
ggsave(plot = p, "results/cellphonedb/analysis/bar_new_sigf_tot_3v2.png", width = 12, height = 10)

p <- plotNewPairs(sigf_means_r_1m_pairs, sigf_means_r_3m_pairs)
ggsave(plot = p, "results/cellphonedb/analysis/bar_new_sigf_r_3v2.png", width = 12, height = 10)

p <- plotNewPairs(sigf_means_n_1m_pairs, sigf_means_n_3m_pairs)
ggsave(plot = p, "results/cellphonedb/analysis/bar_new_sigf_n_3v2.png", width = 12, height = 10)

r_new_pairs_3v2 <- getNewPairs(df1 = sigf_means_r_1m_pairs, df2 = sigf_means_r_3m_pairs) %>% mutate(overall = "R")
n_new_pairs_3v2 <- getNewPairs(df1 = sigf_means_n_1m_pairs, df2 = sigf_means_n_3m_pairs) %>% mutate(overall = "N")

p <- rbind(r_new_pairs_3v2, n_new_pairs_3v2) %>%
  ggplot(aes(overall,value, fill = overall)) + geom_boxplot(outlier.shape = 21) +
  facet_wrap(~cluster2.x, scales = "free") +
  ggsignif::geom_signif(comparisons = list(c("N", "R"))) + scale_fill_manual(values = c("lightgrey", "salmon"))
ggsave(plot = p, "results/cellphonedb/analysis/boxplot_sigf_r_v_n_3v2.png", width = 15, height = 13)








#### Heatmaps
p <- cast(sigf_means_tot_0m_pairs, cluster1 ~ cluster2, sum) %>% pheatmap::pheatmap()
ggsave(plot = p, "results/cellphonedb/analysis/heatmap_sigf_tot_0m.png", width = 8, height = 6)

p <- cast(sigf_means_tot_1m_pairs, cluster1 ~ cluster2, sum) %>% pheatmap::pheatmap()
ggsave(plot = p, "results/cellphonedb/analysis/heatmap_sigf_tot_1m.png", width = 8, height = 6)

p <- cast(sigf_means_tot_3m_pairs, cluster1 ~ cluster2, sum) %>% pheatmap::pheatmap()
ggsave(plot = p, "results/cellphonedb/analysis/heatmap_sigf_tot_3m.png", width = 8, height = 6)




p <- cast(sigf_means_r_0m_pairs, cluster1 ~ cluster2, sum) %>% pheatmap::pheatmap()
ggsave(plot = p, "results/cellphonedb/analysis/heatmap_sigf_r_0m.png", width = 8, height = 6)

p <- cast(sigf_means_r_1m_pairs, cluster1 ~ cluster2, sum) %>% pheatmap::pheatmap()
ggsave(plot = p, "results/cellphonedb/analysis/heatmap_sigf_r_1m.png", width = 8, height = 6)

p <- cast(sigf_means_r_3m_pairs, cluster1 ~ cluster2, sum) %>% pheatmap::pheatmap()
ggsave(plot = p, "results/cellphonedb/analysis/heatmap_sigf_r_3m.png", width = 8, height = 6)




p <- cast(sigf_means_n_0m_pairs, cluster1 ~ cluster2, sum) %>% pheatmap::pheatmap()
ggsave(plot = p, "results/cellphonedb/analysis/heatmap_sigf_n_0m.png", width = 8, height = 6)

p <- cast(sigf_means_n_1m_pairs, cluster1 ~ cluster2, sum) %>% pheatmap::pheatmap()
ggsave(plot = p, "results/cellphonedb/analysis/heatmap_sigf_n_1m.png", width = 8, height = 6)

p <- cast(sigf_means_n_3m_pairs, cluster1 ~ cluster2, sum) %>% pheatmap::pheatmap()
ggsave(plot = p, "results/cellphonedb/analysis/heatmap_sigf_n_3m.png", width = 8, height = 6)







p <- plotDiffHeatmap(df1 = sigf_means_tot_0m_pairs, df2 = sigf_means_tot_1m_pairs)
ggsave(plot = p, "results/cellphonedb/analysis/heatmap_sigf_tot_2v1.png", width = 8, height = 6)

p <- plotDiffHeatmap(df1 = sigf_means_r_1m_pairs, df2 = sigf_means_r_3m_pairs)
ggsave(plot = p, "results/cellphonedb/analysis/heatmap_sigf_tot_3v2.png", width = 8, height = 6)



p <- plotDiffHeatmap(df1 = sigf_means_r_0m_pairs, df2 = sigf_means_r_1m_pairs)
ggsave(plot = p, "results/cellphonedb/analysis/heatmap_sigf_r_2v1.png", width = 8, height = 6)

p <- plotDiffHeatmap(df1 = sigf_means_r_1m_pairs, df2 = sigf_means_r_3m_pairs)
ggsave(plot = p, "results/cellphonedb/analysis/heatmap_sigf_r_3v2.png", width = 8, height = 6)



p <- plotDiffHeatmap(df1 = sigf_means_n_0m_pairs, sigf_means_n_1m_pairs)
ggsave(plot = p, "results/cellphonedb/analysis/heatmap_sigf_n_2v1.png", width = 8, height = 6)

p <- plotDiffHeatmap(df1 = sigf_means_n_1m_pairs, sigf_means_n_3m_pairs)
ggsave(plot = p, "results/cellphonedb/analysis/heatmap_sigf_n_3v2.png", width = 8, height = 6)








## What are the sigf interactions?
sigf_means_r_0m_sigf_pairs <- getSigfPairs(sigf_means_r_0m) %>% mutate(cluster1 = getClusterPhenotypes(cluster1), cluster2 = getClusterPhenotypes(cluster2))
sigf_means_r_1m_sigf_pairs <- getSigfPairs(sigf_means_r_1m) %>% mutate(cluster1 = getClusterPhenotypes(cluster1), cluster2 = getClusterPhenotypes(cluster2))
sigf_means_r_3m_sigf_pairs <- getSigfPairs(sigf_means_r_3m) %>% mutate(cluster1 = getClusterPhenotypes(cluster1), cluster2 = getClusterPhenotypes(cluster2))

sigf_means_n_0m_sigf_pairs <- getSigfPairs(sigf_means_n_0m) %>% mutate(cluster1 = getClusterPhenotypes(cluster1), cluster2 = getClusterPhenotypes(cluster2))
sigf_means_n_1m_sigf_pairs <- getSigfPairs(sigf_means_n_1m) %>% mutate(cluster1 = getClusterPhenotypes(cluster1), cluster2 = getClusterPhenotypes(cluster2))
sigf_means_n_3m_sigf_pairs <- getSigfPairs(sigf_means_n_3m) %>% mutate(cluster1 = getClusterPhenotypes(cluster1), cluster2 = getClusterPhenotypes(cluster2))


## What are the new interactions?
sigf_means_r_sigf_pairs_new <- getNewInteractions(sigf_means_r_0m_sigf_pairs, sigf_means_r_1m_sigf_pairs)
write.table(sigf_means_r_sigf_pairs_new, "results/cellphonedb/analysis/new_pairs_r.txt", sep = "\t", quote = F, row.names = F)

sigf_means_n_sigf_pairs_new <- getNewInteractions(sigf_means_n_0m_sigf_pairs, sigf_means_n_1m_sigf_pairs)
write.table(sigf_means_n_sigf_pairs_new, "results/cellphonedb/analysis/new_pairs_n.txt", sep = "\t", quote = F, row.names = F)

sigf_means_r_sigf_pairs_new %>% filter(cluster1 == "10 CD8 effector/exhausted") %>% write.table("results/cellphonedb/analysis/new_pairs_r_for_exhausted.txt", sep = "\t", quote = F, row.names = F)

