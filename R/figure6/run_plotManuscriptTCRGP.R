
## TCRGP
tcrgp_results <- fread("results/tcrgp/tcrgp_predictions.txt")
colnames(tcrgp_results) <- make.unique(colnames(tcrgp_results))
tcrgp_results <- tcrgp_results %>% dplyr::select(-c(v:clonotype_name))

pr1_tcrb <- fread("tcrb_data/pr1_total_filtered.txt")

meta <- cml_seurat@meta.data %>% left_join(tcrgp_results, by = "barcode")
rownames(meta) <- meta$barcode
cml_seurat2 <- cml_seurat
cml_seurat2@meta.data <- meta
cml_seurat2$pr1_specific <- ifelse(cml_seurat2$trb_cdr3s_aa %in% pr1_tcrb$cdr3aa, "pr1_specific", "other")

cells.to.keep <- cml_seurat2@meta.data %>%
  filter(target.y == "anti-viral" &
           patient %in% c("Batch 5", "Batch 7") &
           pred_epitope.y != "InfA_HA_PKY" &
           cluster %in% c("7 T cells activated", "5 CD8 effectory/exhausted", "3 CD8 effectory", "14 CD8 EM", "11 CD8 CM/naive") |
           pr1_specific == "pr1_specific" & patient %in% c("Batch 5", "Batch 7") & cluster %in% c("7 T cells activated", "5 CD8 effectory/exhausted", "3 CD8 effectory", "14 CD8 EM", "11 CD8 CM/naive")) %>% pull(barcode)
viral_seurat  <- subset(cml_seurat2, cells = cells.to.keep) %>% getLatentUMAP() %>% fixSeurat() %>% getLatentClustering()
viral_seurat$pred_epitope.y[viral_seurat$pr1_specific == "pr1_specific"] <- "PR1"

res                <- c(seq(0.1, 1, 0.1), seq(1.2, 2, 0.2), 2.5, 3)
clustering_columns <- grep("res", colnames(viral_seurat@meta.data), value = T)
clustering_columns <- clustering_columns[order(substr(clustering_columns, 13, nchar(clustering_columns)) %>% as.numeric())]

q <- NULL; i <- 1
for(clustering_column in clustering_columns){
  q[[i]] <- viral_seurat@meta.data[,clustering_column] %>% levels %>% length
  i <- i + 1
}


data.frame(resolution = res, nClusters = q) %>%
  ggplot(aes((resolution),nClusters), label=nClusters) + geom_point(shape = 21) + theme_bw()
# ggsave("results/public/scatter_res_nClusters.png", width = 4, height = 3)

DimPlot(viral_seurat, group.by = "RNA_snn_res.0.4", label = T, repel = T, ncol = 3) + theme_bw(base_size = 12) + labs(x = "UMAP1", y = "UMAP2") + theme(legend.position = "none") + scale_color_manual(values = getPalette(8))

DimPlot(viral_seurat, split.by = "patient", label = F, repel = T, ncol = 3) + theme_bw(base_size = 12) + labs(x = "UMAP1", y = "UMAP2") + theme(legend.position = "none")
DimPlot(viral_seurat, split.by = "orig.ident", label = F, repel = T, ncol = 3) + theme_bw(base_size = 12) + labs(x = "UMAP1", y = "UMAP2") + theme(legend.position = "none")

viral_deg <- FindAllMarkers(viral_seurat, test.use = "t") %>% filter(p_val_adj < 0.05) %>% mutate(dir = ifelse(avg_logFC > 0, "up", "down")) %>% filter(dir == "up")
clonality_genes <- getClonalityGenes(viral_seurat)
fwrite(viral_deg, "results/manuscript/epitope/viral_deg.txt", sep = "\t", quote = F, row.names = F)

viral_deg %>% filter(!gene %in% c(clonality_genes, unwanted_genes, grep("^MT\\-", rownames(viral_seurat), value = T))) %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC) %>% View
viral_deg %>% filter(!gene %in% c(clonality_genes, unwanted_genes, grep("^MT\\-", rownames(viral_seurat), value = T))) %>% filter(gene %in% inhibitory_long)



Idents(viral_seurat) <- viral_seurat$RNA_snn_res.0.4 %>% getClusterPhenotypesViral()
viral_seurat$clusters <- viral_seurat$RNA_snn_res.0.4 %>% getClusterPhenotypesViral()

DotPlot(viral_seurat, features = rev(big_markers), cols = "RdYlBu") + ggpubr::rotate_x_text(90) + labs(x = "", y = "") + theme(legend.position = "top")
ggsave("results/manuscript/epitope/dot_viral_big.pdf", width = 14, height = 5)

DotPlot(viral_seurat, features = rev(unique(do.call(guo_markers, what = "c"))), cols = "RdYlBu") + ggpubr::rotate_x_text(90) + labs(x = "", y = "") + theme(legend.position = "top")
ggsave("results/manuscript/epitope/dot_viral_guo.pdf", width = 14, height = 5)

zhang_cd8_genes <- do.call(zhang_cd8_markers, what = "c") %>% unique()
DotPlot(viral_seurat, features = rev(zhang_cd8_genes), cols = "RdYlBu") + ggpubr::rotate_x_text(90) + labs(x = "", y = "") + theme(legend.position = "top")
ggsave("results/manuscript/epitope/dot_viral_zhang_cd8.pdf", width = 14, height = 5)

DimPlot(viral_seurat, label = T, repel = T, ncol = 3) + theme_bw(base_size = 12) + labs(x = "UMAP1", y = "UMAP2") + theme(legend.position = "none") + scale_color_manual(values = getPalette(8))
ggsave("results/manuscript/epitope/umap.png", width = 5, height = 4)

DimPlot(viral_seurat, label = T, repel = T, ncol = 3, group.by = "cluster") + theme_bw(base_size = 12) + labs(x = "UMAP1", y = "UMAP2") + theme(legend.position = "none") + scale_color_manual(values = getPalette(20))
ggsave("results/manuscript/epitope/umap_cluster.png", width = 5, height = 4)

DimPlot(viral_seurat, split.by = "pred_epitope.y", label = F, repel = T, ncol = 3) + theme_bw(base_size = 12) + labs(x = "UMAP1", y = "UMAP2") + theme(legend.position = "none") + scale_color_manual(values = getPalette(8))
ggsave("results/manuscript/epitope/umap_epitope.png", width = 7, height = 6)

DimPlot(viral_seurat, split.by = "patient", label = F, repel = T, ncol = 3) + theme_bw(base_size = 12) + labs(x = "UMAP1", y = "UMAP2") + theme(legend.position = "none")
ggsave("results/manuscript/epitope/umap_patient.png", width = 7, height = 4)


viral_seurat@meta.data %>% filter(patient == "Batch 5") %>% group_by(pred_epitope.y, clusters) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = "",prop,fill=clusters)) + geom_bar(stat = "identity") + coord_polar("y", start=0) + facet_wrap(~pred_epitope.y) + facets_nice + labs(x = "", y = "") + scale_fill_manual(values = getPalette(8))
ggsave("results/manuscript/epitope/pie_batch5.pdf", width = 7, height = 6)

viral_seurat@meta.data %>% filter(patient == "Batch 7") %>% group_by(pred_epitope.y, clusters) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = "",prop,fill=clusters)) + geom_bar(stat = "identity") + coord_polar("y", start=0) + facet_wrap(~pred_epitope.y) + facets_nice  + labs(x = "", y = "") + scale_fill_manual(values = getPalette(8))
ggsave("results/manuscript/epitope/pie_batch7.pdf", width = 7, height = 6)

saveRDS(viral_seurat, "results/viral_seurat.rds")


getClusterPhenotypesViral <- function(clusters){

  clusters <- plyr::revalue(clusters, replace   = c("0"  = "0 CD8 effector / exhausted",
                                                    "1"  = "1 CD8 effector TIM3+",
                                                    "2"  = "2 CD8 CM",
                                                    "3"  = "3 CD8 EM",
                                                    "4"  = "4 CD8 IFNg"))

  return(clusters)

}

