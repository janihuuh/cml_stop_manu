
library(ComplexHeatmap)
library(circlize)

## Define the genes
naive_markers         <- c("TCF7", "SELL", "LEF1", "CCR7")
cytotoxic_markers     <- c("IFNG", "GZMA", "GZMB", "GZMM", "GZMH", "PRF1", "GNLY")
costimulatory_markers <- c("CD28", "TNFRSF14", "ICOS", "TNFRSF9")
inhibitory_markers    <- c("PDCD1", "CTLA4", "LAG3", "HAVCR2")
cd8_em                <- c("GZMK", "CXCR4", "CXCR3", "CD44")
proliferation         <- c("MKI67")
tumor_reactivity      <- c("ENTPD1", "ITGAE")

marker_list           <- list(cytotoxic_markers, costimulatory_markers, inhibitory_markers, proliferation, tumor_reactivity)

marker_genes          <- do.call("c", marker_list)
marker_genes          <- marker_genes[marker_genes %in% rownames(li_seurat)]

marker_amount         <- do.call("c", lapply(marker_list, length))
# split_known           <- rep(c("6 Naive", "5 Cytotoxic", "4 Costimulatory", "3 Inhibitory", "2 EM", "2 proliferation", "1 tumor_reactivity", "1 activity"), marker_amount)
split_known           <- rep(c("5 Cytotoxic", "2 Costimulatory", "4 Inhibitory", "1 Prolif", "3 Tumor \nreactivity"), marker_amount)
# split_known           <- rep(c("Cytotoxic", "Costimulatory", "Inhibitory", "Proliferation", "Tumor reactivity", "Activity"), marker_amount)

# Count mean by cluster
li_exhausted <- pr1_seurat2
hcc_df <- li_exhausted@assays$RNA@data[marker_genes, ] %>% as.matrix()
hcc_df <- rbind("Orig.ident" = as.factor(li_exhausted@meta.data$orig.ident), hcc_df)
colnames(hcc_df) <- NULL


count_mean_exp_known <- function(df, i){

  cluster_mean <- rowMeans(df[ ,which(df[1,] == i)])
  return(cluster_mean)

}

cluster_means_known <- lapply(1:5, count_mean_exp_known, df = hcc_df)
cluster_means_known <- do.call(rbind, cluster_means_known)

## Z-score
cluster_means_known   <- cluster_means_known[,-1]
cluster_means_known_z <- t(scale(cluster_means_known))


## Annotation
epitopes <- li_exhausted$orig.ident %>% unique() %>% as.character()
ha       <- HeatmapAnnotation(names = anno_text(epitopes, rot = 90))

pdf("results/tcr/heatmap_pr1_timepoint_new.pdf", width = 3.5, height = 10)

Heatmap(cluster_means_known_z,
        bottom_annotation = ha,
        #  bottom_annotation = ha2,
        col = colorRamp2(c(-4, 0, 4), c("dodgerblue", "white", "salmon")),
        cluster_columns = F,
        cluster_rows = F,
        row_names_gp = gpar(fontface = 3),
        row_names_side = c("right"),
        show_heatmap_legend = F,
        split = split_known,
        gap = unit(5, "mm"),
        # bottom_annotation_height = unit(2.5, "cm"),
        #  bottom_annotation_height =  unit(6, "cm"),
        name = "Z-score \n(of scaled log-counts)")

dev.off()




pdf("results/tcr/heatmap_pr1_timepoint_scale.pdf", width = 3.5, height = 10)

Heatmap(cluster_means_known_z,
        bottom_annotation = ha,
        #  bottom_annotation = ha2,
        col = colorRamp2(c(-4, 0, 4), c("dodgerblue", "white", "salmon")),
        cluster_columns = F,
        cluster_rows = F,
        row_names_gp = gpar(fontface = 3),
        row_names_side = c("right"),
        show_heatmap_legend = T,
        split = split_known,
        gap = unit(5, "mm"),
        # bottom_annotation_height = unit(2.5, "cm"),
        #  bottom_annotation_height =  unit(6, "cm"),
        name = "Z-score \n(of scaled log-counts)")

dev.off()


