
getDEGbyCluster <- function(seurat_object, cluster, min_cells = 50){

  message(paste0("===== ", cluster, " ====="))

  ## If under 50 cells to begin with
  if(table(Idents(seurat_object)) %>% as.data.frame() %>% filter(Var1 == cluster) %>% pull(Freq) <= min_cells) return(NULL)

  ## Subet to only cluster
  seurat_cluster         <- subset(seurat_object, ident = cluster)
  Idents(seurat_cluster) <- seurat_cluster$timepoint

  ## Calculate DEG only if there's at least 5 cells per time point
  n_df <- table(Idents(seurat_cluster)) %>% as.data.frame()

  n1 <- n_df %>% filter(Var1 == "baseline") %>% pull(Freq) >= min_cells
  n2 <- n_df %>% filter(Var1 == "6m")       %>% pull(Freq) >= min_cells
  n3 <- n_df %>% filter(Var1 == "12m")      %>% pull(Freq) >= min_cells
  n4 <- n_df %>% filter(Var1 == "relapse")  %>% pull(Freq) >= min_cells

  if(length(n1) == 0) n1 <- FALSE
  if(length(n2) == 0) n2 <- FALSE
  if(length(n3) == 0) n3 <- FALSE
  if(length(n4) == 0) n4 <- FALSE


  cluster_markers_2v1 <- NULL
  cluster_markers_3v1 <- NULL
  cluster_markers_3v2 <- NULL

  cluster_markers_4v1 <- NULL
  cluster_markers_4v2 <- NULL
  cluster_markers_4v3 <- NULL

  if(n1 & n2) cluster_markers_2v1 <- FindMarkers(object = seurat_cluster, ident.1 = "baseline",  ident.2 = "6m",  only.pos = F, min.pct = 0.05, logfc.threshold = 0.1, return.thresh = 1, do.print = T, test.use = "t") %>% add_rownames(var = "gene") %>% mutate(timepoint = "2v1") %>% mutate(p_val_adj = p.adjust(p_val, method = "BH"))
  if(n3 & n1) cluster_markers_3v1 <- FindMarkers(object = seurat_cluster, ident.1 = "baseline",  ident.2 = "12m", only.pos = F, min.pct = 0.05, logfc.threshold = 0.1, return.thresh = 1, do.print = T, test.use = "t") %>% add_rownames(var = "gene") %>% mutate(timepoint = "3v1") %>% mutate(p_val_adj = p.adjust(p_val, method = "BH"))
  if(n3 & n2) cluster_markers_3v2 <- FindMarkers(object = seurat_cluster, ident.1 = "6m", ident.2 = "12m", only.pos = F, min.pct = 0.05, logfc.threshold = 0.1, return.thresh = 1, do.print = T, test.use = "t") %>% add_rownames(var = "gene") %>% mutate(timepoint = "3v2") %>% mutate(p_val_adj = p.adjust(p_val, method = "BH"))

  if(n4 & n1) cluster_markers_4v1 <- FindMarkers(object = seurat_cluster, ident.1 = "baseline", ident.2 = "relapse", only.pos = F, min.pct = 0.05, logfc.threshold = 0.1, return.thresh = 1, do.print = T, test.use = "t") %>% add_rownames(var = "gene") %>% mutate(timepoint = "4v1") %>% mutate(p_val_adj = p.adjust(p_val, method = "BH"))
  if(n4 & n2) cluster_markers_4v2 <- FindMarkers(object = seurat_cluster, ident.1 = "6m", ident.2 = "relapse", only.pos = F, min.pct = 0.05, logfc.threshold = 0.1, return.thresh = 1, do.print = T, test.use = "t") %>% add_rownames(var = "gene") %>% mutate(timepoint = "4v2") %>% mutate(p_val_adj = p.adjust(p_val, method = "BH"))
  if(n4 & n3) cluster_markers_4v3 <- FindMarkers(object = seurat_cluster, ident.1 = "12m", ident.2 = "relapse", only.pos = F, min.pct = 0.05, logfc.threshold = 0.1, return.thresh = 1, do.print = T, test.use = "t") %>% add_rownames(var = "gene") %>% mutate(timepoint = "4v3") %>% mutate(p_val_adj = p.adjust(p_val, method = "BH"))

  df <- rbind(cluster_markers_2v1, cluster_markers_3v1, cluster_markers_3v2, cluster_markers_4v1, cluster_markers_4v2, cluster_markers_4v3)

  if(!is.null(df)) df <- df %>% filter(p_val_adj < 0.05) %>% mutate(cluster = cluster, direction = ifelse(avg_logFC > 0, "up", "down"))
  return(df)

}
