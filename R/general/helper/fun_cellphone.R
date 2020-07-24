

initCellphonedb <- function(df, seurat_object, name = "", sample_n = 50, folder = "results/cellphonedb/input_files/", seed = 123){

  set.seed(seed)

  ## If one wants to use subsampling; remove clusters with less than sample_n cells
  if(!is.null(sample_n)){

    clusters.to.rm  <- df %>% group_by(cluster) %>% summarise(n = n()) %>% filter(n < sample_n) %>% pull(cluster)
    message(paste("remove", clusters.to.rm))
    clusters.to.use <- df %>% group_by(cluster) %>% summarise(n = n()) %>% filter(n >= sample_n) %>% pull(cluster)
    df              <- df %>% filter(cluster %in% clusters.to.use) %>% group_by(cluster)
    df              <- df %>% sample_n(size = sample_n) %>% ungroup()
  }

  ## Remove empty clusters
  df <- df %>% filter(cluster != "")

  counts <- seurat_object@assays$RNA@counts[,df$barcode]
  counts[Matrix::rowSums(counts) != 0, ] %>% as.data.frame %>% add_rownames(var = "Gene") %>%
    dplyr::select(Gene, everything()) %>%
    fwrite(paste0(folder, name, "_counts.txt"), sep = "\t", quote = F, row.names = F)

  df %>% dplyr::select(barcode, cluster) %>% dplyr::rename(Cell = barcode, cell_type = cluster) %>%
    fwrite(paste0(folder, name, "_meta.txt"), sep = "\t", quote = F, row.names = F)

  return(NULL)

}

getPairAmounts <- function(df){

  ## Detect columns with pairs
  cols.to.use    <- grep("\\|", colnames(df), value = T)
  # cols.to.use    <- c("interacting_pair", grep("\\|", colnames(df), value = T))
  df             <- df %>% dplyr::select(cols.to.use)
  pairs.complete <- lapply(cols.to.use, function(x) strsplit(x, "[|]")[[1]][1]) %>% do.call(what = "c")
  pairs.unique   <- unique(pairs.complete)
  # pairs.unique   <- unique(pairs.complete)[-1]

  i <- 1
  df_list <- NULL

  for(pair in pairs.unique){

    # message(pair)
    # pair <- c("interacting_pair", pair)
    df_pair <- df %>% dplyr::select(colnames(df)[pairs.complete %in% pair])

    df_list[[i]] <- apply(df_pair[,], 2, function(x) table(!is.na(x))[2]) %>% # bind_cols(data.frame(interacting_pair = df_pair$interacting_pair)) %>%
      as.data.frame() %>% add_rownames(var = "pair") %>% dplyr::rename(value = ".") %>% mutate(cluster1 = lapply(pair, function(x) strsplit(x, "[|]")[[1]][1]) %>% do.call(what = "c"),
                                                                                               cluster2 = lapply(pair, function(x) strsplit(x, "[|]")[[1]][2]) %>% do.call(what = "c"))

    i <- i + 1

  }

  df1 <- df_list %>% rbindlist()


}

filterRedundantPairs <- function(df1){

  ## Make all pairs
  expand.grid.unique <- function(x, y, incl.eq = TRUE){
    g <- function(i){
      z <- setdiff(y, x[seq_len(i - incl.eq)])
      if(length(z)) cbind(x[i], z, deparse.level = 0)
    }
    do.call(rbind, lapply(seq_along(x), g))
  }


  pairs <- expand.grid.unique(unique(df1$cluster1), unique(df1$cluster1)) %>% as.data.frame() %>% mutate(pairs = paste0(V1, "|", V2))
  df1 <- df1 %>% filter(pair %in% pairs$pairs)
  return(df1)

}

getSigfPairs <- function(df){

  # df = sigf_means_r_pre

  ## Detect columns with pairs
  # cols.to.use    <- grep("\\|", colnames(df), value = T)
  cols.to.use    <- c("interacting_pair", grep("\\|", colnames(df), value = T))
  df             <- df %>% dplyr::select(cols.to.use)
  pairs.complete <- lapply(cols.to.use, function(x) strsplit(x, "[|]")[[1]][1]) %>% do.call(what = "c")
  # pairs.unique   <- unique(pairs.complete)
  pairs.unique   <- unique(pairs.complete)[-1]

  i <- 1
  df_list <- NULL

  for(pair in pairs.unique){

    message(pair)

    pair <- c("interacting_pair", pair)
    df_pair <- df %>% dplyr::select(colnames(df)[pairs.complete %in% pair])

    df_list[[i]]  <- melt(df_pair, id = "interacting_pair") %>% filter(!is.na(value)) %>% mutate(cluster1 = lapply(as.character(variable), function(x) strsplit(x, "[|]")[[1]][1]) %>% do.call(what = "c"),
                                                                                                 cluster2 = lapply(as.character(variable), function(x) strsplit(x, "[|]")[[1]][2]) %>% do.call(what = "c"))

    i <- i + 1

  }

  df1 <- df_list %>% rbindlist()

  ## Make all pairs
  expand.grid.unique <- function(x, y, incl.eq = TRUE){
    g <- function(i){
      z <- setdiff(y, x[seq_len(i - incl.eq)])
      if(length(z)) cbind(x[i], z, deparse.level = 0)
    }
    do.call(rbind, lapply(seq_along(x), g))
  }


  pairs <- expand.grid.unique(unique(df1$cluster1), unique(df1$cluster1)) %>% as.data.frame() %>% mutate(pairs = paste0(V1, "|", V2))
  df1 <- df1 %>% filter(variable %in% pairs$pairs)
  return(df1)

}

getNewInteractions <- function(df1,df2){

  df1 <- df1 %>% mutate(temp = paste(interacting_pair, variable))
  df2 <- df2 %>% mutate(temp = paste(interacting_pair, variable))

  new_df <- df2[!df2$temp %in% df1$temp, ]

}

plotPairs <- function(df1, df2){

  # sigf_r <- sigf_r %>% dplyr::rename(group1 = median_pre, group2 = median_post, p = p.value) %>% mutate(y.position = 30)
  # sigf_r <- sigf_r %>% select(group1, group2, p, y.position)

  df <- merge(df1,df2,by="pair") %>% dplyr::rename(pre = value.x, post = value.y)

  nClustersTemp <- df$cluster1.x %>% unique() %>% length()
  df <- df %>% dplyr::select(pair, pre, post, cluster1.x, cluster2.x) %>% melt(id = c("pair", "cluster1.x", "cluster2.x"))
  df <- df %>% mutate(phenotype = extractCoarsePhenotype(cluster2.x))

  ggplot(df, aes(variable,value,group=pair,fill=phenotype)) +
    geom_point(shape = 21, size = 3) + geom_path() +
    # ggpubr::stat_pvalue_manual(sigf_r) +
    # ggsignif::geom_signif(comparisons = list(c("pre", "post")), test.args = c(paired = T)) +
    # ggrepel::geom_label_repel(data = subset(df, variable == "post"), aes(variable,value,label=cluster2.x,color=cluster2.x), fill = NA) +
    facet_wrap(~cluster1.x, scales = "free") + facets_nice + labs(fill = "", x = "", y = "# of sigf receptor-ligand interactions") +
    # scale_fill_manual(values = getPalette(nClustersTemp))
    scale_fill_manual(values = brewer.pal(9, "Pastel1"))

}

plotNewPairs <- function(df1, df2){

  df <- merge(df1,df2,by="pair") %>% dplyr::rename(pre = value.x, post = value.y)

  nClustersTemp <- df$cluster1.x %>% unique() %>% length()
  df <- df %>% dplyr::select(pair, pre, post, cluster1.x, cluster2.x)

  df$post <- c(df$post - df$pre) %>% as.numeric
  df$pre  <- 0

  nClusters <- df$cluster2.x %>% unique %>% length


  df <- df %>% melt(id = c("pair", "cluster1.x", "cluster2.x")) %>% mutate(phenotype = extractCoarsePhenotype(cluster2.x))

  df %>% filter(variable == "post") %>%
    ggplot(aes(cluster2.x, value, fill = cluster2.x)) +
    geom_bar(stat = "identity") + facet_wrap(~cluster1.x) + scale_fill_manual(values = getPalette(nClusters)) +
    theme(axis.text.x = element_blank())
  # theme(axis.text.x = element_text(angle = 90, hjust = 1))

}

getNewPairs <- function(df1, df2){

  df <- merge(df1,df2,by="pair") %>% dplyr::rename(pre = value.x, post = value.y)

  nClustersTemp <- df$cluster1.x %>% unique() %>% length()
  df <- df %>% dplyr::select(pair, pre, post, cluster1.x, cluster2.x)

  df$post <- c(df$post - df$pre) %>% as.numeric
  df$pre  <- 0

  nClusters <- df$cluster2.x %>% unique %>% length

  df <- df %>% melt(id = c("pair", "cluster1.x", "cluster2.x")) %>% mutate(phenotype = extractCoarsePhenotype(cluster2.x))


}

testSigfPairs <- function(df1, df2, paired = T){

  df <- merge(df1,df2,by="pair") %>% dplyr::rename(pre = value.x, post = value.y)

  df_list <- NULL
  i <- 1

  for(cluster in unique(df$cluster1.x)){

    # cluster = unique(df$cluster1.x)[1]
    df_temp <- df %>% filter(cluster1.x == cluster)

    df_list[[i]] <- wilcox.test(df_temp$pre, df_temp$post, paired = paired) %>% broom::tidy() %>%
      mutate(cluster = cluster, median_pre = median(df_temp$pre, na.rm = T), median_post = median(df_temp$post, na.rm = T),
             direction = ifelse(median_post > median_pre, "up", "down"))
    i <- i + 1

  }

  df_list %>% rbindlist()

}


expand.grid.unique <- function(x, y, incl.eq = TRUE){
  g <- function(i){
    z <- setdiff(y, x[seq_len(i - incl.eq)])
    if(length(z)) cbind(x[i], z, deparse.level = 0)
  }
  do.call(rbind, lapply(seq_along(x), g))
}

f <- function(m) {
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  m
}


plotDiffHeatmap <- function(df1, df2){

  require(igraph)
  require(ComplexHeatmap)

  ## Make all pairs
  expand.grid.unique <- function(x, y, incl.eq = TRUE){
    g <- function(i){
      z <- setdiff(y, x[seq_len(i - incl.eq)])
      if(length(z)) cbind(x[i], z, deparse.level = 0)
    }
    do.call(rbind, lapply(seq_along(x), g))
  }

  f <- function(m) {
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    m
  }


  pairs <- expand.grid.unique(unique(df1$cluster1), unique(df1$cluster1)) %>% as.data.frame() %>% mutate(pairs = paste0(V1, "|", V2))

  df1 <- df1 %>% filter(pair %in% pairs$pairs)
  df2 <- df2 %>% filter(pair %in% pairs$pairs)



  # df2 <- df %>%
  #   mutate(cluster1 = gsub("\\?", "na", cluster1), cluster2 = gsub("\\?", "na", cluster2)) %>%
  #   mutate(cluster1 = gsub("\\ ", "_", cluster1), cluster2 = gsub("\\ ", "_", cluster2)) %>%
  #   mutate(cluster1 = substr(cluster1, 4,nchar(cluster1)), cluster2 = substr(cluster2, 4,nchar(cluster2))) %>% dplyr::select(cluster1, cluster2, value)

  df_mat1 <- cast(df1, cluster1 ~ cluster2, sum)
  rownames(df_mat1) <- df_mat1$cluster1

  df_mat2 <- cast(df2, cluster1 ~ cluster2, sum)
  rownames(df_mat2) <- df_mat2$cluster1

  to.use1 <- rownames(df_mat1) %in% rownames(df_mat2)
  to.use2 <- rownames(df_mat2) %in% rownames(df_mat1)

  df_mat2 <- df_mat2 %>% as.data.frame()
  df_mat2[to.use2,c(F,to.use2)]


  df_diff <- log2(df_mat2[to.use2,c(F,to.use2)] / df_mat1[to.use1,c(F,to.use1)])
  # df_diff <- do.call(data.frame,lapply(df_diff, function(x) replace(x, is.infinite(x),NA)))


  df_diff <- apply(df_diff, 1, function(x){ x[!is.finite(x)] <- NA; return(x)} )


  anno_df <- data.frame(phenotype = extractCoarsePhenotype(rownames(df_mat1)))
  rownames(anno_df) <- rownames(df_mat1)

  nAnno <- anno_df$phenotype %>% unique %>% length()
  # anno_col <- data.frame(phenotype = unique(anno_df$phenotype),
  #                        col       = getPalette3(nAnno))
  # anno_col <- split(anno_col, seq(nrow(anno_col)))

  col = getPalette(nAnno)
  col = brewer.pal(nAnno, "Pastel1")

  # anno_col <- list(phenotype = c("NK" = col[1], "CD4+" = col[2], "CD8+" = col[3], "Monocytes" = col[4], "B-cell" = col[5]), "?" = col[6]))

  # pheatmap::pheatmap(df_diff, cluster_rows = T, cluster_cols = T, color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(20),
  #                    annotation_row = anno_df, annotation_col = anno_df, annotation_colors = anno_col)

  pheatmap::pheatmap(f(t(df_diff)),
                     cluster_rows = F, cluster_cols = F,
                     color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(20),
                     annotation_row = anno_df,
                     annotation_col = anno_df,

                     # annotation_colors = anno_col,
                     na_col = "lightgrey")









  # pheatmap::pheatmap(df_mat[,-1], cluster_rows = T, cluster_cols = T, color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(10))

  # write.table(df_mat, "results/cellphonedb/out/sigf_means_r_0m_mat.txt", sep = "\t", quote = F, row.names = T)
  # write.table(  as.matrix(df[,3:4]), "results/cellphonedb/out/sigf_means_r_0m_el.txt", sep = "\t", quote = F, row.names = T)
  # g <- graph_from_data_frame(df_mat, directed = FALSE, vertices = NULL)
  # g <- graph_from_adjacency_matrix(df_mat)
  # g <- graph_from_edgelist(as.matrix(df[,3:4]), directed = FALSE)


}




plotDiffHeatmap <- function(df1, df2){

  require(igraph)
  require(ComplexHeatmap)

  ## Make all pairs
  expand.grid.unique <- function(x, y, incl.eq = TRUE){
    g <- function(i){
      z <- setdiff(y, x[seq_len(i - incl.eq)])
      if(length(z)) cbind(x[i], z, deparse.level = 0)
    }
    do.call(rbind, lapply(seq_along(x), g))
  }

  f <- function(m) {
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    m
  }


  pairs <- expand.grid.unique(unique(df1$cluster1), unique(df1$cluster1)) %>% as.data.frame() %>% mutate(pairs = paste0(V1, "|", V2))

  # df1 <- df1 %>% filter(pair %in% pairs$pairs)
  # df2 <- df2 %>% filter(pair %in% pairs$pairs)

  df_mat1 <- cast(df1, cluster1 ~ cluster2, sum)
  rownames(df_mat1) <- df_mat1$cluster1

  df_mat2 <- cast(df2, cluster1 ~ cluster2, sum)
  rownames(df_mat2) <- df_mat2$cluster1

  to.use1 <- rownames(df_mat1) %in% rownames(df_mat2)
  to.use2 <- rownames(df_mat2) %in% rownames(df_mat1)

  df_mat2 <- df_mat2 %>% as.data.frame()
  df_mat2[to.use2,c(F,to.use2)]

  df_diff <- log2(df_mat2[to.use2,c(F,to.use2)] / df_mat1[to.use1,c(F,to.use1)])
  df_diff <- apply(df_diff, 1, function(x){ x[is.nan(x)] <- 0; return(x)} )
  df_diff <- apply(df_diff, 1, function(x){ x[!is.finite(x)] <- NA; return(x)} )

  # df_diff[is.na(df_diff)] <- 2

  anno_df <- data.frame(phenotype = getClusterCoarsePhenotype(rownames(df_mat1)))
  rownames(anno_df) <- rownames(df_mat1)

  nAnno <- anno_df$phenotype %>% unique %>% length()
  col = getPalette(nAnno)
  col = brewer.pal(nAnno, "Pastel1")

  anno_col <- list(phenotype = c("NK" = col[1], "CD4" = col[2], "CD8" = col[3], "CD4CD8" = col[4], "B" = col[5], "other" = col[6]))

  # ComplexHeatmap::Heatmap(f(t(df_diff)))

  pheatmap::pheatmap(f(t(df_diff)),
                     cluster_rows = T,
                     cluster_cols = T,
                     color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(seq(-2, 2, 0.1))),
                     # annotation_row = anno_df,
                     # annotation_col = anno_df,
                     # annotation_colors = anno_col,
                     breaks = seq(-2, 2, 0.1),
                     na_col = "lightgrey")

}
