
#### Bulk tcrb
## Prepare files for individual gliph

vdjToGliph <- function(vdj_df){

  ## Write gliph files to the vdj files
  # @ param
  # input: df from vdj

  df <- vdj_df %>% dplyr::select(cdr3aa, v, j, name, freq) %>%
    dplyr::rename("CDR3b"   = "cdr3aa",
                  "TRBV"    = "v",
                  "TRBJ"    = "j",
                  "Patient" =  name,
                  "Counts"  = freq)

  return(df)

}


  files <- list.files(paste0("tcrb_data/unsorted/dasastop/"), full.names = T);
  x <- lapply(files, FUN = function(x){
    message(x);
    x %>% fread %>% filter(count > 1) %>%
      mutate(name = extractFileName(x)) %>%
      vdjToGliph
  }) %>% rbindlist()


## For GLIPH2
df_no_singletons <- x %>% mutate(CDR3a = NA) %>% dplyr::select(CDR3b, TRBV, TRBJ, CDR3a, Patient, Counts) %>% dplyr::rename("subject:condition" = "Patient",	"count" = "Counts")
df_no_singletons <- rbind(df_no_singletons, df)
fwrite(df_no_singletons, "/Users/janihuuh/Dropbox/applications/gliph2/data/cml_total_tcr.txt", sep = "\t", quote = F, row.names = F)

hla_no_singletons <- data.frame(`subject:condition` = unique(df_no_singletons$`subject:condition`)) %>% mutate(allele1 = NA, allele2 = NA)
fwrite(hla_no_singletons, "/Users/janihuuh/Dropbox/applications/gliph2/data/cml_hla.txt", sep = "\t", quote = F, row.names = F)



## Follow the epitope-specific clusters
tcell_clusters <- cml_seurat@meta.data %>% group_by(cluster, tcr = !is.na(new_clonotypes_id)) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% filter(tcr == T) %>% filter(prop > 0.4) %>% pull(cluster) %>% as.character()
cells.to.keep <- cml_seurat@meta.data %>% filter(!is.na(new_clonotypes_id) & cluster %in% tcell_clusters) %>% pull(barcode)
tcr_seurat <- subset(cml_seurat, cells = cells.to.keep)

df <- tcr_seurat@meta.data %>% group_by(trb_cdr3s_aa, v_trb, j_trb, tra_cdr3s_aa, orig.ident) %>% summarise(count = n()) %>% ungroup() %>%
  mutate(orig.ident = paste0(orig.ident, ":", "cml")) %>%
  dplyr::select(trb_cdr3s_aa, v_trb, j_trb, tra_cdr3s_aa,orig.ident,count) %>%
  dplyr::rename("CDR3b" = trb_cdr3s_aa, "TRBV" = v_trb, "TRBJ" = j_trb, "CDR3a" = tra_cdr3s_aa, "subject:condition" = "orig.ident")
df <- df %>% filter(CDR3b != "" & TRBV != "")
df[df == ""] <- NA
fwrite(df, "/Users/janihuuh/Dropbox/applications/gliph2/data/cml_total_tcr.txt", sep = "\t", quote = F, row.names = F, col.names = F)

hla_no_singletons <- data.frame(`subject:condition` = unique(df$`subject:condition`)) %>% mutate(allele1 = NA, allele2 = NA)
fwrite(hla_no_singletons, "/Users/janihuuh/Dropbox/applications/gliph2/data/cml_hla.txt", sep = "\t", quote = F, row.names = F, col.names = F)









## Analyze gliph2
gliph_df <- readGliphFile("results/manuscript/epitope/gliph2_tot.csv") # %>% filter(number_subject > 2 & number_unique_cdr3 >= 5)
gliph_df <- gliph_df %>% filter(vb_score < 0.05 & number_unique_cdr3 >= 3)
colnames(gliph_df) <- make.names(colnames(gliph_df))
gliph_df <- gliph_df %>% dplyr::select(-c(HLA.A:HLA.DRB5))


unique(gliph_df$pattern) %>% length


df <- gliph_df %>% group_by(timepoint, pattern) %>% summarise(n = n())
ggplot(df, aes(timepoint,n,group=pattern, color = pattern)) + geom_point() + geom_path() + scale_color_manual(values = getPalette(28)) + ggpubr::rotate_x_text(angle = 45) + labs(y = "unique TCRs in epitope-specific cluster") +
  ggrepel::geom_text_repel(data = subset(df, n>5), aes(label = pattern), nudge_y = 0.5, nudge_x = 0)
ggsave("results/manuscript/gliph/line_gliph.pdf", width = 7, height = 4)


### co-occurance of clusters
join_df <- gliph_df %>% mutate(key =  paste(TcRa, TcRb, V, J, sep = "_")) %>% select(pattern,key)
gvhd_seurat_tcr$key = paste(gvhd_seurat_tcr$tra_cdr3s_aa, gvhd_seurat_tcr$trb_cdr3s_aa, gvhd_seurat_tcr$v_trb, gvhd_seurat_tcr$j_trb, sep = "_")

addGliph <- function(df, join_df){

  for(epitope in unique(join_df$pattern)){
    df <- df %>% mutate((!!as.name(epitope)) := ifelse(key %in% subset(join_df, pattern == epitope)$key, 1, 0))
  }
  colnames(df) <- make.names(colnames(df))
  return(df)
}

df <- addGliph(df = gvhd_seurat_tcr@meta.data, join_df)
rownames(df) <- df$barcode

epitope_df <- df %>% select(SPRT.G:SL.ET)
epitope_df <- epitope_df[rowSums(epitope_df) != 0, ]
rownames(epitope_df) <- NULL

p <- pheatmap::pheatmap(epitope_df, show_rownames = F, legend = F)
ggsave(plot = print(p), "results/manuscript/gliph/heatmap_cooccurance.pdf", width = 5, height = 6)

data.frame(n = colSums(epitope_df)) %>% add_rownames(var = "pattern") %>%
  ggplot(aes(reorder(pattern,n),n)) + geom_bar(stat = "identity") + ggpubr::rotate_x_text(45) + labs(x = "", y = "nCells")
ggsave("results/manuscript/gliph/bar_tcr.pdf", width = 5, height = 4)
fwrite(df, "results/manuscript/gliph/pattern_occurance.txt", sep = "\t", quote = F, row.names = F)


## Plot gliph results
cols.to.cluster <- data.frame(n = colSums(epitope_df)) %>% add_rownames(var = "pattern") %>% arrange(desc(n)) %>% pull(pattern)
colors <- getPalette(length(cols.to.cluster))

p <- NULL; i <- 1
for(clustering_column in cols.to.cluster){

  message(clustering_column)
  cells   <- gvhd_seurat_tcr@meta.data %>% filter(!!as.name(clustering_column) == 1) %>% pull(barcode)


  p[[i]] <- DimPlot(gvhd_seurat_tcr, reduction = "umap", cells.highlight = cells, cols.highlight = c(colors[i], "gray20"), label = F) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none") +
    labs(title = clustering_column) + labs(x = "UMAP 1", y = "UMAP 2")
  i <- i + 1

}

png("results/manuscript/gliph/umap_gliph.png", width = 1024*c(3/5), height = 1024)
do.call(grid.arrange, c(p, ncol = 4))
dev.off()





DimPlot(gvhd_seurat_tcr, reduction = "umap", split.by = "S.DFQET", cols = getPalette3(12), label = T, repel = T) + theme_bw(base_size = 12) + theme(legend.position = "none") + labs(x = "UMAP1", y = "UMAP2")
ggsave("results/manuscript/gliph/latent_umap.png", width = 7, height = 6)

gvhd_seurat_tcr@meta.data %>% filter(S.DFQET == 1) %>% group_by(timepoint, cluster_new) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>%
  ggplot(aes(timepoint, prop, group = cluster_new, color = cluster_new)) + geom_path() + geom_point()
