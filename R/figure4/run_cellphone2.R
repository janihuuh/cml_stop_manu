


## Visualize; dg-phase
p <- cast(sigf_means_dg_hsc_pairs, cluster1 ~ cluster2, sum)
rownames(p) <- colnames(p)[-1]
p <- p[,-1]
colnames(p) <- substr(colnames(p), 3, nchar(colnames(p)))
rownames(p) <- substr(rownames(p), 3, nchar(rownames(p)))
rownames(p) <- gsub("\\_", "\\/", rownames(p))
colnames(p) <- gsub("\\_", "\\/", colnames(p))

dg_df <- p[!rownames(p) %in% getClusterPhenotypesLSC(as.character(0:7)), colnames(p) %in% getClusterPhenotypesLSC(as.character(0:7))]



## Visualize; tki_hsc
p <- cast(sigf_means_tki_hsc_pairs, cluster1 ~ cluster2, sum)
rownames(p) <- colnames(p)[-1]
p <- p[,-1]
colnames(p) <- substr(colnames(p), 3, nchar(colnames(p)))
rownames(p) <- substr(rownames(p), 3, nchar(rownames(p)))
rownames(p) <- gsub("\\_", "\\/", rownames(p))
colnames(p) <- gsub("\\_", "\\/", colnames(p))

tki_df <- p[!rownames(p) %in% getClusterPhenotypesLSC(as.character(0:7)), colnames(p) %in% getClusterPhenotypesLSC(as.character(0:7))]


## Visualize; tfr
p <- cast(sigf_means_tfr_hsc_pairs, cluster1 ~ cluster2, sum)
rownames(p) <- colnames(p)[-1]
p <- p[,-1]
colnames(p) <- substr(colnames(p), 3, nchar(colnames(p)))
rownames(p) <- substr(rownames(p), 3, nchar(rownames(p)))
rownames(p) <- gsub("\\_", "\\/", rownames(p))
colnames(p) <- gsub("\\_", "\\/", colnames(p))

tfr_df <- p[!rownames(p) %in% getClusterPhenotypesLSC(as.character(0:7)), colnames(p) %in% getClusterPhenotypesLSC(as.character(0:7))]



## Visualize; rp
p <- cast(sigf_means_rp_hsc_pairs, cluster1 ~ cluster2, sum)
rownames(p) <- colnames(p)[-1]
p <- p[,-1]
colnames(p) <- substr(colnames(p), 3, nchar(colnames(p)))
rownames(p) <- substr(rownames(p), 3, nchar(rownames(p)))
rownames(p) <- gsub("\\_", "\\/", rownames(p))
colnames(p) <- gsub("\\_", "\\/", colnames(p))

rp_df <- p[!rownames(p) %in% getClusterPhenotypesLSC(as.character(0:7)), colnames(p) %in% getClusterPhenotypesLSC(as.character(0:7))]





#####

dg_df
tki_df
tfr_df
rp_df

df <- data.frame(name = rownames(dg_df), dg = rowSums(dg_df)) %>% left_join(data.frame(name = rownames(tki_df), tki = rowSums(tki_df))) %>% left_join(data.frame(name = rownames(tfr_df), tfr = rowSums(tfr_df))) %>% left_join(data.frame(name = rownames(rp_df), rp = rowSums(rp_df))) %>% melt(id = "name") 

ggplot(df, aes(variable,value,group=name,color=name)) + geom_path() + ggrepel::geom_text_repel(data = subset(df, variable == "rp"), aes(label = name)) + theme(legend.position = "none")




### only bcr-abl pos cells
dg_bcrabl_df  <- dg_df[,which(grepl("BCRABL1\\+", colnames(dg_df)))]
tki_bcrabl_df <- tki_df[,which(grepl("BCRABL1\\+", colnames(tki_df)))]
tfr_bcrabl_df <- tfr_df[,which(grepl("BCRABL1\\+", colnames(tfr_df)))]
rp_bcrabl_df  <- rp_df[,which(grepl("BCRABL1\\+", colnames(rp_df)))]

df <- data.frame(name = rownames(dg_bcrabl_df), dg = rowSums(dg_bcrabl_df)) %>% 
  left_join(data.frame(name = rownames(tki_bcrabl_df), tki = rowSums(tki_bcrabl_df))) %>% 
  left_join(data.frame(name = rownames(tfr_bcrabl_df), tfr = rowSums(tfr_bcrabl_df))) %>% 
  left_join(data.frame(name = rownames(rp_bcrabl_df), rp = rowSums(rp_bcrabl_df))) %>% melt(id = "name") 

ggplot(df, aes(variable,value)) + geom_boxplot(outlier.shape = NA) + geom_path(aes(group=name,color=name), alpha = 0.5) + ggrepel::geom_text_repel(data = subset(df, variable == "rp"), aes(label = name)) + theme(legend.position = "none") + ggpubr::stat_compare_means()

df %>% # left_join(df %>% group_by(name) %>% summarise(total = sum(value))) %>% 
  ggplot(aes(variable,value)) + geom_path(aes(group=name,color=name), alpha = 0.5) + theme(legend.position = "none") + facet_wrap(~name, scales = "free_y") + facets_nice + geom_point(color="gray20") + labs(x = "", y = "nInteractions")
ggsave("results/cml_stop/cellphonedb/path_bcrabl_interactions.pdf", width = 11, height = 8)



