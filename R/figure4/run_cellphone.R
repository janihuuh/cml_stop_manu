
getClusterPhenotypesLSC <- function(clusters){
  
  clusters <- plyr::revalue(clusters, replace = c(
    
    "0"  = "0 MEP CD38+",
    "1"  = "1 MEP CD38+ cycling" ,
    "2"  = "2 BCRABL1- CD38- HSC" ,
    "3"  = "3 BCRABL1+ GMP cycling" ,
    "4"  = "4 BCRABL1+ HSC/CMP cycling" ,
    "5"  = "5 BCRABL1- CD38- HSC" ,
    "6"  = "6 BCRABL1+ GMP cycling" ,
    "7"  = "7 MEP low quality" ))
  
  return(clusters)
  
}   

diet_seurat <- readRDS("results/cml_stop/diet_seurat.rds")

## Run cellphone db. Make subselection of 1000 cells from responders at time 0 and time 1
dir.create("results/cml_stop/cellphonedb/", showWarnings = F)
dir.create("results/cml_stop/cellphonedb/input_files", showWarnings = F)
dir.create("results/cml_stop/cellphonedb/analysis", showWarnings = F)

sample_n = 10

################################ Init cellphone input files

Idents(diet_seurat) <- diet_seurat$cluster
diet_seurat$project[is.na(diet_seurat$project)] <- "hsc"
diet_seurat$cluster <- gsub("\\/", "_", diet_seurat$cluster)
diet_seurat$cluster <- gsub("\\\n", "_", diet_seurat$cluster)
diet_seurat$barcode <- colnames(diet_seurat)

df <- data.frame(diet_seurat@meta.data, cluster2 = diet_seurat$cluster) %>% filter(project %in% c("CML dg", "hsc")) 
df1 <- fread("results/cml_stop/cellphonedb/input_files/dg_hsc_counts.txt")
df2 <- fread("results/cml_stop/cellphonedb/input_files/dg_hsc_meta.txt")

data.frame(diet_seurat@meta.data, cluster2 = diet_seurat$cluster) %>% filter(project %in% c("CML dg", "hsc")) %>% initCellphonedb(seurat_object = diet_seurat, name = "dg_hsc", sample_n = sample_n, folder = "results/cml_stop/cellphonedb/input_files/")
data.frame(diet_seurat@meta.data, cluster2 = diet_seurat$cluster) %>% filter(project %in% c("hsc") | timepoint == "baseline") %>% initCellphonedb(seurat_object = diet_seurat, name = "tki_hsc", sample_n = sample_n, folder = "results/cml_stop/cellphonedb/input_files/")
data.frame(diet_seurat@meta.data, cluster2 = diet_seurat$cluster) %>% filter(project %in% c("hsc") | timepoint == "12m") %>% initCellphonedb(seurat_object = diet_seurat, name = "tfr_hsc", sample_n = sample_n, folder = "results/cml_stop/cellphonedb/input_files/")
data.frame(diet_seurat@meta.data, cluster2 = diet_seurat$cluster) %>% filter(project %in% c("hsc") | timepoint == "relapse") %>% initCellphonedb(seurat_object = diet_seurat, name = "rp_hsc", sample_n = sample_n, folder = "results/cml_stop/cellphonedb/input_files/")

data.frame(diet_seurat@meta.data, cluster2 = diet_seurat$cluster) %>% filter(project %in% c("hsc") | timepoint == "baseline" & relapse != "Fast") %>% initCellphonedb(seurat_object = diet_seurat, name = "controller_0m_hsc", sample_n = sample_n, folder = "results/cml_stop/cellphonedb/input_files/")
data.frame(diet_seurat@meta.data, cluster2 = diet_seurat$cluster) %>% filter(project %in% c("hsc") | timepoint == "baseline" & relapse == "Fast") %>% initCellphonedb(seurat_object = diet_seurat, name = "fast_0m_hsc", sample_n = sample_n, folder = "results/cml_stop/cellphonedb/input_files/")
data.frame(diet_seurat@meta.data, cluster2 = diet_seurat$cluster) %>% filter(project %in% c("hsc") | timepoint == "baseline" & relapse == "Slow") %>% initCellphonedb(seurat_object = diet_seurat, name = "slow_0m_hsc", sample_n = sample_n, folder = "results/cml_stop/cellphonedb/input_files/")
data.frame(diet_seurat@meta.data, cluster2 = diet_seurat$cluster) %>% filter(project %in% c("hsc") | timepoint == "baseline" & relapse == "None") %>% initCellphonedb(seurat_object = diet_seurat, name = "none_0m_tki_hsc", sample_n = sample_n, folder = "results/cml_stop/cellphonedb/input_files/")


## sigf means
sigf_means_dg_hsc  <- fread("results/cml_stop/cellphonedb/out/dg_hsc/significant_means.txt")
sigf_means_tki_hsc <- fread("results/cml_stop/cellphonedb/out/tki_hsc/significant_means.txt")
sigf_means_tfr_hsc <- fread("results/cml_stop/cellphonedb/out/tfr_hsc/significant_means.txt")
sigf_means_rp_hsc  <- fread("results/cml_stop/cellphonedb/out/rp_hsc/significant_means.txt")

sigf_means_controller_0m  <- fread("results/cml_stop/cellphonedb/out/controller_0m_hsc//significant_means.txt")
sigf_means_slow_0m        <- fread("results/cml_stop/cellphonedb/out/slow_0m_hsc//significant_means.txt")
sigf_means_fast_0m        <- fread("results/cml_stop/cellphonedb/out/fast_0m_hsc//significant_means.txt")
sigf_means_none_0m        <- fread("results/cml_stop/cellphonedb/out/none_0m_tki_hsc//significant_means.txt")

sigf_means_dg_hsc  %>% getReadableCellphoneDb %>% fwrite("results/cml_stop/cellphonedb/out/dg_hsc/significant_means_human_readable.txt", sep = "\t", quote = F, row.names = F)
sigf_means_tki_hsc %>% getReadableCellphoneDb %>% fwrite("results/cml_stop/cellphonedb/out/tki_hsc/significant_means_human_readable.txt", sep = "\t", quote = F, row.names = F)
sigf_means_tfr_hsc %>% getReadableCellphoneDb %>% fwrite("results/cml_stop/cellphonedb/out/tfr_hsc/significant_means_human_readable.txt", sep = "\t", quote = F, row.names = F)
sigf_means_rp_hsc  %>% getReadableCellphoneDb %>% fwrite("results/cml_stop/cellphonedb/out/rp_hsc/significant_means_human_readable.txt", sep = "\t", quote = F, row.names = F)

sigf_means_controller_0m %>% getReadableCellphoneDb %>% fwrite("results/cml_stop/cellphonedb/out/controller_0m_hsc/significant_means_human_readable.txt", sep = "\t", quote = F, row.names = F)
sigf_means_slow_0m       %>% getReadableCellphoneDb %>% fwrite("results/cml_stop/cellphonedb/out/slow_0m_hsc/significant_means_human_readable.txt", sep = "\t", quote = F, row.names = F)
sigf_means_fast_0m       %>% getReadableCellphoneDb %>% fwrite("results/cml_stop/cellphonedb/out/fast_0m_hsc/significant_means_human_readable.txt", sep = "\t", quote = F, row.names = F)
sigf_means_none_0m       %>% getReadableCellphoneDb %>% fwrite("results/cml_stop/cellphonedb/out/none_0m_tki_hsc/significant_means_human_readable.txt", sep = "\t", quote = F, row.names = F)


## How many interactions are there?
sigf_means_dg_hsc_pairs   <- getPairAmounts(sigf_means_dg_hsc)  %>% filterRedundantPairs()
sigf_means_tki_hsc_pairs  <- getPairAmounts(sigf_means_tki_hsc)  %>% filterRedundantPairs() 
sigf_means_tfr_hsc_pairs  <- getPairAmounts(sigf_means_tfr_hsc)  %>% filterRedundantPairs() 
sigf_means_rp_hsc_pairs   <- getPairAmounts(sigf_means_rp_hsc)  %>% filterRedundantPairs() 

sigf_means_controller_0m_pairs <- getPairAmounts(sigf_means_controller_0m) %>% filterRedundantPairs()
sigf_means_slow_0m_pairs       <- getPairAmounts(sigf_means_slow_0m) %>% filterRedundantPairs() 
sigf_means_fast_0m_pairs       <- getPairAmounts(sigf_means_fast_0m) %>% filterRedundantPairs() 
sigf_means_none_0m_pairs       <- getPairAmounts(sigf_means_none_0m) %>% filterRedundantPairs() 



## Visualize; dg-phase
p <- cast(sigf_means_dg_hsc_pairs, cluster1 ~ cluster2, sum)
rownames(p) <- colnames(p)[-1]
p <- p[,-1]
colnames(p) <- substr(colnames(p), 3, nchar(colnames(p)))
rownames(p) <- substr(rownames(p), 3, nchar(rownames(p)))
rownames(p) <- gsub("\\_", "\\/", rownames(p))
colnames(p) <- gsub("\\_", "\\/", colnames(p))

p <- p[!rownames(p) %in% getClusterPhenotypesLSC(as.character(0:7)), colnames(p) %in% getClusterPhenotypesLSC(as.character(0:7))]
q <- p %>% pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(seq(0, 60, 0.6))), breaks = seq(0, 60, 0.6), display_numbers=T,number_color="gray",number_format = "%.0f")
ggsave(plot = q, "results/cml_stop/cellphonedb/heatmap_dg_hsc.png", width = 6, height = 6)

melt(p) %>% left_join(melt(p) %>% group_by(variable) %>% summarise(median = median(value))) %>% 
  ggplot(aes(reorder(variable,median),value,fill=variable)) + geom_boxplot(outlier.shape = NA) + ggpubr::rotate_x_text(angle = 45) + labs(x = "", y = "nInteractions") + 
  scale_fill_manual(values = getPalette5(8)) + theme(legend.position = "none") + geom_jitter(size = 0.5) + ggpubr::stat_compare_means()
ggsave("results/cml_stop/cellphonedb/box_interactions_dg.pdf", width = 5, height = 4)

melt(p) %>% mutate(bcrabl = ifelse(grepl("BCRABL1+", variable), "BCRABL1+", "BCRABL1-")) %>% 
  ggplot(aes(bcrabl,value,fill=bcrabl)) + geom_boxplot(outlier.shape = NA) + ggpubr::rotate_x_text(angle = 45) + labs(x = "", y = "nInteractions") + 
  scale_fill_manual(values = getPalette3(4)) + theme(legend.position = "none") + geom_jitter(size = 0.5) + ggpubr::stat_compare_means()
ggsave("results/cml_stop/cellphonedb/box_bcrabl_interactions_dg.pdf", width = 5, height = 4)



## Visualize; tki_hsc
p <- cast(sigf_means_tki_hsc_pairs, cluster1 ~ cluster2, sum)
rownames(p) <- colnames(p)[-1]
p <- p[,-1]
colnames(p) <- substr(colnames(p), 3, nchar(colnames(p)))
rownames(p) <- substr(rownames(p), 3, nchar(rownames(p)))
rownames(p) <- gsub("\\_", "\\/", rownames(p))
colnames(p) <- gsub("\\_", "\\/", colnames(p))

p <- p[!rownames(p) %in% getClusterPhenotypesLSC(as.character(0:7)), colnames(p) %in% getClusterPhenotypesLSC(as.character(0:7))]
q <- p %>% pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(seq(0, 60, 0.6))), breaks = seq(0, 60, 0.6), display_numbers=T,number_color="gray",number_format = "%.0f")
ggsave(plot = q, "results/cml_stop/cellphonedb/heatmap_tki.png", width = 6, height = 6)

melt(p) %>% left_join(melt(p) %>% group_by(variable) %>% summarise(median = median(value))) %>% 
  ggplot(aes(reorder(variable,median),value,fill=variable)) + geom_boxplot(outlier.shape = NA) + ggpubr::rotate_x_text(angle = 45) + labs(x = "", y = "nInteractions") + 
  scale_fill_manual(values = getPalette5(8)) + theme(legend.position = "none") + geom_jitter(size = 0.5) + ggpubr::stat_compare_means()
ggsave("results/cml_stop/cellphonedb/box_interactions_tki.pdf", width = 5, height = 4)

melt(p) %>% mutate(bcrabl = ifelse(grepl("BCRABL1+", variable), "BCRABL1+", "BCRABL1-")) %>% 
  ggplot(aes(bcrabl,value,fill=bcrabl)) + geom_boxplot(outlier.shape = NA) + ggpubr::rotate_x_text(angle = 45) + labs(x = "", y = "nInteractions") + 
  scale_fill_manual(values = getPalette3(4)) + theme(legend.position = "none") + geom_jitter(size = 0.5) + ggpubr::stat_compare_means()
ggsave("results/cml_stop/cellphonedb/box_bcrabl_interactions_tki.pdf", width = 5, height = 4)


## Visualize; tfr
p <- cast(sigf_means_tfr_hsc_pairs, cluster1 ~ cluster2, sum)
rownames(p) <- colnames(p)[-1]
p <- p[,-1]
colnames(p) <- substr(colnames(p), 3, nchar(colnames(p)))
rownames(p) <- substr(rownames(p), 3, nchar(rownames(p)))
rownames(p) <- gsub("\\_", "\\/", rownames(p))
colnames(p) <- gsub("\\_", "\\/", colnames(p))

p <- p[!rownames(p) %in% getClusterPhenotypesLSC(as.character(0:7)), colnames(p) %in% getClusterPhenotypesLSC(as.character(0:7))]
q <- p %>% pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(seq(0, 60, 0.6))), breaks = seq(0, 60, 0.6), display_numbers=T,number_color="gray",number_format = "%.0f")
ggsave(plot = q, "results/cml_stop/cellphonedb/heatmap_tfr.png", width = 6, height = 6)

melt(p) %>% left_join(melt(p) %>% group_by(variable) %>% summarise(median = median(value))) %>% 
  ggplot(aes(reorder(variable,median),value,fill=variable)) + geom_boxplot(outlier.shape = NA) + ggpubr::rotate_x_text(angle = 45) + labs(x = "", y = "nInteractions") + 
  scale_fill_manual(values = getPalette5(8)) + theme(legend.position = "none") + geom_jitter(size = 0.5) + ggpubr::stat_compare_means()
ggsave("results/cml_stop/cellphonedb/box_interactions_tfr.pdf", width = 5, height = 4)

melt(p) %>% mutate(bcrabl = ifelse(grepl("BCRABL1+", variable), "BCRABL1+", "BCRABL1-")) %>% 
  ggplot(aes(bcrabl,value,fill=bcrabl)) + geom_boxplot(outlier.shape = NA) + ggpubr::rotate_x_text(angle = 45) + labs(x = "", y = "nInteractions") + 
  scale_fill_manual(values = getPalette3(4)) + theme(legend.position = "none") + geom_jitter(size = 0.5) + ggpubr::stat_compare_means()
ggsave("results/cml_stop/cellphonedb/box_bcrabl_interactions_tfr.pdf", width = 5, height = 4)



## Visualize; rp
p <- cast(sigf_means_rp_hsc_pairs, cluster1 ~ cluster2, sum)
rownames(p) <- colnames(p)[-1]
p <- p[,-1]
colnames(p) <- substr(colnames(p), 3, nchar(colnames(p)))
rownames(p) <- substr(rownames(p), 3, nchar(rownames(p)))
rownames(p) <- gsub("\\_", "\\/", rownames(p))
colnames(p) <- gsub("\\_", "\\/", colnames(p))

p <- p[!rownames(p) %in% getClusterPhenotypesLSC(as.character(0:7)), colnames(p) %in% getClusterPhenotypesLSC(as.character(0:7))]
q <- p %>% pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(seq(0, 60, 0.6))), breaks = seq(0, 60, 0.6), display_numbers=T,number_color="gray",number_format = "%.0f")
ggsave(plot = q, "results/cml_stop/cellphonedb/heatmap_rp.png", width = 6, height = 6)

melt(p) %>% left_join(melt(p) %>% group_by(variable) %>% summarise(median = median(value))) %>% 
  ggplot(aes(reorder(variable,median),value,fill=variable)) + geom_boxplot(outlier.shape = NA) + ggpubr::rotate_x_text(angle = 45) + labs(x = "", y = "nInteractions") + 
  scale_fill_manual(values = getPalette5(8)) + theme(legend.position = "none") + geom_jitter(size = 0.5) + ggpubr::stat_compare_means()
ggsave("results/cml_stop/cellphonedb/box_interactions_rp.pdf", width = 5, height = 4)

melt(p) %>% mutate(bcrabl = ifelse(grepl("BCRABL1+", variable), "BCRABL1+", "BCRABL1-")) %>% 
  ggplot(aes(bcrabl,value,fill=bcrabl)) + geom_boxplot(outlier.shape = NA) + ggpubr::rotate_x_text(angle = 45) + labs(x = "", y = "nInteractions") + 
  scale_fill_manual(values = getPalette3(4)) + theme(legend.position = "none") + geom_jitter(size = 0.5) + ggpubr::stat_compare_means()
ggsave("results/cml_stop/cellphonedb/box_bcrabl_interactions_rp.pdf", width = 5, height = 4)






## Visualize; baseline disease controller
p <- cast(sigf_means_controller_0m_pairs, cluster1 ~ cluster2, sum)
rownames(p) <- colnames(p)[-1]
p <- p[,-1]
colnames(p) <- substr(colnames(p), 3, nchar(colnames(p)))
rownames(p) <- substr(rownames(p), 3, nchar(rownames(p)))
rownames(p) <- gsub("\\_", "\\/", rownames(p))
colnames(p) <- gsub("\\_", "\\/", colnames(p))

p <- p[!rownames(p) %in% getClusterPhenotypesLSC(as.character(0:7)), colnames(p) %in% getClusterPhenotypesLSC(as.character(0:7))]
q <- p %>% pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(seq(0, 60, 0.6))), breaks = seq(0, 60, 0.6), display_numbers=T,number_color="gray",number_format = "%.0f")
ggsave(plot = q, "results/cml_stop/cellphonedb/heatmap_0m_ctrl.png", width = 6, height = 6)

melt(p) %>% left_join(melt(p) %>% group_by(variable) %>% summarise(median = median(value))) %>% 
  ggplot(aes(reorder(variable,median),value,fill=variable)) + geom_boxplot(outlier.shape = NA) + ggpubr::rotate_x_text(angle = 45) + labs(x = "", y = "nInteractions") + 
  scale_fill_manual(values = getPalette5(8)) + theme(legend.position = "none") + geom_jitter(size = 0.5) + ggpubr::stat_compare_means()
ggsave("results/cml_stop/cellphonedb/box_interactions_0m_ctrl.pdf", width = 5, height = 4)






## Visualize; baseline slow
p <- cast(sigf_means_slow_0m_pairs, cluster1 ~ cluster2, sum)
rownames(p) <- colnames(p)[-1]
p <- p[,-1]
colnames(p) <- substr(colnames(p), 3, nchar(colnames(p)))
rownames(p) <- substr(rownames(p), 3, nchar(rownames(p)))
rownames(p) <- gsub("\\_", "\\/", rownames(p))
colnames(p) <- gsub("\\_", "\\/", colnames(p))

p <- p[!rownames(p) %in% getClusterPhenotypesLSC(as.character(0:7)), colnames(p) %in% getClusterPhenotypesLSC(as.character(0:7))]
q <- p %>% pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(seq(0, 60, 0.6))), breaks = seq(0, 60, 0.6), display_numbers=T,number_color="gray",number_format = "%.0f")
ggsave(plot = q, "results/cml_stop/cellphonedb/heatmap_0m_slow.png", width = 6, height = 6)

melt(p) %>% left_join(melt(p) %>% group_by(variable) %>% summarise(median = median(value))) %>% 
  ggplot(aes(reorder(variable,median),value,fill=variable)) + geom_boxplot(outlier.shape = NA) + ggpubr::rotate_x_text(angle = 45) + labs(x = "", y = "nInteractions") + 
  scale_fill_manual(values = getPalette5(8)) + theme(legend.position = "none") + geom_jitter(size = 0.5) + ggpubr::stat_compare_means()
ggsave("results/cml_stop/cellphonedb/box_interactions_0m_slow.pdf", width = 5, height = 4)





## Visualize; baseline fast
p <- cast(sigf_means_fast_0m_pairs, cluster1 ~ cluster2, sum)
rownames(p) <- colnames(p)[-1]
p <- p[,-1]
colnames(p) <- substr(colnames(p), 3, nchar(colnames(p)))
rownames(p) <- substr(rownames(p), 3, nchar(rownames(p)))
rownames(p) <- gsub("\\_", "\\/", rownames(p))
colnames(p) <- gsub("\\_", "\\/", colnames(p))

p <- p[!rownames(p) %in% getClusterPhenotypesLSC(as.character(0:7)), colnames(p) %in% getClusterPhenotypesLSC(as.character(0:7))]
q <- p %>% pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(seq(0, 60, 0.6))), breaks = seq(0, 60, 0.6), display_numbers=T,number_color="gray",number_format = "%.0f")
ggsave(plot = q, "results/cml_stop/cellphonedb/heatmap_0m_fast.png", width = 6, height = 6)

melt(p) %>% left_join(melt(p) %>% group_by(variable) %>% summarise(median = median(value))) %>% 
  ggplot(aes(reorder(variable,median),value,fill=variable)) + geom_boxplot(outlier.shape = NA) + ggpubr::rotate_x_text(angle = 45) + labs(x = "", y = "nInteractions") + 
  scale_fill_manual(values = getPalette5(8)) + theme(legend.position = "none") + geom_jitter(size = 0.5) + ggpubr::stat_compare_means()
ggsave("results/cml_stop/cellphonedb/box_interactions_0m_fast.pdf", width = 5, height = 4)




## Visualize; baseline none
p <- cast(sigf_means_none_0m_pairs, cluster1 ~ cluster2, sum)
rownames(p) <- colnames(p)[-1]
p <- p[,-1]
colnames(p) <- substr(colnames(p), 3, nchar(colnames(p)))
rownames(p) <- substr(rownames(p), 3, nchar(rownames(p)))
rownames(p) <- gsub("\\_", "\\/", rownames(p))
colnames(p) <- gsub("\\_", "\\/", colnames(p))

p <- p[!rownames(p) %in% getClusterPhenotypesLSC(as.character(0:7)), colnames(p) %in% getClusterPhenotypesLSC(as.character(0:7))]
q <- p %>% pheatmap::pheatmap(color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(seq(0, 60, 0.6))), breaks = seq(0, 60, 0.6), display_numbers=T,number_color="gray",number_format = "%.0f")
ggsave(plot = q, "results/cml_stop/cellphonedb/heatmap_0m_none.png", width = 6, height = 6)

melt(p) %>% left_join(melt(p) %>% group_by(variable) %>% summarise(median = median(value))) %>% 
  ggplot(aes(reorder(variable,median),value,fill=variable)) + geom_boxplot(outlier.shape = NA) + ggpubr::rotate_x_text(angle = 45) + labs(x = "", y = "nInteractions") + 
  scale_fill_manual(values = getPalette5(8)) + theme(legend.position = "none") + geom_jitter(size = 0.5) + ggpubr::stat_compare_means()
ggsave("results/cml_stop/cellphonedb/box_interactions_0m_none.pdf", width = 5, height = 4)










## Compare fast and ctrl at base line
fast <- cast(sigf_means_fast_0m_pairs, cluster1 ~ cluster2, sum)
rownames(fast) <- colnames(fast)[-1]
fast <- fast[,-1]
colnames(fast) <- substr(colnames(fast), 3, nchar(colnames(fast)))
rownames(fast) <- substr(rownames(fast), 3, nchar(rownames(fast)))
rownames(fast) <- gsub("\\_", "\\/", rownames(fast))
colnames(fast) <- gsub("\\_", "\\/", colnames(fast))
# fast <- fast %>% mutate("20 CD4 Acute Act" = 0, "23 low quality" = 0) 


# ctrl <- cast(sigf_means_controller_0m_pairs, cluster1 ~ cluster2, sum)
ctrl <- cast(sigf_means_none_0m_pairs, cluster1 ~ cluster2, sum)
# ctrl <- cast(sigf_means_slow_0m_pairs, cluster1 ~ cluster2, sum)

rownames(ctrl) <- colnames(ctrl)[-1]
ctrl <- ctrl[,-1]
colnames(ctrl) <- substr(colnames(ctrl), 3, nchar(colnames(ctrl)))
rownames(ctrl) <- substr(rownames(ctrl), 3, nchar(rownames(ctrl)))
rownames(ctrl) <- gsub("\\_", "\\/", rownames(ctrl))
colnames(ctrl) <- gsub("\\_", "\\/", colnames(ctrl))


cmn.names <- rownames(ctrl)[rownames(ctrl) %in% rownames(fast)]
fast <- fast[cmn.names,cmn.names]
ctrl <- ctrl[cmn.names,cmn.names]

colnames(ctrl) == colnames(fast)
rownames(ctrl) == rownames(fast)

log2fc <- log2(ctrl/fast)
log2fc[is.na(log2fc)] <- 0
p <- log2fc
p <- p[!rownames(p) %in% getClusterPhenotypesLSC(as.character(0:7)), colnames(p) %in% getClusterPhenotypesLSC(as.character(0:7))]

pheatmap::pheatmap(p, color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(seq(-2, 2, 0.1))), breaks = seq(-2, 2, 0.1))




cols.to.keep <- colnames(sigf_means_fast_0m)[as.numeric(c(1:12, grep("NK adaptive", colnames(sigf_means_fast_0m))))]
cols.to.keep <- grep("BCRABL", cols.to.keep, value = T)
cols.to.keep <- c(colnames(sigf_means_fast_0m)[as.numeric(c(1:12))], cols.to.keep[grep("\\|1 15 NK adaptive", cols.to.keep, invert = F)])
fast_nk <- sigf_means_fast_0m %>% select(cols.to.keep)
fast_nk <- fast_nk[rowSums(fast_nk[,-c(1:12)], na.rm = T) > 0, ]

cols.to.keep <- colnames(sigf_means_none_0m)[as.numeric(c(1:12, grep("NK adaptive", colnames(sigf_means_none_0m))))]
cols.to.keep <- grep("BCRABL", cols.to.keep, value = T)
cols.to.keep <- c(colnames(sigf_means_none_0m)[as.numeric(c(1:12))], cols.to.keep[grep("\\|1 15 NK adaptive", cols.to.keep, invert = F)])
none_nk <- sigf_means_none_0m %>% select(cols.to.keep)
none_nk <- none_nk[rowSums(none_nk[,-c(1:12)], na.rm = T) > 0, ]

fast_nk_uniq <- fast_nk[!fast_nk$interacting_pair %in% none_nk$interacting_pair, ]
View(fast_nk_uniq)





#### Focus on NK cells
cols.to.keep <- colnames(sigf_means_fast_0m)[intersect(grep("BCRABL1\\-", colnames(sigf_means_fast_0m)), grep("NK", colnames(sigf_means_fast_0m)))]
cols.to.keep <- cols.to.keep[substr(cols.to.keep,5,5) == "B"]
cols.to.keep <- c(colnames(sigf_means_fast_0m)[as.numeric(c(1:12))], cols.to.keep)

dg_nk <- sigf_means_dg_hsc %>% select(cols.to.keep)
dg_nk <- dg_nk[rowSums(dg_nk[,-c(1:12)], na.rm = T) > 0, ]

# dg_nk <- dg_nk[,colSums(dg_nk[,-c(1:12)], na.rm = T) > 0]

dg_nk_heatmap <- dg_nk %>% select(-c(1:12))
dg_nk_heatmap[is.na(dg_nk_heatmap)] <- 0
dg_nk_heatmap_neg <- log10(dg_nk_heatmap+1) %>% as.data.frame()
rownames(dg_nk_heatmap_neg) <- dg_nk$interacting_pair
p <- pheatmap::pheatmap(dg_nk_heatmap_neg, color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(seq(0, 1.5, 0.01))), breaks = seq(0, 1.5, 0.01))
ggsave(plot = print(p), "results/cml_stop/cellphonedb/heatmap_dg_nk_bcrablneg.pdf", width = 5, height = 14)


cols.to.keep <- colnames(sigf_means_fast_0m)[intersect(grep("BCRABL1\\+", colnames(sigf_means_fast_0m)), grep("NK", colnames(sigf_means_fast_0m)))]
cols.to.keep <- cols.to.keep[substr(cols.to.keep,5,5) == "B"]
cols.to.keep <- c(colnames(sigf_means_fast_0m)[as.numeric(c(1:12))], cols.to.keep)

dg_nk <- sigf_means_dg_hsc %>% select(cols.to.keep)
dg_nk <- dg_nk[rowSums(dg_nk[,-c(1:12)], na.rm = T) > 0, ]

# dg_nk <- dg_nk[,colSums(dg_nk[,-c(1:12)], na.rm = T) > 0]

dg_nk_heatmap <- dg_nk %>% select(-c(1:12))
dg_nk_heatmap[is.na(dg_nk_heatmap)] <- 0
dg_nk_heatmap_pos <- log10(dg_nk_heatmap+1) %>% as.data.frame()
rownames(dg_nk_heatmap_pos) <- dg_nk$interacting_pair
p <- pheatmap::pheatmap(dg_nk_heatmap_pos, color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(seq(0, 1.5, 0.01))), breaks = seq(0, 1.5, 0.01))
ggsave(plot = print(p), "results/cml_stop/cellphonedb/heatmap_dg_nk_bcrablpos.pdf", width = 5, height = 14)



cmn.names <- intersect(rownames(dg_nk_heatmap_pos), rownames(dg_nk_heatmap_neg))
a <- dg_nk_heatmap_pos[cmn.names,]
b <- dg_nk_heatmap_neg[cmn.names,]

a <- pheatmap::pheatmap(a, cluster_rows = F, color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(seq(0, 1, 0.01))), breaks = seq(0, 1, 0.01))
b <- pheatmap::pheatmap(b, cluster_rows = F, color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(seq(0, 1, 0.01))), breaks = seq(0, 1, 0.01))

ggsave(plot = print(a), "results/cml_stop/cellphonedb/heatmap_dg_nk_bcrablpos_cmn.pdf", width = 5, height = 14)
ggsave(plot = print(b), "results/cml_stop/cellphonedb/heatmap_dg_nk_bcrablneg_cmn.pdf", width = 5, height = 14)
