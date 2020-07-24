
## Run slingshot
require(slingshot)
require(SingleCellExperiment)
require(SummarizedExperiment)

diet_nk_seurat <- readRDS("results/diet_nk_seurat.rds")
diet_nk_seurat$orig.clusters <- Idents(diet_nk_seurat)
diet_nk_sce <- as.SingleCellExperiment(diet_nk_seurat)

diet_nk_sling <- runSlingshot(diet_nk_sce, cluster_with_earliset_timepoint = "5 CD56 bright", reducedDim = "LATENT_UMAP")
diet_nk_sling <- readRDS("results/diet_nk_sling.rds")

slingshot_object=diet_nk_sling
reducedDim="LATENT_UMAP"

cololors                <- slingshot_object$orig.clusters %>% extractClusterNumber() %>% as.numeric()
cololors[cololors == 0] <- getPalette5(8)[1]
cololors[cololors == 1] <- getPalette5(8)[2]
cololors[cololors == 2] <- getPalette5(8)[3]
cololors[cololors == 3] <- getPalette5(8)[4]
cololors[cololors == 4] <- getPalette5(8)[5]
cololors[cololors == 5] <- getPalette5(8)[6]
cololors[cololors == 6] <- getPalette5(8)[7]
cololors[cololors == 7] <- getPalette5(8)[8]

curve_cols <- c("black", "darkred", "darkblue", "darkolivegreen4")
sling_curve <- SlingshotDataSet(slingshot_object)

pdf("results/figure1/umap_sling_nk.pdf", width = 5, height = 5)
plot(reducedDims(slingshot_object)[[reducedDim]], col = cololors, pch = 16, asp = 1)
for(i in 1:2){
  lines(sling_curve@curves[[i]], lwd = 4, col = curve_cols[i])
}
dev.off()


png("results/figure1/umap_sling_nk.png", width = 5, height = 5, units = "in", res = 1024)
plot(reducedDims(slingshot_object)[[reducedDim]], col = cololors, pch = 16, asp = 1)
for(i in 1:2){
  lines(sling_curve@curves[[i]], lwd = 4, col = curve_cols[i])
}
dev.off()

### Trajectories in different diseases
diet_nk_seurat$orig.clusters <- Idents(diet_nk_seurat)
diet_nk_sce <- as.SingleCellExperiment(diet_nk_seurat)

diet_nk_sling <- runSlingshot(diet_nk_sce, cluster_with_earliset_timepoint = "5 CD56 bright", reducedDim = "LATENT_UMAP")
diet_nk_sling <- readRDS("results/diet_nk_sling.rds")

pdf("results/figure1/umap_sling_nk.pdf", width = 5, height = 5)
plot(reducedDims(diet_nk_sling)[[reducedDim]], col = cololors, pch = 16, asp = 1)
for(i in 1:2){
  lines(sling_curve@curves[[i]], lwd = 4, col = curve_cols[i])
}
dev.off()


## In differet diseases
getSlingNk <- function(project_temp){
  
  cells.to.keep       <- diet_nk_seurat@meta.data %>% filter(project == project_temp) %>% pull(barcode)
  project_seurat      <- subset(diet_nk_seurat, cells = cells.to.keep) %>% as.SingleCellExperiment()
  proje_nk_sling      <- runSlingshot(project_seurat, cluster_with_earliset_timepoint = "5 CD56 bright", reducedDim = "LATENT_UMAP")
  return(proje_nk_sling)
  
}

cml_dg_sling      <- getSlingNk(project_temp = "CML dg")
cml_on_tki_sling  <- getSlingNk(project_temp = "CML on TKI")
cml_off_tki_sling <- getSlingNk(project_temp = "CML off TKI")
nsclc_sling       <- getSlingNk(project_temp = "NSCLC from blood")
rcc_sling         <- getSlingNk(project_temp = "RCC from blood")
cll_sling         <- getSlingNk(project_temp = "CLL")

cml_dg_curve      <- SlingshotDataSet(cml_dg_sling)
cml_on_tki_curve  <- SlingshotDataSet(cml_on_tki_sling)
cml_off_tki_curve <- SlingshotDataSet(cml_off_tki_sling)
nsclc_curve       <- SlingshotDataSet(nsclc_sling)
rcc_curve         <- SlingshotDataSet(rcc_sling)
cll_curve         <- SlingshotDataSet(cll_sling)

png("results/manuscript/public/umap_sling_nk_disease.png", width = 1064*2/3, height = 1064*2/3*2/3, res = 80)

par(mfrow=c(2,3))

plot(reducedDims(diet_nk_sling)[[reducedDim]], col = cololors, pch = 16, asp = 1, main = "CML dg")
for(i in 1:3){lines(cml_dg_curve@curves[[i]], lwd = 4, col = curve_cols[1])}

plot(reducedDims(diet_nk_sling)[[reducedDim]], col = cololors, pch = 16, asp = 1, main = "CML on-TKI")
for(i in 1:3){lines(cml_on_tki_curve@curves[[i]], lwd = 4, col = curve_cols[2])}

plot(reducedDims(diet_nk_sling)[[reducedDim]], col = cololors, pch = 16, asp = 1, main = "CML off-TKI")
for(i in 1:3){lines(cml_off_tki_curve@curves[[i]], lwd = 4, col = curve_cols[3])}

plot(reducedDims(diet_nk_sling)[[reducedDim]], col = cololors, pch = 16, asp = 1, main = "NSCLC")
for(i in 1:1){lines(nsclc_curve@curves[[i]], lwd = 4, col = curve_cols[4])}

plot(reducedDims(diet_nk_sling)[[reducedDim]], col = cololors, pch = 16, asp = 1, main = "RCC")
for(i in 1:3){lines(rcc_curve@curves[[i]], lwd = 4, col = curve_cols[5])}

plot(reducedDims(diet_nk_sling)[[reducedDim]], col = cololors, pch = 16, asp = 1, main = "CLL")
for(i in 1:3){lines(cll_curve@curves[[i]], lwd = 4, col = curve_cols[6])}

dev.off()

curve_cols <- c("black", "darkred", "darkblue", "darkolivegreen4", "salmon", "dodgerblue")
pdf("results/figure1/umap_sling_nk.pdf", width = 5, height = 5)

plot(reducedDims(diet_nk_sling)[[reducedDim]], col = cololors, pch = 16, asp = 1)
for(i in 1){lines(cml_dg_curve@curves[[i]], lwd = 4, col = curve_cols[1])}
for(i in 1){lines(cml_on_tki_curve@curves[[i]], lwd = 4, col = curve_cols[2])}
for(i in 1){lines(cml_off_tki_curve@curves[[i]], lwd = 4, col = curve_cols[3])}
for(i in 1){lines(nsclc_curve@curves[[i]], lwd = 4, col = curve_cols[4])}
for(i in 1){lines(rcc_curve@curves[[i]], lwd = 4, col = curve_cols[5])}
for(i in 1){lines(cll_curve@curves[[i]], lwd = 4, col = curve_cols[6])}

dev.off()

##### Nk cells following cessation
cells.to.keep <- diet_nk_seurat@meta.data %>% filter(project == "CML on TKI") %>% pull(barcode)
diet_nk_tki_seurat <- subset(diet_nk_seurat, cells = cells.to.keep) %>% getLatentUMAP() %>% fixSeurat()

diet_nk_tki_seurat@meta.data %>%
  mutate(relapse   = factor(relapse, levels = c("None", "Slow", "Fast"))) %>%
  mutate(timepoint = factor(as.character(timepoint), levels = c("baseline", "6m", "12m", "relapse"))) %>%
  group_by(timepoint, relapse, cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% # ungroup() %>%
  
  ggplot(aes(as.factor(timepoint),prop,group=relapse,color=relapse)) + geom_path(lwd = 1.5) + facet_wrap(~cluster, scales = "free_y") + theme_bw(base_size = 12) + facets_nice + scale_color_manual(values = getPalette3(4)[c(1,3,2)]) +
  labs(x = "", y = "prop") + ggpubr::rotate_x_text(45)
ggsave("results/figure1/line_nk.pdf", width = 5, height = 5)
