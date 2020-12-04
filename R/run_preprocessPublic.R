

####### Public data sets
folders        <- list.dirs("data/scRNAseq/", recursive = T)[-c(1,17)]
scrnaseq_files <- lapply(folders, function(x){message(getSeuratName(x)); Read10X(data.dir = x) %>% CreateSeuratObject(project = getSeuratName(x), min.cells = 3, min.features = 200)})
public_seurat     <- merge(scrnaseq_files[[1]], scrnaseq_files[-1], add.cell.ids = getSeuratName(folders))

## Basic QC
public_seurat  <- PercentageFeatureSet(public_seurat, pattern = "^MT-", col.name = "percent.mt")
public_seurat  <- PercentageFeatureSet(public_seurat, pattern = "^RP", col.name = "percent.ribo")
public_seurat  <- PercentageFeatureSet(public_seurat, features = cycle.genes, col.name = "percent.cycle")
public_seurat@meta.data$barcode   <- colnames(public_seurat)

## Preprocess data
clonality_genes <- getClonalityGenes(public_seurat)
unwanted_genes  <- getUnwantedGenes(public_seurat)

public_seurat <- public_seurat %>% getQC()
public_seurat <- public_seurat %>% preprocessSeurat(cells.to.use = colnames(public_seurat))

## Do scVI
public_seurat %>% getScviInput(folder = "results/scvi/input_files/")

## Get scVI results
latents    <- fread("results/scvi/results/cml_latents.csv")
public_seurat <- public_seurat %>% putLatentsSeurat(latent = latents)
public_seurat <- public_seurat %>% getLatentClustering() %>% fixSeurat()

## Decide on clustering
public_seurat %>% plotClustering()
Idents(public_seurat) <- public_seurat$RNA_snn_res.0.5 %>% as.character() %>% as.numeric() %>% as.factor()
public_seurat$cluster <- Idents(public_seurat)

## Add TCRab
tot_barcode <- fread("data/scRNAseq+TCRseq/preprocessed/tot_barcode.txt")
public_seurat <- mergeTCRtoSeurat(seurat_object = public_seurat, tcr_df = tot_barcode)

## Get SingleR predictions; omit predictions from cell types rare than 10 cells
public_seurat             <- public_seurat %>% getSingler()
relevant_hpca_clusters <- public_seurat@meta.data %>% group_by(singler_hpca_pred) %>% summarise(n=n()) %>% filter(n>=10) %>% pull(singler_hpca_pred)
relevant_blue_clusters <- public_seurat@meta.data %>% group_by(singler_blueprint_pred) %>% summarise(n=n()) %>% filter(n>=10) %>% pull(singler_blueprint_pred)

public_seurat$singler_hpca_pred      <- ifelse(public_seurat$singler_hpca_pred %in% relevant_hpca_clusters, public_seurat$singler_hpca_pred, "rare")
public_seurat$singler_blueprint_pred <- ifelse(public_seurat$singler_blueprint_pred %in% relevant_blue_clusters, public_seurat$singler_blueprint_pred, "rare")
saveRDS(public_seurat, "results/public_seurat.rds")


## Analyze

## Find DE-genes between all the clusters
all_markers <- FindAllMarkers(public_seurat, test.use = "t", verbose = T, max.cells.per.ident = 1e3)
fwrite(all_markers, "results/manuscript/public/all_markers_1e3.txt", sep = "\t", quote = F, row.names = F)









### Find project specific clusters
patient_df <- public_seurat@meta.data %>% filter(cluster != "5 CLL cells") %>% group_by(orig.ident, project) %>% summarise(n = n()) %>% dplyr::select(-n)
df         <- public_seurat@meta.data %>% filter(cluster != "5 CLL cells") %>% group_by(orig.ident, cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% left_join(patient_df)

p.df <- lapply(unique(df$cluster), FUN = function(x){
  y <- df %>% filter(cluster == x)
  if(length(unique(y$project)) > 2){
    kruskal.test(prop~project, data = y) %>% broom::tidy() %>% mutate(cluster = x)
  }
}) %>% rbindlist() %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% arrange(p.adj)

fwrite(p.df, "results/manuscript/public/cml_kruskal_p_df.txt", sep = "\t", quote = F, row.names = F)

patient_df <- public_seurat@meta.data %>% filter(cluster != "5 CLL cells") %>% mutate(project = ifelse(project %in% c("CML on TKI", "CML dg", "CML off TKI"), "CML", "other")) %>% group_by(orig.ident, project) %>% summarise(n = n()) %>% dplyr::select(-n)
df         <- public_seurat@meta.data %>% filter(cluster != "5 CLL cells") %>% group_by(orig.ident, cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% left_join(patient_df)

p.df <- lapply(unique(df$cluster), FUN = function(x){
  message(x)
  y <- df %>% filter(cluster == x)
  if(length(unique(y$project)) == 2){
    median.x = median(subset(y, project == "CML")$prop)
    median.y = median(subset(y, project != "CML")$prop)
    wilcox.test(prop~project, data = y) %>% broom::tidy() %>% mutate(cluster = x, median.cml = median.x, median.other = median.y, log2fc = log2(median.x/median.y))
  }
}) %>% rbindlist() %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% arrange(p.adj)

fwrite(p.df, "results/manuscript/public/cml_p_df.txt", sep = "\t", quote = F, row.names = F)

#### CMV specific clusters
patient_df <- public_seurat@meta.data %>% filter(cluster != "5 CLL cells") %>% group_by(orig.ident, cmv_status) %>% summarise(n = n()) %>% dplyr::select(-n)
df <- public_seurat@meta.data %>% filter(cluster != "5 CLL cells") %>% group_by(orig.ident, cluster) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% left_join(patient_df)

p.df <- lapply(unique(df$cluster), FUN = function(x){
  message(x)
  y <- df %>% filter(cluster == x)# %>% filter(cmv_status != "CMV unknown")
  if(length(unique(y$cmv_status)) > 2){
    kruskal.test(prop~cmv_status, data = y) %>% broom::tidy() %>% mutate(cluster = x)
  }
}) %>% rbindlist() %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% arrange(p.adj)

fwrite(p.df, "results/manuscript/public/cml_cmv_p_df.txt", sep = "\t", quote = F, row.names = F)
