
### Diagnostic and stop-TKI CML samples
message("Create seurat-object from individual seurat-objects ...")

folders        <- list.dirs("data/scRNAseq/", recursive = T) %>% grep(pattern = "filtered_feature_bc_matrix", value = T) # %>% grep(pattern = "petti", value = T)
scrnaseq_files <- lapply(folders, function(x){message(x); Read10X(data.dir = x) %>% CreateSeuratObject(project = extractSeuratName(x), min.cells = 3, min.features = 200)})
cml_seurat     <- merge(scrnaseq_files[[1]], scrnaseq_files[-1], add.cell.ids = extractSeuratName(folders))

## Basic QC
cycle.genes  <- c("ANLN", "ASPM","BIRC5","CCNA2","CCNB1","CCNB2","CCND1","CD63","CDC20","CDCA8","CDKN3","CENPE","CENPF",
                  "CEP55","CKAP2L","DLGAP5","FOXM1","GTSE1","H2AFZ","HIST1H1B", "HIST1H1C", "HIST1H1D", "HIST1H1E", "HIST1H2AJ",
                  "HIST1H4C", "HJURP", "HMGB1", "HMGB2", "HMMR", "KIF11", "KIF14", "KIF15", "KIF2C", "LMNA",
                  "MCM3", "MKI67", "NCAPG", "NUSAP1", "PCNA", "PLK1", "PRC1", "RRM2", "SMC4", "STMN1", "TK1", "TOP2A", "TPX2", "TUBA1B",
                  "TUBB", "TYMS", "UBE2C")

cml_seurat  <- PercentageFeatureSet(cml_seurat, pattern = "^MT-", col.name = "percent.mt")
cml_seurat  <- PercentageFeatureSet(cml_seurat, pattern = "^RP", col.name = "percent.ribo")
cml_seurat  <- PercentageFeatureSet(cml_seurat, features = cycle.genes, col.name = "percent.cycle")
cml_seurat@meta.data$barcode   <- colnames(cml_seurat)
cml_seurat  <- getQC(cml_seurat)



## Wu (Nature 2020), 5' prime end solution on blood CD45+ human cancers
n <- 17
x = folders

folders        <- list.dirs("public_data/Wu/", recursive = T)[-1] # %>% grep(pattern = "filtered_feature_bc_matrix", value = T) # %>% grep(pattern = "petti", value = T)
scrnaseq_files <- lapply(folders, function(x){message(substr(x, n, nchar(x))); Read10X(data.dir = x) %>% CreateSeuratObject(project = substr(x, n, nchar(x)), min.cells = 3, min.features = 200)})
wu_seurat      <- merge(scrnaseq_files[[1]], scrnaseq_files[-1], add.cell.ids = substr(folders, n, nchar(folders)))
wu_seurat      <- getQC(wu_seurat)

wu_meta <- fread("public_data/Wu/wu_all_metadata.txt") %>% dplyr::rename(barcode = V1) %>% mutate(barcode = gsub(barcode, pattern = "-1", replacement = ""))
wu_meta <- wu_meta[wu_meta$barcode %in% colnames(wu_seurat), ]
wu_meta <- wu_meta[match(colnames(wu_seurat), wu_meta$barcode), ]

wu_meta_seurat           <- cbind(wu_seurat@meta.data, wu_meta)
rownames(wu_meta_seurat) <- wu_meta_seurat$barcode
wu_seurat@meta.data      <- wu_meta_seurat


## Rendeiro (Nat Commun 2020), 3' end solution on blood CD45+ on CLL-patients
rendeiro_seurat <- Read10X(data.dir = "public_data/Rendeiro/") %>% CreateSeuratObject(project = "rendeiro", min.cells = 3, min.features = 200)

rendeiro_seurat$orig.ident <- substr(colnames(rendeiro_seurat), 28, nchar(colnames(rendeiro_seurat)))
rendeiro_seurat$timepoint  <- substr(rendeiro_seurat$orig.ident, 6, nchar(rendeiro_seurat$orig.ident))
rendeiro_seurat$patient    <- substr(rendeiro_seurat$orig.ident, 1, 4)

## Take only d0 samples
Idents(rendeiro_seurat) <- rendeiro_seurat$timepoint
rendeiro_d0_seurat      <- subset(rendeiro_seurat, idents = "d0")
Idents(rendeiro_d0_seurat) <- rendeiro_d0_seurat$orig.ident
rendeiro_d0_seurat <- getQC(rendeiro_d0_seurat)


## 10X human cells
healthy_seurat <- Read10X(data.dir = "public_data/10X/") %>% CreateSeuratObject(project = "healthy_10x", min.cells = 3, min.features = 200)
healthy_seurat <- getQC(healthy_seurat)



## === Merge these data sets together
comparison_seurat     <- merge(wu_seurat, rendeiro_d0_seurat)
comparison_seurat     <- merge(comparison_seurat, healthy_seurat)
comparison_cml_seurat <- merge(cml_seurat, comparison_seurat)

comparison_cml_seurat$project <- "none"
comparison_cml_seurat$project <- ifelse(grepl(comparison_cml_seurat$orig.ident, pattern = "^CLL"), "CLL", comparison_cml_seurat$project)
comparison_cml_seurat$project <- ifelse(grepl(comparison_cml_seurat$orig.ident, pattern = "^Batch"), "CML on TKI", comparison_cml_seurat$project)
comparison_cml_seurat$project <- ifelse(grepl(comparison_cml_seurat$orig.ident, pattern = "^7"), "CML dg", comparison_cml_seurat$project)
comparison_cml_seurat$project <- ifelse(grepl(comparison_cml_seurat$orig.ident, pattern = "^LB"), "Lung", comparison_cml_seurat$project)
comparison_cml_seurat$project <- ifelse(grepl(comparison_cml_seurat$orig.ident, pattern = "^RB"), "Renal", comparison_cml_seurat$project)
comparison_cml_seurat$project <- ifelse(grepl(comparison_cml_seurat$orig.ident, pattern = "^healthy"), "Healthy", comparison_cml_seurat$project)
comparison_cml_seurat$project <- factor(as.character(comparison_cml_seurat$project), levels = c("CML on TKI", "CML dg", "CLL", "Lung", "Renal", "Healthy"))

saveRDS(cml_seurat, "results/cml_seurat_raw.rds")
saveRDS(rendeiro_d0_seurat, "results/rendeiro_d0_seurat.rds")
saveRDS(wu_seurat, "results/wu_seurat.rds")
saveRDS(healthy_seurat, "results/healthy_seurat.rds")
saveRDS(comparison_cml_seurat, "results/comparison_cml_seurat.rds")