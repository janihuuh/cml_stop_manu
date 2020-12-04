
## Run SlingShot on NK cells
require(slingshot)
require(SingleCellExperiment)
require(SummarizedExperiment)

nk_seurat <- readRDS("results/nk_seurat.rds")
nk_seurat$orig.clusters <- Idents(nk_seurat)
nk_sce <- as.SingleCellExperiment(nk_seurat)

nk_sling <- runSlingshot(nk_sce, cluster_with_earliset_timepoint = "5 CD56 bright", reducedDim = "LATENT_UMAP")
nk_sling_curve <- SlingshotDataSet(nk_sling)

### Trajectories in different diseases
cml_dg_sling      <- getSlingNk(project_temp = "CML dg")
cml_on_tki_sling  <- getSlingNk(project_temp = "CML on TKI")
cml_off_tki_sling <- getSlingNk(project_temp = "CML off TKI")
nsclc_sling       <- getSlingNk(project_temp = "NSCLC")
rcc_sling         <- getSlingNk(project_temp = "RCC")
cll_sling         <- getSlingNk(project_temp = "CLL")
hc_sling          <- getSlingNk(project_temp = "Healthy")

cml_dg_curve      <- SlingshotDataSet(cml_dg_sling)
cml_on_tki_curve  <- SlingshotDataSet(cml_on_tki_sling)
cml_off_tki_curve <- SlingshotDataSet(cml_off_tki_sling)
nsclc_curve       <- SlingshotDataSet(nsclc_sling)
rcc_curve         <- SlingshotDataSet(rcc_sling)
cll_curve         <- SlingshotDataSet(cll_sling)
hc_curve          <- SlingshotDataSet(hc_sling)
