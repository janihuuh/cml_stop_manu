
## Run cellphone db. Make subselection of 1000 cells from responders at time 0 and time 1
sample_n = 1000

################################ Init cellphone input files

### no relapse
data.frame(cml_seurat@meta.data, cluster = Idents(cml_seurat)) %>%
  initCellphonedb(seurat_object = public_seurat, name = "public_seurat", sample_n = sample_n)
