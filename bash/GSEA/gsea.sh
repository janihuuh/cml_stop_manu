#!/bin/bash
me=whoami

## pDC
cd /Users/$me/Dropbox/cml_stop

RNKFILE=results/gsea/pDC/pDCs.rnk
OUTDIR=results/gsea/pDC/
LABEL=pDC_reactome

## 28_pDC time point
cd /Users/$me/Dropbox/cml_stop

RNKFILE=results/gsea/28_pDC/baseline_vs_6m.rnk
OUTDIR=results/gsea/28_pDC/
LABEL=reactome

GMT=/Users/$me/Dropbox/applications/GSEA/h.all.v6.2.symbols.gmt # hallmark
GMT=/Users/$me/Dropbox/applications/GSEA/c2.cp.reactome.v7.0.symbols.gmt # reactome


## nk cd56dim over time
cd /Users/$me/Dropbox/cml_stop

RNKFILE=results/gsea/nk_cd56_dim/baseline_vs_6m.rnk
OUTDIR=results/gsea/nk_cd56_dim/
LABEL=nk_cd56_dim_baseline_vs_6m_reactome

GMT=/Users/$me/Dropbox/applications/GSEA/h.all.v6.2.symbols.gmt # hallmark
GMT=/Users/$me/Dropbox/applications/GSEA/c2.cp.reactome.v7.0.symbols.gmt # reactome


## nk cd56bright over time
cd /Users/$me/Dropbox/cml_stop

RNKFILE=results/gsea/nk_cd56_bright/baseline_vs_6m.rnk
OUTDIR=results/gsea/nk_cd56_bright/
LABEL=reactome

GMT=/Users/$me/Dropbox/applications/GSEA/h.all.v6.2.symbols.gmt # hallmark
GMT=/Users/$me/Dropbox/applications/GSEA/c2.cp.reactome.v7.0.symbols.gmt # reactome


### ============= no need to modify
folder=/Users/$me/Dropbox/cml_stop
gsea_path=/Users/$me/Dropbox/applications/gsea-3.0.jar

cd $folder
mdkir $OUTDIR

java -cp $gsea_path \
    -Xmx5g xtools.gsea.GseaPreranked \
    -rpt_label $LABEL \
    -rnk $folder/$RNKFILE \
    -gmx $GMT \
    -out $folder/$OUTDIR \
    -plot_top_x 250 \
    -collapse false \
    -mode Max_probe \
    -norm meandiv \
    -scoring_scheme weighted \
    -include_only_symbols true \
    -make_sets true \
    -rnd_seed 149 \
    -zip_report false \
    -gui false \
    -nperm 1000 \
    -set_min 5 \
    -set_max 500
