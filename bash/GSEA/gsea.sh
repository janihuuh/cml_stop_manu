#!/bin/bash
me=whoami

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
