
## Activate virtual env
me=$(whoami)
source /Users/$me/Dropbox/AML_TIM3/applications/cpdb-venv/bin/activate
cd /Users/$me/Dropbox/cml_stop/results/cellphonedb/

## Use cellphonedb

######## no relapse 0m
cellphonedb method statistical_analysis --iterations=1000 --threads=31 \
  --counts-data hgnc_symbol \
  --project-name public_total \
  input_files/public_total_meta.txt input_files/public_total_counts.txt
