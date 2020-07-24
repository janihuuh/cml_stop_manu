
## Activate virtual env

## Activate virtual env
me=$(whoami)
# source /Users/$me/Dropbox/lag3/applications/cellphone_db_venv/bin/activate
source /Users/janihuuh/Dropbox/AML_TIM3/applications/cpdb-venv/bin/activate
cd /Users/$me/Dropbox/cml_stop/results/cellphonedb/

## Use cellphonedb

######## no relapse 0m
cellphonedb method statistical_analysis --iterations=100 --threads=31 \
  --counts-data hgnc_symbol \
  --project-name no_0m \
  input_files/no_0m_meta.txt input_files/no_0m_counts.txt

######## no relapse 6m
cellphonedb method statistical_analysis --iterations=100 --threads=31 \
  --counts-data hgnc_symbol \
  --project-name no_6m \
  input_files/no_6m_meta.txt input_files/no_6m_counts.txt

######## no relapse 12m
cellphonedb method statistical_analysis --iterations=100 --threads=31 \
  --counts-data hgnc_symbol \
  --project-name no_12m \
  input_files/no_12m_meta.txt input_files/no_12m_counts.txt

######## slow relapse 0m
cellphonedb method statistical_analysis --iterations=100 --threads=31 \
  --counts-data hgnc_symbol \
  --project-name slow_0m \
  input_files/slow_0m_meta.txt input_files/slow_0m_counts.txt

######## no relapse 6m
cellphonedb method statistical_analysis --iterations=100 --threads=31 \
  --counts-data hgnc_symbol \
  --project-name slow_6m \
  input_files/slow_6m_meta.txt input_files/slow_6m_counts.txt

######## no relapse 12m
cellphonedb method statistical_analysis --iterations=100 --threads=31 \
  --counts-data hgnc_symbol \
  --project-name slow_relapse \
  input_files/slow_relapse_meta.txt input_files/slow_relapse_counts.txt

######## fast relapse 0m
cellphonedb method statistical_analysis --iterations=100 --threads=31 \
  --counts-data hgnc_symbol \
  --project-name fast_0m \
  input_files/fast_0m_meta.txt input_files/fast_0m_counts.txt

######## fast relapse relapse
cellphonedb method statistical_analysis --iterations=100 --threads=25 \
  --counts-data hgnc_symbol \
  --project-name fast_relapse \
  input_files/fast_relapse_meta.txt input_files/fast_relapse_counts.txt

######## disease_ctrl_0m
cellphonedb method statistical_analysis --iterations=100 --threads=25 \
  --counts-data hgnc_symbol \
  --project-name disease_ctrl_0m \
  input_files/disease_ctrl_0m_meta.txt input_files/disease_ctrl_0m_counts.txt

  ######## tki
  cellphonedb method statistical_analysis --iterations=100 --threads=25 \
    --counts-data hgnc_symbol \
    --project-name cml_tki_0m \
    input_files/cml_tki_0m_meta.txt input_files/cml_tki_0m_counts.txt

    ######## dg

    cellphonedb method statistical_analysis --iterations=100 --threads=25 \
      --counts-data hgnc_symbol \
      --project-name cml_dg_0m \
      input_files/cml_dg_0m_meta.txt input_files/cml_dg_0m_counts.txt
