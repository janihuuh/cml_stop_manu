
## Fetch data from FIMM cluster
## scRNAseq-data
me=$(whoami)

scrna_files=outs/filtered_feature_bc_matrix
molecule_file=outs/molecule_info.h5

bam_file=outs/possorted_genome_bam.bam
bai_file=outs/possorted_genome_bam.bam.bai
barcode_file=outs/filtered_feature_bc_matrix/barcodes.tsv.gz

local_folder=/Users/$me/Dropbox/cml_stop_clean/data/scRNAseq


## CML stop dg; part I/II
scrna_path=/fas/NGS/pipes/cellranger/fimm_sca_mustjoki/NordCML007/Batch1-5_151019-121119/count_191203_A00464_0133_BHMCGHDMXX
names=$(echo 716_dg 720_dg)
names=$(echo 720_dg)

for name in $names; do
  mkdir $local_folder/$name/;
  scp -r jhuuhtan@atlas.genome.helsinki.fi:$scrna_path/$name/$scrna_files $local_folder/$name/
done

## CML stop dg; part II/II
scrna_path=/fas/NGS/pipes/cellranger/fimm_sca_mustjoki/NordCML007/Batch1-5_151019-121119/count_191203_A00464_0134_AHM7FHDMXX
names=$(echo 706_dg 730_dg)

for name in $names; do
  mkdir $local_folder/$name/;
  scp -r jhuuhtan@atlas.genome.helsinki.fi:$scrna_path/$name/$scrna_files $local_folder/$name/
done


## CML stop TKI; part I/VII
scrna_path=/fas/NGS/pipes/cellranger/fimm_sca_mustjoki/CML_STOP_TKI/Batch1_310119/count_190313_A00464_0050_AHJ3G3DMXX/
names=$(echo SI-GA-C8 SI-GA-D8)

for name in $names; do
  mkdir $local_folder/$name/;
  scp -r jhuuhtan@atlas.genome.helsinki.fi:$scrna_path/$name/$scrna_files $local_folder/$name/
done

## CML stop TKI; part II/VII
scrna_path=/fas/NGS/pipes/cellranger/fimm_sca_mustjoki/CML_STOP_TKI/Batch2_210219/count_190313_A00464_0050_AHJ3G3DMXX/
names=$(echo SI-GA-F8)

for name in $names; do
  mkdir $local_folder/$name/;
  scp -r jhuuhtan@atlas.genome.helsinki.fi:$scrna_path/$name/$scrna_files $local_folder/$name/
done

## CML stop TKI; part III/VII
scrna_path=/fas/NGS/pipes/cellranger/fimm_sca_mustjoki/CML_STOP_TKI/Batch2_210219/count_HJ3G3DMXX_HJMGWDMXX/
names=$(echo SI-GA-A9 SI-GA-H8)

for name in $names; do
  mkdir $local_folder/$name/;
  scp -r jhuuhtan@atlas.genome.helsinki.fi:$scrna_path/$name/$scrna_files $local_folder/$name/
done

## CML stop TKI; part IV/VII
scrna_path=/fas/NGS/pipes/cellranger/fimm_sca_mustjoki/CML_STOP_TKI/Batch3_180319/count_190423_A00464_0057_AHJMGWDMXX/
names=$(echo SI-GA-C12 SI-GA-D12 SI-GA-E12)

for name in $names; do
  mkdir $local_folder/$name/;
  scp -r jhuuhtan@atlas.genome.helsinki.fi:$scrna_path/$name/$scrna_files $local_folder/$name/
done

## CML stop TKI; part V/VII
scrna_path=/fas/NGS/pipes/cellranger/fimm_sca_mustjoki/CML_STOP_TKI/Batch4_260319/count_190509_A00464_0063_AHJTCMDMXX/
names=$(echo 407_12months  407_baseline)

for name in $names; do
  mkdir $local_folder/$name/;
  scp -r jhuuhtan@atlas.genome.helsinki.fi:$scrna_path/$name/$scrna_files $local_folder/$name/
done

## CML stop TKI; part VI/VII
scrna_path=/fas/NGS/pipes/cellranger/fimm_sca_mustjoki/CML_STOP_TKI/Batch5_270319/count_190506_A00464_0062_BH7LJVDRXX/
names=$(echo 822_12months_270319  822_6months_270319  822_baseline_270319)

for name in $names; do
  mkdir $local_folder/$name/;
  scp -r jhuuhtan@atlas.genome.helsinki.fi:$scrna_path/$name/$scrna_files $local_folder/$name/
done

## CML stop TKI; part VII/VII
scrna_path=/fas/NGS/pipes/cellranger/fimm_sca_mustjoki/CML_STOP_TKI/Batch5_270319/count_190506_A00464_0062_BH7LJVDRXX/
names=$(echo 1618_baseline  1624_baseline 1618_relapse   1624_relapse)

for name in $names; do
  mkdir $local_folder/$name/;
  scp -r jhuuhtan@atlas.genome.helsinki.fi:$scrna_path/$name/$scrna_files $local_folder/$name/
done






## ============= TCRab data
file1=all_contig_annotations.csv
file2=clonotypes.csv
file3=consensus_annotations.csv
file4=filtered_contig_annotations.csv
file5=metrics_summary.csv
tcr_files=$(echo $file1 $file2 $file3 $file4 $file5)

local_folder=/Users/$me/Dropbox/cml_stop_clean/data/scRNAseq+TCRseq/
mkdir $local_folder


## Batch1
remote_folder=/fas/NGS/pipes/cellranger/fimm_sca_mustjoki/TIM3/Batch1_220219/vdj_190423_A00464_0056_BH7JCVDRXX
names=$(echo FH6753_TCR FH6088_TCR)

for name in $names; do
  echo $name
  mkdir $local_folder/$name/;
  for file in $tcr_files; do
    echo $file
    scp jhuuhtan@atlas.genome.helsinki.fi:$remote_folder/$name/outs/$file $local_folder/$name;
  done;
done
