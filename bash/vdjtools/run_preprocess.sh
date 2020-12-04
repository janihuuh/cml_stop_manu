
## Init objects to vdjtools scripts

#############
me=$(whoami)
application_files=/Users/$me/Dropbox/aplastic_anemia_tcr
files=/Users/$me/Dropbox/cml_stop/

vdj=$application_files/applications/vdjtools-1.2.1/vdjtools-1.2.1.jar
vdjdb=$application_files/applications/vdjdb-1.1.5/vdjdb-1.1.5.jar
vdjdb_new=$application_files/data/selected_tcrb/databases/vdjdb_new


#############

cd $files/tcrb_data/unsorted/cml_sc/
cml_unsort_sc=$(ls -d "$PWD"/*);

cd $files/tcrb_data/unsorted/cml/
cml_unsort_cml=$(ls -d "$PWD"/*);

cd /Users/$me/Dropbox/Emerson/data/vdjt/
emerson=$(ls -d "$PWD"/*);

###############
cd $files
###############
clear


## Preprocess; convert into more human readable format
java -Xmx4G -jar $vdj Convert -S ImmunoSeqV2 $cml_unsort_sc $files/tcrb_data/unsorted/cml_sc/
java -Xmx4G -jar $vdj Convert -S ImmunoSeqV2 $cml_unsort_cml $files/tcrb_data/unsorted/cml/

## Remove unproductive clonotypes
java -Xmx4G -jar $vdj FilterNonFunctional $cml_unsort_sc $files/tcrb_data/unsorted/cml_sc/
java -Xmx4G -jar $vdj FilterNonFunctional $cml_unsort_cml $files/tcrb_data/unsorted/cml/

## Subsample emerson to 40k
java -Xmx4G -jar $vdj DownSample -x	40000 $emerson $files/tcrb_data/unsorted/emerson40k/

## Calc diversity stats for each sample
unsort=$(echo $cml_unsort_sc $cml_unsort_cml)
unsort_healthy=$(echo $cml_unsort_sc $cml_unsort_dasastop $emerson)

java -Xmx4G -jar $vdj CalcDiversityStats --intersect-type aa $unsort $files/tcrb_data/results/diversity/cml_unsampled
java -Xmx4G -jar $vdj CalcDiversityStats --intersect-type aa --downsample-to 10000 $unsort $files/tcrb_data/results/diversity/10k_sampled
java -Xmx4G -jar $vdj CalcDiversityStats --intersect-type aa --downsample-to 10000 $unsort_healthy $files/tcrb_data/results/diversity/10k_sampled_w_healthy
