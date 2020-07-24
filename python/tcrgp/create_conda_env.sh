
path_to_venv=/scratch/cs/csb/projects/aplastic_anemia/tcrgp_env
module load anaconda3

# create environment with package pip in it
module load teflon
conda create --prefix $path_to_venv python pip ipython
module unload teflon

# this must be run in each shell to set up the environment variables properly.
source activate $path_to_venv

# install more packages, either conda or pip
pip install tensorflow==2.0.0-beta1 
conda install tensorflow

# pip install tensorflow
pip install gpflow
pip install -U scikit-learn

# install as a kernel
pip install ipykernel
python -m ipykernel install --user --name=python-tcrgp --display-name="tcrgp-ENV"

## Finish
module unload teflon
conda deactivate
# source deactivate $path_to_venv
