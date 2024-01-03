#PBS -l walltime=7:59:00
#PBS -l select=1:ncpus=1:mem=16gb

module load anaconda3/personal
python get_fragments_and_bootstrap_PEL_gnomix_0.1cM_maxanc_X_onlyfemales.py 
