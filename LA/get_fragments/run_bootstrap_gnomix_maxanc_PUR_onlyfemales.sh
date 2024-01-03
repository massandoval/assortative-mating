#PBS -l walltime=7:59:00
#PBS -l select=1:ncpus=1:mem=16gb
module load anaconda3/personal
python bootstrap_gnomix_maxanc_onlyfemales.py PUR
