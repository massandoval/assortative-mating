#PBS -l walltime=23:59:00
#PBS -l select=1:ncpus=1:mem=16gb
#PBS -J 1-1000
module load anaconda3/personal
python run_AutX_pulses_false_fully_connected_scale10E3_rfmix_gnomix_maxanc_onlyfemales.py $PBS_ARRAY_INDEX 1 0 21
python run_AutX_pulses_false_fully_connected_scale10E3_rfmix_gnomix_maxanc_onlyfemales.py $PBS_ARRAY_INDEX 1 1 21
python run_AutX_pulses_false_fully_connected_scale10E3_rfmix_gnomix_maxanc_onlyfemales.py $PBS_ARRAY_INDEX 1 2 21
python run_AutX_pulses_false_fully_connected_scale10E3_rfmix_gnomix_maxanc_onlyfemales.py $PBS_ARRAY_INDEX 1 0 22
python run_AutX_pulses_false_fully_connected_scale10E3_rfmix_gnomix_maxanc_onlyfemales.py $PBS_ARRAY_INDEX 1 1 22
python run_AutX_pulses_false_fully_connected_scale10E3_rfmix_gnomix_maxanc_onlyfemales.py $PBS_ARRAY_INDEX 1 2 22
