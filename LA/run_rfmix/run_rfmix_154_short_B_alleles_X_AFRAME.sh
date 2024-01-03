#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=32:mem=62gb
module load anaconda3/personal
python2 RunRFMix.py PopPhased /rds/general/project/human-popgen-datasets/live/HGDP_1000g/X/HGDP_1000gNYC_3pop_target_ref_rfmix154_short_B_X_AFRAME.phased.alleles /rds/general/project/human-popgen-datasets/live/HGDP_1000g/RFMix/classes_rfmix1.5.4_AFRAME.txt /rds/general/project/human-popgen-datasets/live/HGDP_1000g/RFMix/markerslocation_short_B_chrX_AFRAME.txt -o /rds/general/project/human-popgen-datasets/live/HGDP_1000g/RFMix/HGDP_1000gNYC_3pop_target_ref_output_rfmix154_short_B_X_AFRAME -w 0.1 -G 19 --num-threads 32 -e 3 --forward-backward

