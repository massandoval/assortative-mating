#PBS -l walltime=71:59:00
#PBS -l select=1:ncpus=32:mem=124gb
#$ python3 gnomix.py <query_file> <output_folder> <chr_nr> <phase> <genetic_map_file> <reference_file> <sample_map_file>
module load anaconda3/personal
path_genetic_map="/rds/general/project/human-popgen-datasets/live/HGDP_1000g/gnomix"
path_gnomix="/rds/general/project/human-popgen-datasets/live/HGDP_1000g/gnomix/gnomix"
path_input="/rds/general/project/human-popgen-datasets/live/HGDP_1000g/X"
path_output="/rds/general/project/human-popgen-datasets/ephemeral/output_0.1/X"
python3	$path_gnomix/gnomix.py $path_input/HGDP_1000gNYC_3pop_target_rename_rfmix154_X.phased.vcf.gz $path_output X True $path_genetic_map/chrX_gnomix.b38.gmap $path_input/HGDP_1000gNYC_3pop_ref_rename_rfmix154_X.phased.vcf.gz $path_genetic_map/HGDP_1000g_X.samples
