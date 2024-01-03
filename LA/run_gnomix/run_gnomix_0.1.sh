#PBS -l walltime=7:59:00
#PBS -l select=1:ncpus=32:mem=256gb
PBS_ARRAY_INDEX=12
#$ python3 gnomix.py <query_file> <output_folder> <chr_nr> <phase> <genetic_map_file> <reference_file> <sample_map_file>
module load anaconda3/personal
path_genetic_map="/rds/general/project/human-popgen-datasets/live/HGDP_1000g/gnomix"
path_gnomix="/rds/general/project/human-popgen-datasets/live/HGDP_1000g/gnomix/gnomix"
path_input="/rds/general/project/human-popgen-datasets/live/HGDP_1000g/${PBS_ARRAY_INDEX}"
path_output="/rds/general/project/human-popgen-datasets/ephemeral/output_0.1/${PBS_ARRAY_INDEX}"
python3	$path_gnomix/gnomix.py $path_input/HGDP_1000gNYC_3pop_target_rename_rfmix154_${PBS_ARRAY_INDEX}.phased.vcf.gz $path_output $PBS_ARRAY_INDEX True $path_genetic_map/chr${PBS_ARRAY_INDEX}_gnomix.b38.gmap $path_input/HGDP_1000gNYC_3pop_ref_rename_rfmix154_${PBS_ARRAY_INDEX}.phased.vcf.gz $path_genetic_map/HGDP_1000g.samples