# vContact2

# activate vcontact2 conda environment
module purge
module load Miniconda3
source activate vContact2

# Load dependencies
export PATH="/nesi/nobackup/massey03212/Uganda/full_co_assembly/Human/Healthy/host_filtered_reads/megahit_out/vibrant/vContact2:$PATH"
module load DIAMOND/0.9.32-GCC-9.2.0
module load MCL/14.137-gimkl-2020a

# run vcontact2_gene2genome
vcontact2_gene2genome -p viral_taxonomy/checkv_combined.faa -o viral_taxonomy/viral_genomes_g2g.csv -s 'Prodigal-FAA'

# deactivate conda environment
conda deactivate
