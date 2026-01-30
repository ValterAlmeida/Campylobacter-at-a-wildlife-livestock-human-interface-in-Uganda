# CheckV

#!/bin/bash -e
#SBATCH -A massey03212
#SBATCH -J checkv
#SBATCH --time 08:40:00
#SBATCH --mem 128GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH -e checkv.err
#SBATCH -o checkv.out
#SBATCH --profile=task
#SBATCH --export NONE

export SLURM_EXPORT_ENV=ALL

# Load the module
module purge
module load CheckV/0.7.0-gimkl-2020a-Python-3.8.2

# Set up working directories
cd /nesi/nobackup/massey03212/Uganda/full_co_assembly/Human/Healthy/host_filtered_reads/megahit_out/vibrant/checkv2
mkdir -p checkv_out

# Run main analyses
checkv_in="vibrant/VIBRANT_healthy_vibrant.m1000/VIBRANT_phages_rna_human_unwell.m1000/rna_human_unwell.m1000.phages_combined.fna"

srun checkv end_to_end ${checkv_in} checkv_out -t 10 --quiet
