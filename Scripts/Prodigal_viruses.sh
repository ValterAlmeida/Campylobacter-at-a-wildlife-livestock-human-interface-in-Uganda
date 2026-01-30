# Prodigal - Viruses

Prodigal

#!/bin/bash -e
#SBATCH -A massey03212
#SBATCH -J prodigal
#SBATCH --time 00:40:00
#SBATCH --mem=10GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH -e prodigal.err
#SBATCH -o prodigal.out
#SBATCH --export NONE

export SLURM_EXPORT_ENV=ALL

# Load dependencies
module purge
module load Prodigal/2.6.3-GCC-9.2.0

# Set up working directories
cd /nesi/nobackup/massey03212/Uganda/full_co_assembly/Human/Healthy/host_filtered_reads/megahit_out/vibrant

mkdir -p viral_taxonomy

# Run main analyses
srun prodigal -p meta -q \
-i checkv_combined.fna \
-a viral_taxonomy/checkv_combined.faa
