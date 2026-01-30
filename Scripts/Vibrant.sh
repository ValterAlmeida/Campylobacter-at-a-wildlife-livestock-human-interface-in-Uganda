# Vibrant - Identifying viral contigs

#!/bin/bash -e
#SBATCH -A massey03212
#SBATCH -J vibrant
#SBATCH --time 00:40:00
#SBATCH --mem=10GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH -e vibrant.err
#SBATCH -o vibrant.out
#SBATCH --profile=task
#SBATCH --export NONE

export SLURM_EXPORT_ENV=ALL

# Load VIBRANT module
module purge
module load VIBRANT/1.2.1-gimkl-2020a

# Set up working directories
cd /nesi/nobackup/massey03212/Uganda/RNA/human_unwell/megahit_out
mkdir -p vibrant

# Run VIBRANT
srun VIBRANT_run.py -t 16 \
-i rna_human_unwell.m1000.fna \
-d $DB_PATH \
-folder vibrant/
