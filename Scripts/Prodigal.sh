# Prodigal - Annotation

#!/bin/bash
#SBATCH -A massey03212
#SBATCH -J prodigal
#SBATCH --time 2-00:00:00
#SBATCH --mem=200MB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --requeue
#SBATCH -e slurm-%x.%A-%a.err
#SBATCH -o slurm-%x.%A-%a.out

module purge
module load Prodigal/2.6.3-GCC-9.2.0
export EXTRACT_DIR=/nesi/nobackup/massey03212/Treponema_chapter/02_pooled_non_dereplicated_genomes
export PROD_DIR=/nesi/nobackup/massey03212/Treponema_chapter/07_annotation/01_Prodigal
