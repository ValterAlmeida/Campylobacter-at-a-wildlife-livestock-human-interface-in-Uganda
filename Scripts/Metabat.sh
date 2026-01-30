# Metabat - Binning

#!/bin/bash -e
#SBATCH --account       massey03212
#SBATCH --job-name      Metabat_hh
#SBATCH --time          01:30:00
#SBATCH --mem           12GB
#SBATCH --cpus-per-task 16
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out

module purge
module load MetaBAT/2.15-GCC-11.3.0

runMetaBat.sh -t $SLURM_CPUS_PER_TASK -m 1500 hh.m1000.fna human.sorted.bam
