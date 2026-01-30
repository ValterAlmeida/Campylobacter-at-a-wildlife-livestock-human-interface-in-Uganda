#!/bin/bash -e
#SBATCH --account       massey03212
#SBATCH --job-name      CheckM_Campy
#SBATCH --time          24:00:00
#SBATCH --mem           48GB
#SBATCH --cpus-per-task 12
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out

module purge
module load CheckM2/1.0.1-Miniconda3

checkm2 predict -t $SLURM_CPUS_PER_TASK -x fa \
                --input /nesi/nobackup/massey03212/Nat_comms/Genomic_analysis/10_Metabolic_functions/Mags_refs \
                --output-directory checkm_out/ \
                --stdout > checkm2.txt
