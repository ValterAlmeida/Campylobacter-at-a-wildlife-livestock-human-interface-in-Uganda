#!/bin/bash -e

#SBATCH --account       massey03212
#SBATCH --job-name      annotate_DRAM
#SBATCH --time          48:00:00
#SBATCH --mem           32GB
#SBATCH --cpus-per-task 24
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out

# Load modules
module purge
module load DRAM/1.3.5-Miniconda3

# Set working directory
cd /nesi/nobackup/massey03212/Nat_comms/Genomic_analysis/10_Metabolic_functions/04_gene_annotation_and_coverage

# Run DRAM
DRAM.py annotate \
  -i 'Mags_refs/*.fa' \
  --checkm_quality DRAM_input_files/checkm2.tsv \
  --gtdb_taxonomy DRAM_input_files/gtdbtk.bac120.summary.tsv \
  -o dram_annotations \
  --threads $SLURM_CPUS_PER_TASK
