#!/bin/bash -e

#SBATCH --account       massey03212
#SBATCH --job-name      annotate_pfam
#SBATCH --time          01:00:00
#SBATCH --mem           8GB
#SBATCH --cpus-per-task 10
#SBATCH --array         0-121
#SBATCH --error         %x_%A_%a.err
#SBATCH --output        %x_%A_%a.out

# Load modules
module purge
module load HMMER/3.4-GCC-12.3.0

# Working directory
cd /nesi/nobackup/massey03212/Nat_comms/Genomic_analysis/10_Metabolic_functions/04_gene_annotation_and_coverage

# Get file for this task
prot_file=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" faa_list.txt)
out_file=$(basename "$prot_file" .faa)

# HMM database
db=Pfam-A.hmm

# Run HMMER
hmmsearch --tblout gene_annotations/${out_file}.pfam.txt -E 0.001 \
          --cpu $SLURM_CPUS_PER_TASK \
          $db $prot_file
