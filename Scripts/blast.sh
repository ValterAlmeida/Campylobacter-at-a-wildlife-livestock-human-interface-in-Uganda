#!/bin/bash -e
#SBATCH --account massey03212
#SBATCH --job-name Blast_Campy
#SBATCH --time 6:00:00
#SBATCH --mem 4GB
#SBATCH --cpus-per-task 6
#SBATCH --error %x_%j.err
#SBATCH --output %x_%j.out

# Load necessary modules
module purge
module load BLASTDB/2025-08
module load BLAST/2.16.0-GCC-12.3.0

blastn \
  -query blast_contigs/all_samples.contigs.fa \
  -db nt \
  -out blast_contigs/contigs_vs_refseq.blastn.tsv \
  -outfmt "6 qseqid sseqid pident length qcovs evalue bitscore staxids sscinames" \
  -evalue 1e-10 \
  -max_target_seqs 5 \
  -num_threads 6
