#!/bin/bash -e
#SBATCH --account    massey03212
#SBATCH --job-name   Abricate_Campy
#SBATCH --time       4:00:00
#SBATCH --mem        2GB
#SBATCH --cpus-per-task 04
#SBATCH --error      %x_%j.err
#SBATCH --output     %x_%j.out

# Load modules
module purge
module load ABRicate/1.0.0-GCC-11.3.0-Perl-5.34.1

# Set working and output directories
WORKDIR=/nesi/nobackup/massey03212/Nat_comms/Genomic_analysis/06_Abricate/Campylobacter_infans
OUTDIR=/nesi/nobackup/massey03212/Nat_comms/Genomic_analysis/06_Abricate/Campylobacter_infans

# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"

# Move to working directory
cd "$WORKDIR"

# Run ABRicate and save outputs to OUTDIR
abricate --db resfinder *.fa > "$OUTDIR/resfinder.tsv"
abricate --db vfdb *.fa > "$OUTDIR/vfdb.tsv"

