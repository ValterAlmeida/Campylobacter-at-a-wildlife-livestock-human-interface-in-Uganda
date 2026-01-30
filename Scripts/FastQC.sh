#FastQC - Quality Control

#!/bin/bash -e
#SBATCH --account       massey03212
#SBATCH --job-name      FastQC
#SBATCH --time          4:00:00
#SBATCH --mem           16GB
#SBATCH --cpus-per-task 16
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out

# Set up working directories
cd cd /nesi/nobackup/massey03212/

# Load Fastqc
module purge
module load FastQC/0.11.9

# Run FastQC

fastqc *.gz -o /nesi/nobackup/massey03212/fastqc_output
