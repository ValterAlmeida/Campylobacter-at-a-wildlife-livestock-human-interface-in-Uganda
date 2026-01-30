#!/bin/bash -e
#SBATCH --account       massey03212
#SBATCH --job-name      Megahit
#SBATCH --time          36:00:00
#SBATCH --mem           128GB
#SBATCH --cpus-per-task 16
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out

# Load MEGAHIT module
module purge
module load MEGAHIT/1.1.4-gimkl-2018b-Python-2.7.16

for file in *_R1_hostFilt.fastq.gz
  do
    base=$(basename ${file} _R1_hostFilt.fastq.gz)

# Run MEGAHIT

megahit -1 ${file} -2 ${base}_R2_hostFilt.fastq.gz -o ${base}_megahit_out

done
