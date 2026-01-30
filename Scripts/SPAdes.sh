# SPAdes - Test


#!/bin/bash -e
#SBATCH -A massey03212
#SBATCH -J SPAdes_gorwild
#SBATCH --time 08:00:00
#SBATCH --mem 27GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH -e sapdes_gorillawild.err
#SBATCH -o spades_gorillawild.out
#SBATCH --export NONE

export SLURM_EXPORT_ENV=ALL

# Set up working directories
cd /nesi/project/massey03212/AGRF/AGRF_CAGRF21098330_H3CGTDSX3/Trimmomatic/bbmap/Gorilla_wild/host_filtered_reads/spades

# Load SPAdes module
module purge
module load SPAdes/3.15.3-gimkl-2020a

for file in *_R1_hostFilt.fastq.gz
  do
    base=$(basename ${file} _R1_hostFilt.fastq.gz)

# Run SPAdes

spades.py --meta -k 33,55,77,99,121 -t 16 -1 ${file} -2 ${base}_R2_hostFilt.fastq.gz -o spades_assembly/

done
