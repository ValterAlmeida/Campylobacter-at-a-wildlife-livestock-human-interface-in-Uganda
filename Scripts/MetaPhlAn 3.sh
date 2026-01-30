MetaPhlAn 3

#!/bin/bash -e
#SBATCH -A massey03212
#SBATCH -J HH_mpa3
#SBATCH --time 07:50:00
#SBATCH --mem 30GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH -e metaphlan3.err
#SBATCH -o metaphlan3.out
#SBATCH --export NONE

export SLURM_EXPORT_ENV=ALL

# Set up working directories
cd /nesi/nobackup/massey03212/Uganda/full_co_assembly/Human/Healthy/host_filtered_reads/metaphlan3

# Load Metaphlan3 and Python
module purge
module load MetaPhlAn/3.0.9-gimkl-2020a-Python-3.8.2

mkdir output

for file in *.fastq;
do
base=$(basename ${file} .fastq) 
# Run Metaphlan3

metaphlan ${base}.fastq --input_type fastq > ${base}.txt

metaphlan ${base}.fastq --bowtie2out ${base}.bowtie2.bz2 --nproc 16 --input_type fastq > profiled_${base}.txt

metaphlan ${base}.bowtie2.bz2 --nproc 16 --input_type bowtie2out > metaphlan_${base}.txt

merge_metaphlan_tables.py *profiled* > output/merged_abundance_table.txt

done