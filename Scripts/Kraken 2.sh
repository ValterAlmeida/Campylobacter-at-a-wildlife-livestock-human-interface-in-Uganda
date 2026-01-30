# Kraken 2

#!/bin/bash -e
#SBATCH --account       massey03212
#SBATCH --job-name      Krak_Brac
#SBATCH --time          10:00:00
#SBATCH --mem           08GB
#SBATCH --cpus-per-task 08
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out

#module load
module purge
module load Kraken2/2.1.2-GCC-9.2.0
module load Bracken/2.7-GCC-11.3.0

#First step

for file in *.fastq
do
base=$(basename ${file} .fastq)
kraken2 --db /opt/nesi/db/Kraken2/standard-2022-07/ --threads $SLURM_CPUS_PER_TASK --use-names ${base}.fastq --output ${base}.txt --report ${base}.kreport
done

#Second step

for file in *.kreport
do
base=$(basename ${file} .kreport)
bracken -d /opt/nesi/db/Kraken2/standard-2022-07/ -i ${base}.kreport -l S -o ${base}_bracken.report -t 10 -r 100
done

# Combine reports:

# For kraken2 combine - Python/2.7.18-gimkl-2020a

python combine_bracken_outputs.py --files *_bracken_species_bracken.report -o bracken_all.report.tsv
