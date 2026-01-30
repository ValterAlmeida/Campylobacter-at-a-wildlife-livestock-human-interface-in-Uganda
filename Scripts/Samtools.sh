#Samtools

#!/bin/bash -e
#SBATCH --account       massey03212
#SBATCH --job-name      Samtools_HH
#SBATCH --time          04:00:00
#SBATCH --mem           24GB
#SBATCH --cpus-per-task 16
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out

module purge
module load SAMtools/1.15.1-GCC-11.3.0

samtools sort -@ $SLURM_CPUS_PER_TASK hh.bam -o hh.sorted.bam

samtools index -@ $SLURM_CPUS_PER_TASK hh.sorted.bam
