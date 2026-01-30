#Bowtie2 - index

#!/bin/bash -e
#SBATCH --account       massey03212
#SBATCH --job-name      BT2_index_HH
#SBATCH --time          00:50:00
#SBATCH --mem           16GB
#SBATCH --cpus-per-task 16
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out

module load Bowtie2/2.4.5-GCC-11.3.0

bowtie2-build --threads $SLURM_CPUS_PER_TASK hh.m1000.fna hh_bt2
__________________________________________________________________________________________

#Bowtie2 - mapping

#!/bin/bash -e
#SBATCH --account       massey03212
#SBATCH --job-name      bt2_HH
#SBATCH --time          08:00:00
#SBATCH --mem           06GB
#SBATCH --cpus-per-task 16
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out

module load Bowtie2/2.4.5-GCC-11.3.0

bowtie2 --minins 200 --maxins 800 --threads $SLURM_CPUS_PER_TASK --sensitive -x hh_bt2 -1 cattle_R1_hostFilt.fastq.gz -2 cattle_R2_hostFilt.fastq.gz -S cattle.bam
