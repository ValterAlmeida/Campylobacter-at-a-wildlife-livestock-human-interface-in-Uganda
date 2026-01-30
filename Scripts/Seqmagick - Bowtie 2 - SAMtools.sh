#Binning
__________________________________________________________________________________________

Seqmagick

module load seqmagick/0.7.0-gimkl-2018b-Python-3.7.3

seqmagick convert --min-length 1000 final.contigs.fa cattle.m1000.fna

seqmagick convert --min-length 1000 contigs.fasta gw2.m1000.fna

seqmagick convert --min-length 1000 contig.fa ug405.m1000.fna
__________________________________________________________________________________________

Bowtie2 - index

#!/bin/bash -e
#SBATCH --account       massey03212
#SBATCH --job-name      BT2_index_CA
#SBATCH --time          00:50:00
#SBATCH --mem           16GB
#SBATCH --cpus-per-task 16
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out

module load Bowtie2/2.4.5-GCC-11.3.0

bowtie2-build --threads $SLURM_CPUS_PER_TASK cattle.m1000.fna cattle_bt2
__________________________________________________________________________________________

Bowtie2 - mapping

#!/bin/bash -e
#SBATCH --account       massey03212
#SBATCH --job-name      bt2_cattle
#SBATCH --time          08:00:00
#SBATCH --mem           08GB
#SBATCH --cpus-per-task 16
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out

module load Bowtie2/2.4.5-GCC-11.3.0

bowtie2 --minins 200 --maxins 800 --threads $SLURM_CPUS_PER_TASK --sensitive -x cattle_bt2 -1 cattle_R1_hostFilt.fastq.gz -2 cattle_R2_hostFilt.fastq.gz -S cattle.bam

__________________________________________________________________________________________

Samtools

#!/bin/bash -e
#SBATCH --account       massey03212
#SBATCH --job-name      Samtools_GW
#SBATCH --time          02:00:00
#SBATCH --mem           24GB
#SBATCH --cpus-per-task 16
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out

module purge
module load SAMtools/1.15.1-GCC-11.3.0

samtools sort -@ $SLURM_CPUS_PER_TASK gw.bam -o gw.sorted.bam

samtools index -@ $SLURM_CPUS_PER_TASK gw.sorted.bam

__________________________________________________________________________________________