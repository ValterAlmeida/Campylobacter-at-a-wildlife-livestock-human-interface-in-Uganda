#BBMap

# Build index reference - Human

#!/bin/bash -e
#SBATCH -A massey03212
#SBATCH -J host_filt_bbmap_index
#SBATCH -J Human_ref
#SBATCH --time 00:40:00
#SBATCH --mem 23GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH -e host_filt_bbmap_index.err
#SBATCH -o host_filt_bbmap_index.out
#SBATCH --export NONE

export SLURM_EXPORT_ENV=ALL

cd /nesi/project/massey03212/AGRF/AGRF_CAGRF21098330_H3CGTDSX3/Trimmomatic/bbmap/Human_ref

# Load BBMap module
module purge
module load BBMap/38.95-gimkl-2020a

# Build indexed reference file via BBMap
srun bbmap.sh ref=GCF_000001405.39_GRCh38.p13_genomic.fna -Xmx23g

________________________________________________________________________________

# Removing host contamination - Human:

#!/bin/bash -e
#SBATCH -A massey03212
#SBATCH -J bbmap_Hum_Com
#SBATCH --time 14:00:00
#SBATCH --mem 27GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20
#SBATCH -e bbmap_human_Com.err
#SBATCH -o bbmap_human_Com.out
#SBATCH --export NONE

export SLURM_EXPORT_ENV=ALL

# Set up working directories
cd /nesi/project/massey03212/AGRF/AGRF_CAGRF21098330_H3CGTDSX3/Trimmomatic/bbmap/Human_community
mkdir -p host_filtered_reads/

# Load BBMap module
module purge
module load BBMap/38.95-gimkl-2020a

for file in *_R1.trim.fastq.gz
  do
    base=$(basename ${file} _R1.trim.fastq.gz)

# Run bbmap

srun bbmap.sh in1=${file} in2=${base}_R2.trim.fastq.gz -Xmx27g -t=20 \
minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 qtrim=rl trimq=10 untrim \
path=/nesi/project/massey03212/AGRF/AGRF_CAGRF21098330_H3CGTDSX3/Trimmomatic/bbmap/Human_ref \
out1=/nesi/project/massey03212/AGRF/AGRF_CAGRF21098330_H3CGTDSX3/Trimmomatic/bbmap/Human_community/host_filtered_reads/${base}_R1_hostFilt.fastq.gz out2=/nesi/project/massey03212/AGRF/AGRF_CAGRF21098330_H3CGTDSX3/Trimmomatic/bbmap/Human_community/host_filtered_reads/${base}_R2_hostFilt.fastq.gz
done
