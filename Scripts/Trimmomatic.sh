#Trimmomatic - Trimming

#!/bin/bash -e
#SBATCH --account       massey03212
#SBATCH --job-name      Trimmomatic
#SBATCH --time          24:00:00
#SBATCH --mem           32GB
#SBATCH --cpus-per-task 16
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out

export SLURM_EXPORT_ENV=ALL

# Set up working directories
cd /nesi/nobackup/massey03212/

# Load Trimmomatic
module purge
module load Trimmomatic/0.39-Java-1.8.0_144

for file in *_1.fq.gz
  do
    base=$(basename ${file} _1.fq.gz)
    trimmomatic PE -threads 16 -phred33 ${file} ${base}_2.fq.gz \
                 ${base}_1.trim.fq.gz ${base}_1un.trim.fq.gz \
                 ${base}_2.trim.fq.gz ${base}_2un.trim.fq.gz \
                 ILLUMINACLIP:/scale_wlg_persistent/filesets/opt_nesi/CS400_centos7_bdw/nullarbor/2.0.20191013/conf/trimmomatic.fa:1:30:11 LEADING:10 TRAILING:10 MINLEN:80 HEADCROP:10
done
