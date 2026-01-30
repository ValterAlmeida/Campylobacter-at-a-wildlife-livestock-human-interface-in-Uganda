# IDBA-UD - Test

#!/bin/bash -e
#SBATCH -A massey03212
#SBATCH -J IDBA_Hum_Heal
#SBATCH --time 01:20:00
#SBATCH --mem 20GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 12
#SBATCH -e idba.err
#SBATCH -o idba.out
#SBATCH --export NONE

export SLURM_EXPORT_ENV=ALL

# Set up working directories
cd /nesi/nobackup/massey03212/Uganda/full_co_assembly/Human/Healthy/individual/host_filtered_reads/mega_spades/IDBA

# Load IDBA module
module purge
module load IDBA-UD/1.1.3-gimkl-2018b

for file in *.fastq;
do
base=$(basename ${file} .fastq)
# Run IDBA
idba_ud --num_threads 8 --mink 33 --maxk 99 --step 22 -l ${base}.fastq -o ${base}
done
