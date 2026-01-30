# Roary

#!/bin/bash -e
#SBATCH --account       massey03212
#SBATCH --job-name      Roary_Prevotella
#SBATCH --time          02:20:00
#SBATCH --mem           1GB
#SBATCH --cpus-per-task 08
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out

module purge

module load Roary/3.13.0-gimkl-2020a

roary -p $SLURM_CPUS_PER_TASK -f roary -e -n -v *.gff
