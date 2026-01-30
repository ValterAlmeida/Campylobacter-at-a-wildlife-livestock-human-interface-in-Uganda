#Prokka

#!/bin/bash -e
#SBATCH --account       massey03212
#SBATCH --job-name      Prokka_Trep
#SBATCH --time          48:00:00
#SBATCH --mem           500MB
#SBATCH --cpus-per-task 04
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out

module purge
module load prokka/1.14.5-GCC-9.2.0

for file in *.fa
do
    base=$(basename ${file} .fa)

    # Run PROKKA with --prefix and -o
    prokka ${file} --prefix ${base} -o ${base} --cpus $SLURM_CPUS_PER_TASK
done
