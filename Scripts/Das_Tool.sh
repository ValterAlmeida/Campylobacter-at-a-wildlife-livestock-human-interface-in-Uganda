#DAS_Tool - Dereplication

#!/bin/bash -e
#SBATCH --account       massey03212
#SBATCH --job-name      Dastool_gh
#SBATCH --time          12:00:00
#SBATCH --mem           08GB
#SBATCH --cpus-per-task 18
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out

module purge
mkdir -p dastool_out/

module load DAS_Tool/1.1.1-gimkl-2018b-R-3.6.1
module load DIAMOND/0.9.25-gimkl-2018b
module load USEARCH/11.0.667-i86linux32

DAS_Tool -i metabat_associations.txt,maxbin_associations.txt,concoct_associations.txt -l MetaBAT,MaxBin,Concoct \
         -t $SLURM_CPUS_PER_TASK --write_bins 1 --search_engine blast \
         -c gh.m1000.fna \
         -o dastool_out/
