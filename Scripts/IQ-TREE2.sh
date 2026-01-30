# IQ-TREE2 - Phylogenetic tree

#!/bin/bash -e
#SBATCH --account       massey03212
#SBATCH --job-name      IQTree_GW
#SBATCH --time          168:00:00
#SBATCH --mem           26GB
#SBATCH --cpus-per-task 48
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out

module purge
module load IQ-TREE/2.2.0.5-gimpi-2022a

iqtree2 -s gtdbtk.bac120.user_msa.fasta.gz -T $SLURM_CPUS_PER_TASK -mem 26G
