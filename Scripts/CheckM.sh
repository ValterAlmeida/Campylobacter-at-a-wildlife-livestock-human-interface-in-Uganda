# CheckM - Genome completeness and contamination

#!/bin/bash -e
#SBATCH --account       massey03212
#SBATCH --job-name      CheckM_GW
#SBATCH --time          05:00:00
#SBATCH --mem           64GB
#SBATCH --cpus-per-task 16
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out

module purge
module load CheckM/1.0.13-gimkl-2018b-Python-2.7.16
export PYTHONNOUSERSITE=1

checkm lineage_wf -t $SLURM_CPUS_PER_TASK --pplacer_threads $SLURM_CPUS_PER_TASK -x fa \
                  --tab_table -f checkm.txt \
                  dastool_out/_DASTool_bins/ checkm_out/
