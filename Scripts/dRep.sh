#!/bin/bash -e
#SBATCH --account       massey03212
#SBATCH --job-name      dRep_gtdb_2_4
#SBATCH --time          24:00:00
#SBATCH --mem           32GB
#SBATCH --cpus-per-task 16
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out

module purge

module load CheckM/1.0.13-gimkl-2018b-Python-2.7.16
module load drep/2.3.2-gimkl-2018b-Python-3.7.3
module load Mash/2.1-gimkl-2018b

dRep dereplicate dRep_output/ -p $SLURM_CPUS_PER_TASK --genomeInfo dRep.genomeInfo -g ../dastool_out/_DASTool_bins/*.fa -str 100 -pa 0.95 -sa 0.99
