# MaxBin - Binning

#!/bin/bash -e
#SBATCH --account       massey03212
#SBATCH --job-name      Maxbin_gw
#SBATCH --time          96:30:00
#SBATCH --mem           27GB
#SBATCH --cpus-per-task 16
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out

module purge
module load MaxBin/2.2.6-gimkl-2018b-Perl-5.28.1

cut -f1,4,6,8,10 metabat.txt > maxbin.txt

mkdir -p maxbin/
run_MaxBin.pl -thread $SLURM_CPUS_PER_TASK -min_contig_length 1500 \
              -contig gw.m1000.fna \
              -abund maxbin.txt \
              -out maxbin
