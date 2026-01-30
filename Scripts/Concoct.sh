#CONCOCT - Binning

#!/bin/bash -e
#SBATCH --account       massey03212
#SBATCH --job-name      Concoct_hh
#SBATCH --time          168:00:00
#SBATCH --mem           08GB
#SBATCH --cpus-per-task 08
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out

module purge
module load CONCOCT/1.1.0-gimkl-2020a-Python-3.8.2

cut_up_fasta.py hh.m1000.fna -c 10000 -o 0 --merge_last -b hh_10K.bed > hh_10K.fa

concoct_coverage_table.py hh_10K.bed hh.sorted.bam > coverage_table_hh.tsv

concoct --composition_file hh_10K.fa --coverage_file coverage_table_hh.tsv -b concoct_hh_output/

merge_cutup_clustering.py concoct_hh_output/clustering_gt1000.csv > concoct_hh_output/clustering_merged.csv

mkdir concoct_hh_output/fasta_bins
extract_fasta_bins.py hh.m1000.fna concoct_hh_output/clustering_merged.csv --output_path concoct_hh_output/fasta_bins
