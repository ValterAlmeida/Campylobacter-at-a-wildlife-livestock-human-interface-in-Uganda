# Concoct - Metabat - Maxbin - DasTool

#Binning

#!/bin/bash -e
#SBATCH -A massey03212
#SBATCH -J Conct_go
#SBATCH --time 180:00:00
#SBATCH --mem 32GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 32
#SBATCH -e concoct.err
#SBATCH -o concoct.out
#SBATCH --export NONE

module load CONCOCT/1.1.0-gimkl-2020a-Python-3.8.2

# Run Concoct

cut_up_fasta.py GO.m1000.fna -c 10000 -o 0 --merge_last -b go_10K.bed > go_10K.fa

concoct_coverage_table.py go_10K.bed GO.sorted.bam > coverage_table_go.tsv

concoct --composition_file go_10K.fa --coverage_file coverage_table_go.tsv -b concoct_go_output/

merge_cutup_clustering.py concoct_go_output/clustering_gt1000.csv > concoct_go_output/clustering_merged.csv

mkdir concoct_go_output/fasta_bins
extract_fasta_bins.py GO.m1000.fna concoct_go_output/clustering_merged.csv --output_path concoct_go_output/fasta_bins

#New:

#!/bin/bash -e
#SBATCH --account       massey03212
#SBATCH --job-name      Concoct_gw
#SBATCH --time          168:00:00
#SBATCH --mem           08GB
#SBATCH --cpus-per-task 08
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out

module purge
module load CONCOCT/1.1.0-gimkl-2020a-Python-3.8.2

cut_up_fasta.py gw1.m1000.fna -c 10000 -o 0 --merge_last -b gw1_10K.bed > gw1_10K.fa

concoct_coverage_table.py gw1_10K.bed gw1.sorted.bam > coverage_table_gw1.tsv

concoct --composition_file gw1_10K.fa --coverage_file coverage_table_gw1.tsv -b concoct_gw1_output/

merge_cutup_clustering.py concoct_gw1_output/clustering_gt1000.csv > concoct_gw1_output/clustering_merged.csv

mkdir concoct_gw1_output/fasta_bins
extract_fasta_bins.py gw1.m1000.fna concoct_gw1_output/clustering_merged.csv --output_path concoct_gw1_output/fasta_bins

__________________________________________________________________________________________

#!/bin/bash -e
#SBATCH -A massey03212
#SBATCH -J Metabat_CA
#SBATCH --time 02:00:00
#SBATCH --mem 32GB
#SBATCH --cpus-per-task 20
#SBATCH -e metabat.err
#SBATCH -o metabat.out
#SBATCH --export NONE

module purge
module load MetaBAT/2.15-GCC-11.3.0

jgi_summarize_bam_contig_depths --outpqutDepth metabat.txt cattle.sorted.bam

metabat2 -t $SLURM_CPUS_PER_TASK -m 1500 \
         -i cattle.m1000.fna \
         -a metabat.txt \
         -o metabat
         
#If it doesnt work: runMetaBat.sh -t $SLURM_CPUS_PER_TASK -m 1500 CA.m1000.fna cattle.sorted.bam

#!/bin/bash -e
#SBATCH --account       massey03212
#SBATCH --job-name      Metabat_hh
#SBATCH --time          01:30:00
#SBATCH --mem           12GB
#SBATCH --cpus-per-task 16
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out

module purge
module load MetaBAT/2.15-GCC-11.3.0

runMetaBat.sh -t $SLURM_CPUS_PER_TASK -m 1500 CA.m1000.fna cattle.sorted.bam

__________________________________________________________________________________________

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
              
__________________________________________________________________________________________

for bin_path in metabat/*.fa;
do
    bin_name=$(basename ${bin_path} .fa)

    grep ">" ${bin_path} | sed 's/>//g' | sed "s/$/\t${bin_name}/g" >> metabat_associations.txt
done


for bin_path in maxbin/*.fasta;
do
    bin_name=$(basename ${bin_path} .fasta)
    grep ">" ${bin_path} | sed 's/>//g' | sed "s/$/\t${bin_name}/g" >> maxbin_associations.txt
done


for bin_path in concoct_cattle_output/fasta_bins/*.fa;
do 
    bin_name=$(basename ${bin_path} .fa)
    grep ">" ${bin_path} | sed 's/>//g' | sed "s/$/\t${bin_name}/g" >> concoct_associations.txt 
done

__________________________________________________________________________________________

#Dereplication
__________________________________________________________________________________________

#DasTool

#!/bin/bash -e
#SBATCH --account       massey03212
#SBATCH --job-name      Dastool_gw
#SBATCH --time          12:00:00
#SBATCH --mem           16GB
#SBATCH --cpus-per-task 16
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out

module purge
mkdir -p dastool_out/

module load DAS_Tool/1.1.5-gimkl-2022a-R-4.2.1
module load DIAMOND/2.0.15-GCC-11.3.0
module load USEARCH/11.0.667-i86linux32

DAS_Tool -i metabat_associations.txt,maxbin_associations.txt -l MetaBAT,MaxBin \
         -t $SLURM_CPUS_PER_TASK --write_bins 1 --search_engine blast \
         -c cattle.m1000.fna \
         -o dastool_out/


__________________________________________________________________________________________

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