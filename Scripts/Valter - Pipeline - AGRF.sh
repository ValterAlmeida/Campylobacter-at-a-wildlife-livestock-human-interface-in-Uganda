#Valter - Pipeline - AGRF files
__________________________________________________________________________________________

#Trimming
__________________________________________________________________________________________

Trimmomatic

#!/bin/bash -e
#SBATCH -A massey03212
#SBATCH -J Trimmomatic
#SBATCH --time 250:00:00
#SBATCH --mem 32GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH -e Trimmomatic.err
#SBATCH -o Trimmomatic.out
#SBATCH --export NONE

export SLURM_EXPORT_ENV=ALL

# Set up working directories
cd /nesi/project/massey03212/AGRF/AGRF_CAGRF21098330_H3CGTDSX3/Trimmomatic

# Load Trimmomatic
module purge
module load Trimmomatic/0.39-Java-1.8.0_144

for file in *_R1.fastq.gz
  do
    base=$(basename ${file} _R1.fastq.gz)
    trimmomatic PE -threads 16 -phred33 ${file} ${base}_R2.fastq.gz \
                 ${base}_R1.trim.fastq.gz ${base}_R1un.trim.fastq.gz \
                 ${base}_R2.trim.fastq.gz ${base}_R2un.trim.fastq.gz \
                 ILLUMINACLIP:/scale_wlg_persistent/filesets/opt_nesi/CS400_centos7_bdw/nullarbor/2.0.20191013/conf/trimmomatic.fa:1:30:11 LEADING:10 TRAILING:10 MINLEN:80 HEADCROP:10 
done
__________________________________________________________________________________________

#Quality control
__________________________________________________________________________________________

FastQC

#!/bin/bash -e
#SBATCH -A massey03212
#SBATCH -J Fastqc
#SBATCH --time 72:00:00
#SBATCH --mem 64GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 2
#SBATCH -e fastqc.err
#SBATCH -o fastqc.out
#SBATCH --export NONE

export _JAVA_OPTIONS=-Djava.io.tmpdir=/nesi/nobackup/massey03212/javatmp
export _JAVA_OPTIONS='-Xmx12g'

# Set up working directories
cd /nesi/project/massey03212/AGRF/AGRF_CAGRF21098330_H3CGTDSX3/Trimmomatic

# Load Fastqc
module purge
module load FastQC/0.11.9

# Run FastQC

fastqc *.gz -o /nesi/project/massey03212/AGRF/AGRF_CAGRF21098330_H3CGTDSX3/Trimmomatic/fastqc_output

__________________________________________________________________________________________

#Removing the host contamination
__________________________________________________________________________________________

BBMap

Build index reference - Gorilla

#!/bin/bash -e
#SBATCH -A massey03212
#SBATCH -J host_filt_bbmap_index
#SBATCH -J Gorilla_ref
#SBATCH --time 00:40:00
#SBATCH --mem 23GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH -e host_filt_bbmap_index.err
#SBATCH -o host_filt_bbmap_index.out
#SBATCH --export NONE

export SLURM_EXPORT_ENV=ALL

cd /nesi/project/massey03212/AGRF/AGRF_CAGRF21098330_H3CGTDSX3/Trimmomatic/bbmap/Gorilla_wild

# Load BBMap module
module purge
module load BBMap/38.95-gimkl-2020a

# Build indexed reference file via BBMap
srun bbmap.sh ref=GCF_008122165.1_Kamilah_GGO_v0_genomic.fna -Xmx23g

__________________________________________________________________________________________

BBMap

Build index reference - Human

#!/bin/bash -e
#SBATCH -A massey03212
#SBATCH -J host_filt_bbmap_index
#SBATCH -J Human_ref
#SBATCH --time 00:40:00
#SBATCH --mem 23GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH -e host_filt_bbmap_index.err
#SBATCH -o host_filt_bbmap_index.out
#SBATCH --export NONE

export SLURM_EXPORT_ENV=ALL

cd /nesi/project/massey03212/AGRF/AGRF_CAGRF21098330_H3CGTDSX3/Trimmomatic/bbmap/Human_ref

# Load BBMap module
module purge
module load BBMap/38.95-gimkl-2020a

# Build indexed reference file via BBMap
srun bbmap.sh ref=GCF_000001405.39_GRCh38.p13_genomic.fna -Xmx23g

__________________________________________________________________________________________

BBMap

Build index reference - Cattle

#!/bin/bash -e
#SBATCH -A massey03212
#SBATCH -J host_filt_bbmap_index
#SBATCH -J Cattle_ref
#SBATCH --time 00:40:00
#SBATCH --mem 23GB
#SBATCH --ntasks 1
#SBATCH -e host_filt_bbmap_index.err
#SBATCH -o host_filt_bbmap_index.out
#SBATCH --export NONE

export SLURM_EXPORT_ENV=ALL

cd /nesi/project/massey03212/AGRF/AGRF_CAGRF21098330_H3CGTDSX3/Trimmomatic/bbmap/Cattle_ref

# Load BBMap module
module purge
module load BBMap/38.95-gimkl-2020a

# Build indexed reference file via BBMap
srun bbmap.sh ref=GCF_002263795.1_ARS-UCD1.2_genomic.fna -Xmx23g

__________________________________________________________________________________________

BBMap

Build index reference - Goat

#!/bin/bash -e
#SBATCH -A massey03212
#SBATCH -J host_filt_bbmap_index
#SBATCH -J Goat_ref
#SBATCH --time 00:40:00
#SBATCH --mem 23GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH -e host_filt_bbmap_index.err
#SBATCH -o host_filt_bbmap_index.out
#SBATCH --export NONE

export SLURM_EXPORT_ENV=ALL

cd /nesi/project/massey03212/AGRF/AGRF_CAGRF21098330_H3CGTDSX3/Trimmomatic/bbmap/Goat_ref

# Load BBMap module
module purge
module load BBMap/38.95-gimkl-2020a

# Build indexed reference file via BBMap
srun bbmap.sh ref=GCF_001704415.1_ARS1_genomic.fna -Xmx23g

__________________________________________________________________________________________

Removing host contamination - Gorilla:

#!/bin/bash -e
#SBATCH -A massey03212
#SBATCH -J bbmap_gor
#SBATCH --time 14:00:00
#SBATCH --mem 27GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20
#SBATCH -e bbmap_gorilla_wild.err
#SBATCH -o bbmap_gorilla_wild.out
#SBATCH --export NONE

export SLURM_EXPORT_ENV=ALL

# Set up working directories
cd /nesi/project/massey03212/AGRF/AGRF_CAGRF21098330_H3CGTDSX3/Trimmomatic/bbmap/Gorilla_habituated
mkdir -p host_filtered_reads/

# Load BBMap module
module purge
module load BBMap/38.95-gimkl-2020a

for file in *_R1.trim.fastq.gz
  do
    base=$(basename ${file} _R1.trim.fastq.gz)

# Run bbmap

srun bbmap.sh in1=${file} in2=${base}_R2.trim.fastq.gz -Xmx27g -t=20 \
minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 qtrim=rl trimq=10 untrim \
path=/nesi/project/massey03212/AGRF/AGRF_CAGRF21098330_H3CGTDSX3/Trimmomatic/bbmap/Gorilla_ref \
outu1=/nesi/project/massey03212/AGRF/AGRF_CAGRF21098330_H3CGTDSX3/Trimmomatic/bbmap/Gorilla_habituated/host_filtered_reads/${base}_R1_hostFilt.fastq.gz outu2=/nesi/project/massey03212/AGRF/AGRF_CAGRF21098330_H3CGTDSX3/Trimmomatic/bbmap/Gorilla_habituated/host_filtered_reads/${base}_R2_hostFilt.fastq.gz
done

__________________________________________________________________________________________

Removing host contamination - Human:

#!/bin/bash -e
#SBATCH -A massey03212
#SBATCH -J bbmap_Hum_Com
#SBATCH --time 14:00:00
#SBATCH --mem 27GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20
#SBATCH -e bbmap_human_Com.err
#SBATCH -o bbmap_human_Com.out
#SBATCH --export NONE

export SLURM_EXPORT_ENV=ALL

# Set up working directories
cd /nesi/project/massey03212/AGRF/AGRF_CAGRF21098330_H3CGTDSX3/Trimmomatic/bbmap/Human_community
mkdir -p host_filtered_reads/

# Load BBMap module
module purge
module load BBMap/38.95-gimkl-2020a

for file in *_R1.trim.fastq.gz
  do
    base=$(basename ${file} _R1.trim.fastq.gz)

# Run bbmap

srun bbmap.sh in1=${file} in2=${base}_R2.trim.fastq.gz -Xmx27g -t=20 \
minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 qtrim=rl trimq=10 untrim \
path=/nesi/project/massey03212/AGRF/AGRF_CAGRF21098330_H3CGTDSX3/Trimmomatic/bbmap/Human_ref \
out1=/nesi/project/massey03212/AGRF/AGRF_CAGRF21098330_H3CGTDSX3/Trimmomatic/bbmap/Human_community/host_filtered_reads/${base}_R1_hostFilt.fastq.gz out2=/nesi/project/massey03212/AGRF/AGRF_CAGRF21098330_H3CGTDSX3/Trimmomatic/bbmap/Human_community/host_filtered_reads/${base}_R2_hostFilt.fastq.gz
done

__________________________________________________________________________________________

Removing host contamination - Cattle:

#!/bin/bash -e
#SBATCH -A massey03212
#SBATCH -J bbmap_Cattle
#SBATCH --time 14:00:00
#SBATCH --mem 27GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20
#SBATCH -e bbmap_Cattle.err
#SBATCH -o bbmap_Cattle.out
#SBATCH --export NONE

export SLURM_EXPORT_ENV=ALL

# Set up working directories
cd /nesi/project/massey03212/AGRF/AGRF_CAGRF21098330_H3CGTDSX3/Trimmomatic/bbmap/Cattle
mkdir -p host_filtered_reads/

# Load BBMap module
module purge
module load BBMap/38.95-gimkl-2020a

for file in *_R1.trim.fastq.gz
  do
    base=$(basename ${file} _R1.trim.fastq.gz)

# Run bbmap

srun bbmap.sh in1=${file} in2=${base}_R2.trim.fastq.gz -Xmx27g -t=20 \
minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 qtrim=rl trimq=10 untrim \
path=/nesi/project/massey03212/AGRF/AGRF_CAGRF21098330_H3CGTDSX3/Trimmomatic/bbmap/Cattle \
out1=/nesi/project/massey03212/AGRF/AGRF_CAGRF21098330_H3CGTDSX3/Trimmomatic/bbmap/Cattle/host_filtered_reads/${base}_R1_hostFilt.fastq.gz out2=/nesi/project/massey03212/AGRF/AGRF_CAGRF21098330_H3CGTDSX3/Trimmomatic/bbmap/Cattle/host_filtered_reads/${base}_R2_hostFilt.fastq.gz
done

__________________________________________________________________________________________

Removing host contamination - Goat:

#!/bin/bash -e
#SBATCH -A massey03212
#SBATCH -J bbmap_Goat
#SBATCH --time 14:00:00
#SBATCH --mem 27GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20
#SBATCH -e bbmap_Goat.err
#SBATCH -o bbmap_Goat.out
#SBATCH --export NONE

export SLURM_EXPORT_ENV=ALL

# Set up working directories
cd /nesi/project/massey03212/AGRF/AGRF_CAGRF21098330_H3CGTDSX3/Trimmomatic/bbmap/Goat
mkdir -p host_filtered_reads/

# Load BBMap module
module purge
module load BBMap/38.95-gimkl-2020a

for file in *_R1.trim.fastq.gz
  do
    base=$(basename ${file} _R1.trim.fastq.gz)

# Run bbmap

srun bbmap.sh in1=${file} in2=${base}_R2.trim.fastq.gz -Xmx27g -t=20 \
minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 qtrim=rl trimq=10 untrim \
path=/nesi/project/massey03212/AGRF/AGRF_CAGRF21098330_H3CGTDSX3/Trimmomatic/bbmap/Goat \
out1=/nesi/project/massey03212/AGRF/AGRF_CAGRF21098330_H3CGTDSX3/Trimmomatic/bbmap/Goat/host_filtered_reads/${base}_R1_hostFilt.fastq.gz out2=/nesi/project/massey03212/AGRF/AGRF_CAGRF21098330_H3CGTDSX3/Trimmomatic/bbmap/Goat/host_filtered_reads/${base}_R2_hostFilt.fastq.gz
done

__________________________________________________________________________________________

#Assemblies
__________________________________________________________________________________________

Megahit

#!/bin/bash -e
#SBATCH -A massey03212
#SBATCH -J Megahit_HH
#SBATCH --time 144:00:00
#SBATCH --mem 128GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH -e megahit.err
#SBATCH -o megahit.out
#SBATCH --export NONE

export SLURM_EXPORT_ENV=ALL

# Set up working directories
cd /nesi/project/massey03212/AGRF/AGRF_CAGRF21098330_H3CGTDSX3/Trimmomatic/bbmap/Human_hospital/host_filtered_reads

# Load MEGAHIT module
module purge
module load MEGAHIT/1.1.4-gimkl-2018b-Python-2.7.16

# Run MEGAHIT

megahit -1 congo_gorilla_R1_hostFilt.fastq.gz -2 congo_gorilla_R2_hostFilt.fastq.gz -o megahit_assembly/
__________________________________________________________________________________________

Megahit - Multiple files

#!/bin/bash -e
#SBATCH -A massey03212
#SBATCH -J Megahit_GH
#SBATCH --time 01:00:00
#SBATCH --mem 128GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH -e megahit.err
#SBATCH -o megahit.out
#SBATCH --export NONE

export SLURM_EXPORT_ENV=ALL

# Load MEGAHIT module
module purge
module load MEGAHIT/1.1.4-gimkl-2018b-Python-2.7.16

for file in *_R1_hostFilt.fastq.gz
  do
    base=$(basename ${file} _R1_hostFilt.fastq.gz)

# Run MEGAHIT

megahit -1 ${file} -2 ${base}_R2_hostFilt.fastq.gz -o ${base}_megahit_out

done

__________________________________________________________________________________________

#!/bin/bash -e
#SBATCH --account       massey03212
#SBATCH --job-name      Megahit_CA
#SBATCH --time          36:00:00
#SBATCH --mem           128GB
#SBATCH --cpus-per-task 16
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out

# Load MEGAHIT module
module purge
module load MEGAHIT/1.1.4-gimkl-2018b-Python-2.7.16

for file in *_R1_hostFilt.fastq.gz
  do
    base=$(basename ${file} _R1_hostFilt.fastq.gz)

# Run MEGAHIT

megahit -1 ${file} -2 ${base}_R2_hostFilt.fastq.gz -o ${base}_megahit_out

done

__________________________________________________________________________________________

SPAdes

#!/bin/bash -e
#SBATCH -A massey03212
#SBATCH -J spades_GW
#SBATCH --time 168:00:00
#SBATCH --mem 1600GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 32
#SBATCH -e spades_gw.err
#SBATCH -o spades_gw.out
#SBATCH --export NONE

export SLURM_EXPORT_ENV=ALL

module load SPAdes/3.15.0-gimkl-2020a

spades.py --meta --memory 1500 -k 33,55,77,99,121 -t 32 -1 gor_wild_R1.fastq.gz -2 gor_wild_R2.fastq.gz -o spades_GW_out/
__________________________________________________________________________________________

#Binning
__________________________________________________________________________________________

Seqmagick

module load seqmagick/0.7.0-gimkl-2018b-Python-3.7.3

seqmagick convert --min-length 1000 final.contigs.fa cattle.m1000.fna

seqmagick convert --min-length 1000 contigs.fasta gw2.m1000.fna

seqmagick convert --min-length 1000 contig.fa ug405.m1000.fna
__________________________________________________________________________________________

#Mapping

Bowtie2 - index

#!/bin/bash -e
#SBATCH --account       massey03212
#SBATCH --job-name      BT2_index_CA
#SBATCH --time          00:50:00
#SBATCH --mem           16GB
#SBATCH --cpus-per-task 16
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out

module load Bowtie2/2.4.5-GCC-11.3.0

bowtie2-build --threads $SLURM_CPUS_PER_TASK cattle.m1000.fna cattle_bt2
__________________________________________________________________________________________

Bowtie2 - mapping

#!/bin/bash -e
#SBATCH --account       massey03212
#SBATCH --job-name      bt2_cattle
#SBATCH --time          08:00:00
#SBATCH --mem           08GB
#SBATCH --cpus-per-task 16
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out

module load Bowtie2/2.4.5-GCC-11.3.0

bowtie2 --minins 200 --maxins 800 --threads $SLURM_CPUS_PER_TASK --sensitive -x cattle_bt2 -1 cattle_R1_hostFilt.fastq.gz -2 cattle_R2_hostFilt.fastq.gz -S cattle.bam

__________________________________________________________________________________________

Samtools

#!/bin/bash -e
#SBATCH --account       massey03212
#SBATCH --job-name      Samtools_GW
#SBATCH --time          02:00:00
#SBATCH --mem           24GB
#SBATCH --cpus-per-task 16
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out

module purge
module load SAMtools/1.15.1-GCC-11.3.0

samtools sort -@ $SLURM_CPUS_PER_TASK gw.bam -o gw.sorted.bam

samtools index -@ $SLURM_CPUS_PER_TASK gw.sorted.bam

__________________________________________________________________________________________

Binning

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

New:

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
         
If it doesnt work: runMetaBat.sh -t $SLURM_CPUS_PER_TASK -m 1500 CA.m1000.fna cattle.sorted.bam

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

and Maxbin again:

for bin_path in maxbin/*.fasta;
do
    bin_name=$(basename ${bin_path} .fasta)
    grep ">" ${bin_path} | sed 's/>//g' | sed "s/$/\t${bin_name}/g" >> maxbin_associations.txt
done

Concoct

for bin_path in concoct_cattle_output/fasta_bins/*.fa;
do 
    bin_name=$(basename ${bin_path} .fa)
    grep ">" ${bin_path} | sed 's/>//g' | sed "s/$/\t${bin_name}/g" >> concoct_associations.txt 
done

__________________________________________________________________________________________

#Dereplication
__________________________________________________________________________________________

DasTool

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

__________________________________________________________________________________________

#Assessing the quality of genomes recovered
__________________________________________________________________________________________

CheckM

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
                  
__________________________________________________________________________________________

#!/bin/bash -e
#SBATCH -A massey03212
#SBATCH -J checkM_Campy
#SBATCH --time 03:00:00
#SBATCH --mem 128GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 32
#SBATCH -e checkM.err
#SBATCH -o checkM.out
#SBATCH --export NONE

module purge
module load CheckM/1.2.1-gimkl-2022a-Python-3.10.5
export PYTHONNOUSERSITE=1

checkm lineage_wf -t 10 --pplacer_threads 32 -x fa \
                  --tab_table -f checkm.txt \
                  dastool_out/_DASTool_bins/ checkm_out/

__________________________________________________________________________________________

#Taxonomic analysis
__________________________________________________________________________________________

GTDB-Tk

#!/bin/bash -e
#SBATCH --account       massey03212
#SBATCH --job-name      GTDBTk_GW
#SBATCH --time          03:00:00
#SBATCH --mem           128GB
#SBATCH --cpus-per-task 16
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out

module purge
module load GTDB-Tk/2.1.0-gimkl-2020a-Python-3.9.9

gtdbtk classify_wf -x fa --cpus 20 --genome_dir dastool_out/_DASTool_bins/ --out_dir gtdbtk_out/



New:

#!/bin/bash -e
#SBATCH --account       massey03212
#SBATCH --job-name      GTDB-Tk_4
#SBATCH --time          05:00:00
#SBATCH --mem           128GB
#SBATCH --cpus-per-task 16
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out

module purge
module load GTDB-Tk/2.1.0-gimkl-2020a-Python-3.9.9

gtdbtk classify_wf -x fa --cpus 16 --genome_dir _DASTool_bins/ --out_dir gtdbtk_out/
gtdbtk classify_wf -x fa --cpus 20 --genome_dir metabat/ --out_dir gtdbtk_metabat_out/
gtdbtk classify_wf -x fa --cpus 20 --genome_dir dereplicated_genomes/ --out_dir gtdbtk_out/

__________________________________________________________________________________________

#Other analysis strategies
__________________________________________________________________________________________
Phylophlan 3

#!/bin/bash -e
#SBATCH -A massey03212
#SBATCH -J Phylophlan3
#SBATCH --time 00:40:00
#SBATCH --mem 2GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 08
#SBATCH -e phylophlan.err
#SBATCH -o phylophlan.out
#SBATCH --export NONE

module load PhyloPhlAn/3.0.1-gimkl-2020a-Python-3.8.2

phylophlan -i brachyspira -d /opt/nesi/db/PhyloPhlAn --diversity medium -f configuration_file

__________________________________________________________________________________________

Kraken 2

#!/bin/bash -e
#SBATCH -A massey03212
#SBATCH -J Kra_Brac_GW
#SBATCH --time 07:00:00
#SBATCH --mem 64GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 02
#SBATCH -e kraken.err
#SBATCH -o kraken.out
#SBATCH --export NONE

module purge
module load Kraken2/2.1.1-GCC-9.2.0
module load Bracken/2.6.0-GCCcore-9.2.0

for file in *_R1_hostFilt.fastq.gz
  do
    base=$(basename ${file} _R1_hostFilt.fastq.gz)
    
# Run Kraken2

kraken2 --db /opt/nesi/db/Kraken2/standard/ --use-names --paired ${file} ${base}_R2_hostFilt.fastq.gz --output ${base}.txt --report ${base}.kreport

# Run Bracken2

bracken -d /opt/nesi/db/Kraken2/standard/ -i ${base}.kreport -o ${base}_bracken.report -t 16 -r 100

done

__________________________________________________________________________________________

Manual curation
__________________________________________________________________________________________

VizBin:
__________________________________________________________________________________________

#First step for input files:

mkdir original_bins

for file in dastool_out/_DASTool_bins/*;
do
    # Copy and rename bin file
    cp ${file} original_bins/bin_${i}.fna
    # extract bin ID
    binID=$(basename original_bins/bin_${i}.fna .fna)
    # Add bin ID to sequence headers
    sed -i -e "s/>/>${binID}_/g" original_bins/bin_${i}.fna
    # Increment i
    ((i+=1))
done

#Secong step for input files:

module purge
module load CONCOCT/1.1.0-gimkl-2020a-Python-3.8.2

mkdir bins_20k

# loop over .fna files to generate chopped (fragmented) files using CONCONT's cut_up_fasta.py
for bin_file in original_bins/*;
do
    bin_name=$(basename ${bin_file} .fna)
    cut_up_fasta.py -c 20000 -o 0 --merge_last ${bin_file} > bins_20k/${bin_name}.20k.fna
done

#concatenate

cat bins_20k/*.fna > all_bins.fna

#Mapping

module load Bowtie2/2.4.5-GCC-11.3.0

bowtie2-build all_bins.fna gw_brachy_bt2

Bowtie2:

#!/bin/bash -e
#SBATCH --account       massey03212
#SBATCH --job-name      bt2_cattle
#SBATCH --time          04:00:00
#SBATCH --mem           08GB
#SBATCH --cpus-per-task 12
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out

module purge
module load Bowtie2/2.4.5-GCC-11.3.0
module load SAMtools/1.15.1-GCC-11.3.0

bowtie2 --minins 200 --maxins 800 --threads $SLURM_CPUS_PER_TASK --sensitive -x gw_brachy_bt2 -1 gw_R1_hostFilt.fastq.gz -2 gw_R2_hostFilt.fastq.gz -S gw_brachy.sam

samtools sort -@ $SLURM_CPUS_PER_TASK -o gw_brachy.bam gw_brachy.sam

_______________________________________________________

Samtools:

samtools sort -@ 10 -o gw_brachy.bam gw_brachy.sam

#Generate coverage table of subcontigs:

module purge
module load MetaBAT/2.15-GCC-11.3.0

# calculate coverage table
jgi_summarize_bam_contig_depths --outputDepth gw_brachy_cov.txt gw_brachy.bam

#Generate annotation table for VizBin:

# Set up annotation file headers
echo "coverage,label,length" > all_bins.sample1.vizbin.ann

# loop through bin .fna files
for bin_file in *.fna; do
    # extract bin ID
    binID=$(basename ${bin_file} .fna)
    # loop through each sequence header in bin_file
    for header in `grep ">" ${bin_file}`; do
        contigID=$(echo ${header} | sed 's/>//g')
        # identify this line from the coverage table (gw_brachy_cov.txt), and extract contigLen (column 2) and coverage for sample1.bam (column 4)
        contigLen=$(grep -P "${contigID}\t" gw_brachy_cov.txt | cut -f2)
        sample1_cov=$(grep -P "${contigID}\t" gw_brachy_cov.txt | cut -f4)
        # Add to vizbin.ann file
        echo "${sample1_cov},${binID},${contigLen}" >> all_bins.sample1.vizbin.ann
    done
done



# Other scripts
__________________________________________________________________________________________

#!/bin/bash -e
#SBATCH -A massey03212
#SBATCH -J MAFFT
#SBATCH --time 07:00:00
#SBATCH --mem 64GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH -e mafft.err
#SBATCH -o mafft.out
#SBATCH --export NONE

module purge
module load MAFFT/7.487-gimkl-2020a-with-extensions

mafft --auto --anysymbol brachy_hosts_ref.fasta > brachy_hosts_ref.al

module load prodigal/2.6.3-GCC-9.2.0
prodigal -i metabat.294.contigs.fa -o gene.coords.gbk -a protein.translations.faa

module load MAFFT/7.487-gimkl-2020a-with-extensions
mafft protein.translations.faa > brachyspira_proteins.al

module load FastTree/2.1.11-GCCcore-9.2.0
fasttree brachyspira_proteins.al > tree

module load prokka/1.14.5-GCC-9.2.0
prokka C_curvus.fasta --outdir C_curvus

module load MAFFT/7.487-gimkl-2020a-with-extensions
mafft --auto --anysymbol gltA.fasta > gltA.al

module load FastTree/2.1.11-GCCcore-9.2.0
fasttree gtdbtk.bac120.msa.fasta > tree

module load Roary/3.13.0-gimkl-2020a
roary *.gff

__________________________________________________________________________________________

#!/bin/bash -e
#SBATCH -A massey03212
#SBATCH -J Roary_Brach
#SBATCH --time 01:40:00
#SBATCH --mem 64GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH -e roary.err
#SBATCH -o roary.out
#SBATCH --export NONE

module purge

module load Roary/3.13.0-gimkl-2020a

roary *.gff

or

module purge

module load Roary/3.13.0-gimkl-2020a

roary -f roary -e -n -v *.gff

_______________________________________________________

#!/bin/bash -e
#SBATCH --account       massey03212
#SBATCH --job-name      Roary_Brachy
#SBATCH --time          02:20:00
#SBATCH --mem           1GB
#SBATCH --cpus-per-task 08
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out

module purge

module load Roary/3.13.0-gimkl-2020a

roary -p $SLURM_CPUS_PER_TASK -f roary -e -n -v *.gff
_______________________________________________________

module purge

module load FastTree/2.1.11-GCCcore-9.2.0

FastTree -nt -gtr core_gene_alignment.aln > my_tree.newick

module purge

module load Python/3.10.5-gimkl-2022a

python roary_plots.py my_tree.newick gene_presence_absence.csv

__________________________________________________________________________________________

#!/bin/bash -e
#SBATCH -A massey03212
#SBATCH -J IQ-Tree_E.col
#SBATCH --time 01:50:00
#SBATCH --mem 64GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH -e iqtree.err
#SBATCH -o iqtree.out
#SBATCH --export NONE

module purge

module load IQ-TREE/2.2.0.5-gimpi-2022a
iqtree2 -s concatenated_gene_aa_alignment.fasta -T 16 -mem 64G

# With bootstrap:

#!/bin/bash -e
#SBATCH -A massey03212
#SBATCH -J IQTree_dRep_campy
#SBATCH --time 168:00:00
#SBATCH --mem 64GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH -e iqtree.err
#SBATCH -o iqtree.out
#SBATCH --export NONE

module purge
module load IQ-TREE/2.2.0.5-gimpi-2022a

iqtree2 -s gtdbtk.bac120.user_msa.fasta.gz -m TEST -b 1000 -T 16 -mem 64G

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

iqtree2 -s gtdbtk.bac120.user_msa.fasta.gz -T $SLURM_CPUS_PER_TASK -mem 32G

__________________________________________________________________________________________

#!/bin/bash -e
#SBATCH -A massey03212
#SBATCH -J Busco_HH_Campy
#SBATCH --time 07:00:00
#SBATCH --mem 64GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH -e busco.err
#SBATCH -o busco.out
#SBATCH --export NONE

module purge
module load BUSCO/5.2.2-gimkl-2020a

busco -i HH.faa --auto-lineage-prok -o HH_campy -m proteins

__________________________________________________________________________________________

Cognac

#!/bin/bash -e
#SBATCH -A massey03212
#SBATCH -J Cognac
#SBATCH --time 00:40:00
#SBATCH --mem 8GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH -e cognac.err
#SBATCH -o cognac.out
#SBATCH --export NONE

module purge
module load R/4.1.0-gimkl-2020a
module load MAFFT/7.487-gimkl-2020a-with-extensions
module load CD-HIT/4.8.1-gimkl-2017a

R

library(cognac)

algnEnv = cognac(
  fastaDir   = "/nesi/nobackup/massey03212/cognac/E_coli/fasta",
  featureDir = "/nesi/nobackup/massey03212/cognac/E_coli/gff",
  threadVal  = 4
  )

__________________________________________________________________________________________

PROKKA

#!/bin/bash -e
#SBATCH --account	massey03212
#SBATCH --job-name	Prokka_Campy_Heli
#SBATCH --time          01:00:00
#SBATCH --mem           1GB
#SBATCH --cpus-per-task 04
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out

module purge
module load prokka/1.14.5-GCC-9.2.0

for file in *.fa
  do
    base=$(basename ${file} .fa)

# Run PROKKA

prokka ${file} -o ${base} --cpus $SLURM_CPUS_PER_TASK

done


module load prokka/1.14.5-GCC-9.2.0

for file in *.fa
  do
    base=$(basename ${file} .fa)

# Run PROKKA

prokka ${file} -o ${base} --cpus $SLURM_CPUS_PER_TASK

done

___________________________________________________________________________________________

#!/bin/bash
#SBATCH -A uoa00446
#SBATCH -J 2-drep-refined-after-removing-duplicated-contigs
#SBATCH --time 24:00:00
#SBATCH --mem=50GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --requeue
#SBATCH -e slurm.2-drep-refined-after-removing-duplicated-contigs.%A.err
#SBATCH -o slurm.2-drep-refined-after-removing-duplicated-contigs%A.out

source /nesi/project/uoa02469/Software/PipelineSetup/binning.sh

cd /nesi/nobackup/uoa00446/wgscesar/8i.refined-bins-analysis/1b-quality-check-binswithduplicatedcontigs-drep-checkm-gtdbtk/1-checkm-refined-after-duplicates-removal
echo "genome,completeness,contamination" > refined-after-duplicate-removal.dRep.genomeInfo
cut -f1,12 checkm.txt | sed 's/\t/.fa,/g' > p1.txt
cut -f13 checkm.txt > p2.txt
paste -d "," p1.txt p2.txt | tail -n+2 >> refined-after-duplicate-removal.dRep.genomeInfo
rm p*.txt

##Run dRep
module load Python/3.7.3-gimkl-2018b
module load Mash/2.1-gimkl-2018b
dRep dereplicate /nesi/nobackup/uoa00446/wgscesar/8i.refined-bins-analysis/1b-quality-check-binswithduplicatedcontigs-drep-checkm-gtdbtk/2-drep -g /nesi/nobackup/uoa00446/wgscesar/8i.refined-bins-analysis/1a-quality_bins_after_duplicate_check/*.fa -p 6 --genomeInfo /nesi/nobackup/uoa00446/wgscesar/8i.refined-bins-analysis/1b-quality-check-binswithduplicatedcontigs-drep-checkm-gtdbtk/1-checkm-refined-after-duplicates-removal/refined-after-duplicate-removal.dRep.genomeInfo


source:

# Dependency module loads
module load BBMap/37.93-gimkl-2017a
module load BEDTools/2.26.0-gimkl-2017a
module load Bowtie/1.2.0-gimkl-2017a
module load Bowtie2/2.3.2-gimkl-2017a
module load SAMtools/1.8-gimkl-2017a
module load MetaBAT/2.13-GCC-7.4.0
module load MUMmer/3.23-gimkl-2017a

module load DAS_Tool/1.1.1-gimkl-2018b-R-3.6.1
module load MaxBin/2.2.6-gimkl-2018b-Perl-5.28.1

module load CONCOCT/1.0.0-gimkl-2018b-Python-2.7.16
module load Python/3.7.3-gimkl-2018b # Main python session for dRep

module load CheckM/1.0.13-gimkl-2018b-Python-2.7.16