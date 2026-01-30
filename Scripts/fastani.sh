#!/bin/bash -e
#SBATCH --account       massey03212
#SBATCH --job-name      FastANI_Campy
#SBATCH --time          24:00:00
#SBATCH --mem           1GB
#SBATCH --cpus-per-task 5
#SBATCH --error         %x_%j.err
#SBATCH --output        %x_%j.out

module purge
module load FastANI/1.33-intel-2022a

# Paths
GENOME_DIR="/nesi/nobackup/massey03212/Nat_comms/Genomic_analysis/02_gtdbtk/new_dataset_corrected"
OUTPUT="fastani_results.tsv"
GENOME_LIST="genome_list.txt"

# Generate genome list with absolute paths
find "$GENOME_DIR" -type f \( -name "*.fa" -o -name "*.fna" -o -name "*.fasta" \) -exec realpath {} \; > "$GENOME_LIST"

# Read genomes into array
mapfile -t GENOMES < "$GENOME_LIST"
NUM_GENOMES=${#GENOMES[@]}
echo "Processing $NUM_GENOMES genomes"

# Create output directory for temp files
mkdir -p tmp_results

# Main ANI calculation loop
for (( i=0; i<NUM_GENOMES; i++ )); do
  for (( j=i; j<NUM_GENOMES; j++ )); do
    # Create unique output filename
    OUT_FILE="tmp_results/ani_${i}_${j}.txt"
    
    # Run FastANI with error handling
    if ! fastANI -q "${GENOMES[$i]}" \
                 -r "${GENOMES[$j]}" \
                 --minFraction 0.1 \
                 --fragLen 1000 \
                 -o "$OUT_FILE" 2>> fastani_errors.log
    then
      echo "ERROR comparing ${GENOMES[$i]} vs ${GENOMES[$j]}" >> fastani_errors.log
      # Create empty file to prevent missing comparisons
      touch "$OUT_FILE"
    fi
  done
done

# Combine results
cat tmp_results/*.txt > "$OUTPUT"

echo "Completed ANI calculations. Output: $OUTPUT"
echo "Error log: fastani_errors.log"
