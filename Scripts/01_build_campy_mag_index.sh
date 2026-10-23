#!/bin/bash -e
#SBATCH --account=massey03212
#SBATCH --job-name=campy_mag_index
#SBATCH --output=logs/campy_mag_index_%j.out
#SBATCH --error=logs/campy_mag_index_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=4

set -o pipefail

# ==============================================================================
# Build one Bowtie2 index from the 23 dereplicated Campylobacter MAGs used for
# Supplementary Fig. 1 (the bin_id column of bin_coverage_taxonomy.tsv).
#
# Every contig header is rewritten as  MAGID|contig  so that reads can be
# assigned to a MAG by splitting RNAME on "|".
#
# The script stops if any of the 23 MAGs is missing from DEREP_DIR. That is the
# check that DEREP_DIR is the dereplicated set used in the manuscript.
# ==============================================================================

module purge
module load Bowtie2/2.5.4-GCC-12.3.0

BASE_DIR="/nesi/nobackup/massey03212/Nat_comms_05_26/Genomic_analysis/18_MAG_relative_abundance"
DEREP_DIR="${DEREP_DIR:-/nesi/project/massey03212/Nat_comms/Genomic_analysis/10_Metabolic_functions/04_gene_annotation_and_coverage/mags_coverage/mags}"
INDEX_DIR="${BASE_DIR}/index"
MAG_LIST="${INDEX_DIR}/campy_mags_23.txt"
COMBINED_FASTA="${INDEX_DIR}/campy_mags_23_combined.fasta"
INDEX_PREFIX="${INDEX_DIR}/campy_mags_23_index"

mkdir -p "$INDEX_DIR"

echo "=== Step 1: write the list of 23 MAGs (bin_id from bin_coverage_taxonomy.tsv) ==="
cat > "$MAG_LIST" << 'EOF'
11.gh06
160_sub.gw
178_sub.hc18
18.gh
211.go29
228_sub.gh
56.hh12
94_sub.hh13
bin.11.go26
bin.111.go27
bin.156.hh15
bin.22.hc16
bin.304.hh
bin.512.hc
bin.63.hh11
bin.67.hc17
bin.98.hc17
maxbin.010.go28
maxbin.030.go30
maxbin.082_sub.hc18
maxbin.142.hh11
maxbin.144.hc16
maxbin.327.ca
EOF
echo "MAGs in list: $(wc -l < "$MAG_LIST")"

echo "=== Step 2: check that all 23 MAGs exist in the dereplicated set ==="
echo "DEREP_DIR: ${DEREP_DIR}"
if [ ! -d "$DEREP_DIR" ]; then
    echo "ERROR: DEREP_DIR does not exist. Aborting."
    exit 1
fi
echo "Genomes in DEREP_DIR: $(ls "$DEREP_DIR" | wc -l)"

missing=0
while read -r id; do
    if [ ! -f "${DEREP_DIR}/${id}.fa" ]; then
        echo "MISSING: ${id}.fa"
        missing=$((missing + 1))
    fi
done < "$MAG_LIST"

if [ "$missing" -gt 0 ]; then
    echo "ERROR: ${missing} of 23 MAGs not found in this dRep output."
    echo "This is probably not the dereplicated set used for Supplementary Fig. 1. Aborting."
    exit 1
fi
echo "OK: all 23 MAGs found."

ls "$DEREP_DIR" | sed 's/\.fa$//' | sort > "${INDEX_DIR}/derep_dir_genomes.txt"
sort "$MAG_LIST" > "${INDEX_DIR}/campy_mags_23.sorted.txt"
n_extra=$(comm -13 "${INDEX_DIR}/campy_mags_23.sorted.txt" "${INDEX_DIR}/derep_dir_genomes.txt" | wc -l)
echo "Genomes in DEREP_DIR that are NOT in the 23-MAG list: ${n_extra} (they are ignored)"

echo "=== Step 3: build the combined FASTA with MAGID|contig headers ==="
: > "$COMBINED_FASTA"
while read -r id; do
    awk -v id="$id" '/^>/ {split($0, h, " "); print ">" id "|" substr(h[1], 2); next} {print}' \
        "${DEREP_DIR}/${id}.fa" >> "$COMBINED_FASTA"
done < "$MAG_LIST"

n_contigs=$(grep -c "^>" "$COMBINED_FASTA")
n_dup=$(grep "^>" "$COMBINED_FASTA" | sort | uniq -d | wc -l)
echo "Contigs in combined FASTA: ${n_contigs}"
echo "Duplicated contig headers: ${n_dup}"
if [ "$n_dup" -gt 0 ]; then
    echo "ERROR: duplicated headers in combined FASTA. Aborting."
    exit 1
fi

echo "=== Step 4: MAG lengths (bp) for optional length normalisation ==="
awk '/^>/ {split(substr($0, 2), a, "|"); m = a[1]; next} {len[m] += length($0)} END {for (m in len) print m "\t" len[m]}' \
    "$COMBINED_FASTA" | sort > "${INDEX_DIR}/campy_mags_23_lengths.tsv"
cat "${INDEX_DIR}/campy_mags_23_lengths.tsv"

echo "=== Step 5: build the Bowtie2 index ==="
bowtie2-build --threads "$SLURM_CPUS_PER_TASK" "$COMBINED_FASTA" "$INDEX_PREFIX" \
    > "${INDEX_DIR}/bowtie2_build.log"

echo "=== Done. Index: ${INDEX_PREFIX} ==="
