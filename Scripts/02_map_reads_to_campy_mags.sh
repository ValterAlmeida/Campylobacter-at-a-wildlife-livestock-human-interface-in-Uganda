#!/bin/bash -e
#SBATCH --account=massey03212
#SBATCH --job-name=campy_mag_relabund
#SBATCH --output=logs/campy_mag_relabund_%j.out
#SBATCH --error=logs/campy_mag_relabund_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=14

set -o pipefail

# ==============================================================================
# Competitive mapping of host-filtered reads to the 23 dereplicated Campylobacter
# MAGs, one library at a time (for-loop, no array).
#
# Per library, the script writes the number of mapped reads for every MAG
# (all alignments, and alignments with >= 95% identity) plus the library size
# in read pairs, taken from the Bowtie2 summary. No BAM files are kept, to spare
# nobackup space.
#
# Counting unit: individual reads (each mate counted separately). Library size
# in reads = 2 x total_pairs.
#
# Re-running is safe: libraries with an existing counts file are skipped.
#
# Defaults are the 553 individual libraries. For other sets, override at submit:
#   READS_DIR=/path/to/pooled R1_SUFFIX=_R1.fastq.gz R2_SUFFIX=_R2.fastq.gz \
#   TAG=pooled sbatch 02_map_reads_to_campy_mags.sh
# ==============================================================================

module purge
module load Bowtie2/2.5.4-GCC-12.3.0

BASE_DIR="/nesi/nobackup/massey03212/Nat_comms_05_26/Genomic_analysis/18_MAG_relative_abundance"
INDEX_PREFIX="${BASE_DIR}/index/campy_mags_23_index"
MAG_LIST="${BASE_DIR}/index/campy_mags_23.txt"

READS_DIR="${READS_DIR:-/nesi/nobackup/massey03255/Valter_06_26/Orion_fastq}"
R1_SUFFIX="${R1_SUFFIX:-_hostFilt_R1.fastq}"
R2_SUFFIX="${R2_SUFFIX:-_hostFilt_R2.fastq}"
TAG="${TAG:-individual}"

RES_DIR="${BASE_DIR}/results_${TAG}"
mkdir -p "${RES_DIR}/counts" "${RES_DIR}/bt2_logs"

if [ ! -s "${INDEX_PREFIX}.1.bt2" ]; then
    echo "ERROR: index not found at ${INDEX_PREFIX}. Run 01_build_campy_mag_index.sh first."
    exit 1
fi

echo "Reads dir:  ${READS_DIR}"
echo "Suffixes:   ${R1_SUFFIX} / ${R2_SUFFIX}"
echo "Tag:        ${TAG}"
echo "Output dir: ${RES_DIR}"

n_found=0
n_done=0
n_skipped=0

for R1 in "${READS_DIR}"/*"${R1_SUFFIX}"; do
    [ -e "$R1" ] || { echo "ERROR: no files matching *${R1_SUFFIX} in ${READS_DIR}"; exit 1; }

    sample=$(basename "$R1" "$R1_SUFFIX")
    R2="${READS_DIR}/${sample}${R2_SUFFIX}"
    OUT="${RES_DIR}/counts/${sample}.tsv"
    LOG="${RES_DIR}/bt2_logs/${sample}.log"
    n_found=$((n_found + 1))

    if [ -s "$OUT" ]; then
        n_skipped=$((n_skipped + 1))
        continue
    fi
    if [ ! -f "$R2" ]; then
        echo "WARNING: no R2 for ${sample}. Skipping."
        continue
    fi

    echo "[$(date +%H:%M:%S)] ${sample}"

    bowtie2 \
        --minins 200 \
        --maxins 800 \
        --sensitive \
        --no-unal \
        --threads "$SLURM_CPUS_PER_TASK" \
        -x "$INDEX_PREFIX" -1 "$R1" -2 "$R2" 2> "$LOG" | \
    awk -F'\t' -v OFS='\t' -v s="$sample" -v magfile="$MAG_LIST" '
        BEGIN { while ((getline line < magfile) > 0) mags[++n] = line }
        /^@/ { next }
        {
            if (int($2 / 4) % 2 == 1) next          # unmapped mate
            split($3, a, "|"); m = a[1]
            all[m]++
            nm = -1
            for (i = 12; i <= NF; i++) if ($i ~ /^NM:i:/) { nm = substr($i, 6) + 0; break }
            if (nm >= 0 && nm / length($10) <= 0.05) hi[m]++
        }
        END { for (i = 1; i <= n; i++) { m = mags[i]; print s, m, all[m] + 0, hi[m] + 0 } }
    ' > "${OUT}.tmp"

    pairs=$(awk 'NR == 1 {print $1; exit}' "$LOG")
    if [ -z "$pairs" ]; then
        echo "ERROR: could not read the library size from ${LOG}"
        exit 1
    fi
    awk -F'\t' -v OFS='\t' -v p="$pairs" '{print $0, p}' "${OUT}.tmp" > "$OUT"
    rm -f "${OUT}.tmp"
    n_done=$((n_done + 1))
done

echo "Libraries found: ${n_found}; mapped now: ${n_done}; already done: ${n_skipped}"

LONG="${RES_DIR}/mag_counts_long_${TAG}.tsv"
{
    printf 'sample\tmag\tmapped_reads\tmapped_reads_id95\ttotal_pairs\n'
    cat "${RES_DIR}"/counts/*.tsv
} > "$LONG"

echo "Libraries in the final table: $(cut -f1 "$LONG" | tail -n +2 | sort -u | wc -l)"
echo "Rows (libraries x 23 MAGs): $(tail -n +2 "$LONG" | wc -l)"
echo "Final table: ${LONG}"
