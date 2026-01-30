#!/bin/bash
#SBATCH --account=massey03212
#SBATCH --job-name=kofam_scan
#SBATCH --time=6:00:00
#SBATCH --mem=5GB
#SBATCH --cpus-per-task=12
#SBATCH --output=00.slurm.out/slurm-kofam-pipeline-%j.out
#SBATCH --error=00.slurm.out/slurm-kofam-pipeline-%j.err

# --- Define Directories and Variables ---

export KOFAM_PREP_DIR="./kofamscan_ready"
export OUT="./02.genome.kofamscan"
export INPUT_FILE="all_genomes_combined.faa"
export MAPPING_FILE="genome_name_map.txt"

# Create necessary directories
mkdir -p 00.slurm.out
mkdir -p ${KOFAM_PREP_DIR}
mkdir -p ${OUT}

# Define software and database paths (Updated to /nesi/project)
export KOFAM_BIN="/nesi/project/uoa02469/Software/kofam_scan_v1.3.0/bin/exec_annotation" 
export KO_LIST_DB="/nesi/project/uoa02469/Software/kofam_scan_v1.3.0/db/ko_list"
export PROFILES_DB="/nesi/project/uoa02469/Software/kofam_scan_v1.3.0/db/profiles"

# ==============================================================================
### 0. Create Temporary Renaming Map
# The map links the messy filename prefix (key) to the clean KEGGDecoder label (value).
# ==============================================================================

cat <<EOF > ${MAPPING_FILE}
178_sub.hc18          Ug.Hnc.C.sp9
94_sub.hh13           Ug.Hc13.C.sp9
GCA_022723245.1_ASM2272324v1_genomic  GCA022723245.1.C.sp9
GCA_022734655.1_ASM2273465v1_genomic  GCA022734655.1.C.sp9
GCA_037334345.1_ASM3733434v1_genomic  GCA037334345.1.C.sp9
GCA_900539255.1_UMGS113_genomic       GCA900539255.1.C.sp9
GCA_900772045.1_SRS820613_28_genomic  GCA900772045.1.C.sp9
GCA_902479625.1_UHGG-TPA_MGYG-HGUT-02813_genomic  GCA902479625.1.C.sp9
GCA_934667445.1_ERR7738533_bin.68_genomic GCA934667445.1.C.sp9
GCA_934692235.1_ERR7745335_bin.520_genomic GCA934692235.1.C.sp9
GCA_934726075.1_ERR7738211_bin.102_genomic GCA934726075.1.C.sp9
GCA_934726175.1_ERR7745719_bin.258_genomic GCA934726175.1.C.sp9
GCA_938045135.1_ERR2619723_bin.54_CONCOCT_v1.1_MAG_genomic GCA938045135.1.C.sp9
KY206_7               KY206.7.C.sp9
bin.156.hh15          Ug.Hc15.C.sp9
bin.98.hc17           Ug.Hnc17.C.sp9
maxbin.142.hh11       Ug.Hc11.C.sp9
maxbin.144.hc16       Ug.Hnc16.C.sp9
56.hh12               Ug.Hc12.C.inf
BF181_maxbin.011      BF181.C.inf
GCA_009827235.1_ASM982723v1_genomic GCA009827235.1.C.inf
GCA_022739255.1_ASM2273925v1_genomic GCA022739255.1.C.inf
GCA_037302995.1_ASM3730299v1_genomic GCA037302995.1.C.inf
GCA_900766765.1_SRS475627_10_genomic GCA900766765.1.C.inf
GCA_902465685.1_UHGG-TPA_MGYG-HGUT-00755_genomic GCA902465685.1.C.inf
GCA_934717465.1_ERR7738353_bin.56_genomic GCA934717465.1.C.inf
GCA_934728225.1_ERR7745561_bin.49_genomic GCA934728225.1.C.inf
GCA_934728275.1_ERR7738232_bin.112_genomic GCA934728275.1.C.inf
GCA_937975125.1_ERR2619726_bin.91_CONCOCT_v1.1_MAG_genomic GCA937975125.1.C.inf
GCA_958349105.1_ERR10960912_bin.9_MetaWRAP_v1.3_MAG_genomic GCA958349105.1.C.inf
GCA_958370705.1_SRR17382065_bin.37_MetaWRAP_v1.3_MAG_genomic GCA958370705.1.C.inf
GCA_963591955.1_H07622-L1_cleanbin_000012_GreatApes_genomic GCA963591955.1.C.inf
GCF_013416015.1_ASM1341601v1_genomic GCF013416015.1.C.inf
bin.22.hc16             Ug.Hnc16.C.inf
bin.304.hh              Ug.Hc.C.inf
bin.512.hc              Ug.Hnc.C.inf
bin.63.hh11             Ug.Hc11.C.inf
bin.67.hc17             Ug.Hnc17.C.inf
maxbin.082_sub.hc18     Ug.Hnc18.C.inf
EOF

# ==============================================================================
### 1. File Management: Header Renaming and Concatenation
# ==============================================================================

echo "--- 1. Starting File Preparation (Renaming Headers with Clean Prefix) ---"

# 1. Read the mapping file into an associative array for easy lookup
declare -A NAME_MAP
while read -r OLD_NAME NEW_NAME; do
    # Remove the .genes.no_metadata.faa suffix from the OLD_NAME for lookup
    # Only map if both names are non-empty
    if [[ -n "$OLD_NAME" && -n "$NEW_NAME" ]]; then
        BASE_NAME=$(echo "$OLD_NAME" | sed 's/\.genes\.no_metadata\.faa//g')
        NAME_MAP["$BASE_NAME"]="$NEW_NAME"
    fi
done < "${MAPPING_FILE}"

# 2. Renaming Headers
# List of all input files
INPUT_FILES=(./*.genes.no_metadata.faa)

for IN_FILE in "${INPUT_FILES[@]}"; do
    # Get the original base name from the filename
    BASE_NAME=$(basename "${IN_FILE}" | sed 's/\.genes\.no_metadata\.faa//g')
    
    # Look up the new, clean name
    CLEAN_NAME="${NAME_MAP["$BASE_NAME"]}"
    
    # Ensure a mapping exists before proceeding
    if [[ -z "$CLEAN_NAME" ]]; then
        echo "Error: No clean name found for base name: ${BASE_NAME}. Skipping."
        continue
    fi

    OUT_FILE="${KOFAM_PREP_DIR}/${CLEAN_NAME}.faa"

    # AWK: Prepend the CLEAN_NAME to every FASTA header
    awk -v name="${CLEAN_NAME}" '
        /^>/ {
            # New Header: >CLEAN_NAME_ORIGINAL_HEADER
            print ">" name "_" substr($0, 2);
            next;
        }
        {
            print;
        }
    ' "${IN_FILE}" > "${OUT_FILE}"
done

# 3. Concatenation
rm -f ${INPUT_FILE} # Clean up any previous combined file
cat ${KOFAM_PREP_DIR}/*.faa > ${INPUT_FILE}

echo "Total sequences in combined file:" $(grep -c ">" ${INPUT_FILE})


# ==============================================================================
### 2. KofamScan Execution
# ==============================================================================

echo -e "\n--- 2. Starting KofamScan Execution ---"

## Purge all loaded modules, and load new software:
module purge
module load Ruby/2.6.1-gimkl-2018b
module load HMMER/3.2.1-gimkl-2018b
module load Parallel/20190622

# First run: -f mapper (concise output for KEGGDecoder)
${KOFAM_BIN} -o ${OUT}/kofam.combined.ortholog.txt --cpu 12 -E 0.0001 -k ${KO_LIST_DB} -p ${PROFILES_DB} -f mapper ${INPUT_FILE}

# Second run: -f detail (detailed output)
${KOFAM_BIN} -o ${OUT}/kofam.combined.ortholog.detail.txt --cpu 12 -E 0.0001 -k ${KO_LIST_DB} -p ${PROFILES_DB} -f detail ${INPUT_FILE}

echo "KofamScan finished."


# ==============================================================================
### 3. KEGGDecoder Execution
# ==============================================================================

echo -e "\n--- 3. Starting KEGGDecoder Analysis ---"

## Load modules for Python/KEGGDecoder:
module purge
module load Python/3.7.3-gimkl-2018b

# Install KEGGDecoder (You can comment this out after the first successful run)
pip install KEGGDecoder --user

# Define Input File (Output from the previous KofamScan step)
export KOFAM_OUTPUT_FILE=${OUT}/kofam.combined.ortholog.txt

## Run KEGGDecoder on kofamscan output:
# Static output
~/.local/bin/KEGG-decoder -i ${KOFAM_OUTPUT_FILE} -o ${OUT}/kegg.kofam.combined.ortholog.static.txt -v static

# Interactive output
~/.local/bin/KEGG-decoder -i ${KOFAM_OUTPUT_FILE} -o ${OUT}/kegg.kofam.combined.ortholog.interactive.txt -v interactive

echo "KEGGDecoder finished."


# ==============================================================================
### 4. Cleanup
# ==============================================================================

echo -e "\n--- 4. Cleanup ---"
rm -r tmp                   # Remove KofamScan temporary directory
rm -rf ${KOFAM_PREP_DIR}    # Remove the directory with individual prepared files
rm -f ${MAPPING_FILE}       # Remove the temporary mapping file

echo "Full pipeline completed. Results are in ${OUT}"
