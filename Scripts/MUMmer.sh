## MUMmer - Pipeline:

## Step 1:

dnadiff -p C_vicugnae_vs_cjejuni GCA_000009085.1_C_jejuni_genomic.fa maxbin.030.go30.fa

show-coords -rclTH C_vicugnae_vs_cjejuni.delta > C_vicugnae_vs_cjejuni.coords.tab

delta-filter -1 C_vicugnae_vs_cjejuni.delta > C_vicugnae_vs_cjejuni.filtered.delta

show-coords -rclTH C_vicugnae_vs_cjejuni.filtered.delta > filtered.coords.tab


## Step 2:

## Longer mapped regions:

#!/bin/bash

DELTA_FILE="C_vicugnae_vs_cjejuni.filtered.delta"

grep "^>" "$DELTA_FILE" | cut -d' ' -f1,2 | sed 's/^>//' | while read ref query; do
    echo "Processing $ref vs $query"
    outname="${ref}_vs_${query}.txt"
    show-aligns "$DELTA_FILE" "$ref" "$query" > "$outname"
done

## To combine mapped sequences:

mkdir -p combined_sequences

## Create a script:

# combine_alignments.awk
BEGIN {
    RS = "-- BEGIN alignment";
    FS = "\n";
    ORS = "";
    combined_ref = "";
    combined_query = "";
}

{
    ref_seq = "";
    query_seq = "";
    for (i = 3; i <= NF; i++) {
        if ($i ~ /^[0-9]/) {
            split($i, parts, / +/);
            ref_line = parts[2];
            i++;
            split($i, parts, / +/);
            query_line = parts[2];
            for (j = 1; j <= length(ref_line); j++) {
                ref_char = substr(ref_line, j, 1);
                query_char = substr(query_line, j, 1);
                query_seq = query_seq (query_char == "." ? ref_char : query_char);
            }
            ref_seq = ref_seq ref_line;
        }
    }
    combined_ref = combined_ref ref_seq;
    combined_query = combined_query query_seq;
}

END {
    print "Combined Reference (" ref_name "):\n" combined_ref "\n";
    print "Combined Query (" query_name "):\n" combined_query "\n";
}

## Run the script on all files:

#!/bin/bash

# Directory containing alignment files
ALIGNMENTS_DIR="combined_sequences"

# Output directory
OUTPUT_DIR="combined_sequences"

# Process each file
for file in "$ALIGNMENTS_DIR"/AL111168.1_vs_*.txt; do
    # Extract reference and query names from the filename
    base=$(basename "$file" .txt)
    ref_name="AL111168.1"
    query_name=$(echo "$base" | cut -d'_' -f4- | sed 's/_vs_//')  # Extracts k141_XXXXXXX

    # Run the AWK script with variables for ref/query names
    awk -v ref_name="$ref_name" -v query_name="$query_name" \
        -f combine_alignments.awk "$file" > "$OUTPUT_DIR/${base}_combined.txt"

    echo "Processed: $file → $OUTPUT_DIR/${base}_combined.txt"
done
