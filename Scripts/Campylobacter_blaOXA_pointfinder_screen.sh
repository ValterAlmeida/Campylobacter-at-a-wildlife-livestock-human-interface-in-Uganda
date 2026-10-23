#!/bin/bash -e
#SBATCH --account       massey03212
#SBATCH --job-name      campy_oxa_pf
#SBATCH --time          06:00:00
#SBATCH --mem           8GB
#SBATCH --cpus-per-task 4
#SBATCH --error         logs/%x_%j.err
#SBATCH --output        logs/%x_%j.out

# campylobacter_blaOXA_pointfinder_screen.sh
#
# Resistance checks on a directory of genomes (one FASTA per genome):
#   1. ABRicate resfinder at default thresholds (>=80% identity and coverage); genomes with
#      no blaOXA hit are listed as negatives
#   2. relaxed ABRicate resfinder re-screen (50% identity, 50% coverage) on all genomes
#   3. targeted BLASTN of reference blaOXA alleles (taken from ABRicate's resfinder database)
#      against the negative genomes; each best hit is classed as full length, partial at a
#      contig end (possible assembly truncation) or partial internal
#   4. optional: BLASTN of blaOXA sequences cut from your own genomes (OWN_REFS table)
#   5. optional: PointFinder via ResFinder (genus level, --unknown_mut) on every genome
#
# Required environment variables:
#   GENOME_DIR     directory with one FASTA per genome (<name>.fa)
#   OUT_DIR        output directory (created)
# Optional:
#   GENOME_EXT     genome file extension without dot (default: fa)
#   THREADS        threads when not run under Slurm (default: 4)
#   PYTHON         python3 executable for the summary step (default: python3)
#   OXA_ALLELES    space-separated ResFinder allele names
#                  (default: "blaOXA-471_1 blaOXA-193_1 blaOXA-489_1")
#   RESFINDER_SEQS path to ABRicate's resfinder "sequences" file (default: found from the abricate install)
#   BLAST_EVALUE   BLASTN e-value (default: 1e-5)
#   END_DIST       bp from a contig end that counts as "at contig end" (default: 100)
#   OWN_REFS       tab-separated table, columns: genome file name, region, strand (+ or -), label
#                  (see own_mag_oxa_regions.example.tsv); regions are cut from files in GENOME_DIR
#   RESFINDER_DB, POINTFINDER_DB
#                  git clones of the CGE resfinder_db and pointfinder_db; PointFinder runs only if both are set
#   RESFINDER_PYTHON  python with "resfinder" installed (default: python3)
#   MOD_ABRICATE, MOD_BLAST, MOD_SAMTOOLS
#                  module names (NeSI defaults; ignored if "module" is not available)
#
# Usage (create logs/ first):
#   mkdir -p logs
#   GENOME_DIR=/path/genomes OUT_DIR=/path/out sbatch campylobacter_blaOXA_pointfinder_screen.sh
#
# Outputs in OUT_DIR:
#   abricate_resfinder_default.tsv, abricate_resfinder_relaxed.tsv, negative_genomes.txt
#   blast/<genome>.<refset>.tsv, blaOXA_blast_summary.tsv
#   pointfinder/<genome>/            (if enabled)
#   run_info.txt                     tool versions and database commit hashes

set -eu

: "${GENOME_DIR:?set GENOME_DIR}"
: "${OUT_DIR:?set OUT_DIR}"
GENOME_EXT=${GENOME_EXT:-fa}
THREADS=${SLURM_CPUS_PER_TASK:-${THREADS:-4}}
PYTHON=${PYTHON:-python3}
RESFINDER_PYTHON=${RESFINDER_PYTHON:-python3}
OXA_ALLELES=${OXA_ALLELES:-"blaOXA-471_1 blaOXA-193_1 blaOXA-489_1"}
BLAST_EVALUE=${BLAST_EVALUE:-1e-5}
END_DIST=${END_DIST:-100}
MOD_ABRICATE=${MOD_ABRICATE:-ABRicate/1.0.0-GCC-11.3.0-Perl-5.34.1}
MOD_BLAST=${MOD_BLAST:-BLAST/2.16.0-GCC-12.3.0}
MOD_SAMTOOLS=${MOD_SAMTOOLS:-SAMtools/1.23.1-GCC-12.3.0}

load() {
  if command -v module > /dev/null 2>&1; then
    module purge
    for m in "$@"; do module load "$m"; done
  fi
}

shopt -s nullglob
genomes=("$GENOME_DIR"/*."$GENOME_EXT")
shopt -u nullglob
if [ ${#genomes[@]} -eq 0 ]; then echo "ERROR: no .$GENOME_EXT files in $GENOME_DIR" >&2; exit 1; fi
mkdir -p "$OUT_DIR/blast"
echo "Started: $(date)"
echo "Genomes found: ${#genomes[@]}"

load "$MOD_ABRICATE" "$MOD_BLAST" "$MOD_SAMTOOLS"

abricate --db resfinder "${genomes[@]}" > "$OUT_DIR/abricate_resfinder_default.tsv"
abricate --db resfinder --minid 50 --mincov 50 "${genomes[@]}" > "$OUT_DIR/abricate_resfinder_relaxed.tsv"

awk -F'\t' 'NR > 1 && $6 ~ /^blaOXA/ { print $1 }' "$OUT_DIR/abricate_resfinder_default.tsv" | sort -u > "$OUT_DIR/.oxa_positive.txt"
: > "$OUT_DIR/negative_genomes.txt"
for g in "${genomes[@]}"; do
  if ! grep -qxF "$g" "$OUT_DIR/.oxa_positive.txt"; then echo "$g" >> "$OUT_DIR/negative_genomes.txt"; fi
done
rm -f "$OUT_DIR/.oxa_positive.txt"
n_neg=$(wc -l < "$OUT_DIR/negative_genomes.txt")
echo "Genomes with a blaOXA hit at default thresholds: $(( ${#genomes[@]} - n_neg )); negatives: $n_neg"

if [ -n "${RESFINDER_SEQS:-}" ]; then
  DB_SEQS=$RESFINDER_SEQS
else
  ABRICATE_ROOT=$(dirname "$(readlink -f "$(command -v abricate)")")/..
  DB_SEQS=$ABRICATE_ROOT/db/resfinder/sequences
fi
if [ ! -s "$DB_SEQS" ]; then echo "ERROR: ResFinder sequences not found: $DB_SEQS (set RESFINDER_SEQS)" >&2; exit 1; fi
OXA_REFS=$OUT_DIR/oxa_refs.fasta
: > "$OXA_REFS"
for gene in $OXA_ALLELES; do
  awk -v g="$gene" 'BEGIN { p = 0 } /^>/ { p = ($0 ~ (g "(~| |$)")) } p' "$DB_SEQS" >> "$OXA_REFS"
done
n_refs=$(grep -c '^>' "$OXA_REFS" || true)
n_want=$(echo $OXA_ALLELES | wc -w)
if [ "$n_refs" -ne "$n_want" ]; then echo "ERROR: extracted $n_refs reference alleles, expected $n_want" >&2; exit 1; fi

REFSETS="oxa_refs"
if [ -n "${OWN_REFS:-}" ]; then
  OWN_FASTA=$OUT_DIR/own_MAG_oxa_refs.fasta
  : > "$OWN_FASTA"
  while IFS=$'\t' read -r gfile region strand label; do
    if [ -z "$gfile" ] || [ "${gfile:0:1}" = "#" ]; then continue; fi
    if [ "$strand" = "-" ]; then rc="-i"; else rc=""; fi
    samtools faidx $rc "$GENOME_DIR/$gfile" "$region" | sed "1s/^>.*/>${label}/" >> "$OWN_FASTA"
  done < "$OWN_REFS"
  REFSETS="oxa_refs own_refs"
fi

while read -r g; do
  [ -z "$g" ] && continue
  name=$(basename "$g" ."$GENOME_EXT")
  for set in $REFSETS; do
    if [ "$set" = "oxa_refs" ]; then q=$OXA_REFS; else q=$OWN_FASTA; fi
    blastn -query "$q" -subject "$g" -evalue "$BLAST_EVALUE" \
      -outfmt "6 qseqid sseqid pident length qlen slen sstart send evalue bitscore" \
      > "$OUT_DIR/blast/${name}.${set}.tsv"
  done
done < "$OUT_DIR/negative_genomes.txt"

"$PYTHON" - "$OUT_DIR" "$END_DIST" "$GENOME_EXT" << 'PYEOF'
import sys, os, glob, collections
out, end_dist, ext = sys.argv[1], int(sys.argv[2]), sys.argv[3]
rows = []
def short(q):
    return q.split("~~~")[1] if "~~~" in q else q
for path in sorted(glob.glob(os.path.join(out, "blast", "*.tsv"))):
    base = os.path.basename(path)[:-4]
    genome, refset = base.rsplit(".", 1)
    best = {}
    seen_q = set()
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        q, s, pid, aln, qlen, slen, ss, se, ev, bit = p
        seen_q.add(q)
        bit = float(bit)
        if q not in best or bit > best[q][0]:
            best[q] = (bit, s, float(pid), int(aln), int(qlen), int(slen), int(ss), int(se), ev)
    refs = []
    for line in open(os.path.join(out, "oxa_refs.fasta" if refset == "oxa_refs" else "own_MAG_oxa_refs.fasta")):
        if line.startswith(">"):
            refs.append(line[1:].split()[0])
    for q in refs:
        if q not in best:
            rows.append((genome, refset, short(q), "no_hit", "", "", "", "", ""))
            continue
        bit, s, pid, aln, qlen, slen, ss, se, ev = best[q]
        qcov = aln / qlen
        dist = min(min(ss, se) - 1, slen - max(ss, se))
        if qcov >= 0.9:
            cls = "full_length"
        elif dist <= end_dist:
            cls = "partial_at_contig_end"
        else:
            cls = "partial_internal"
        rows.append((genome, refset, short(q), cls, s, f"{pid:.2f}", f"{qcov:.2f}", str(dist), ev))
with open(os.path.join(out, "blaOXA_blast_summary.tsv"), "w") as fh:
    fh.write("genome\treference_set\tquery\tclass\tcontig\tpident\tquery_coverage\tdist_to_contig_end\tevalue\n")
    for r in rows:
        fh.write("\t".join(r) + "\n")
counts = collections.Counter(r[3] for r in rows)
print("BLAST classes:", dict(counts))
PYEOF

{
  echo "ABRicate:"
  abricate --version 2>&1 || true
  abricate --list 2>&1 || true
  echo "BLAST:"
  blastn -version 2>&1 | head -n 1 || true
  echo "SAMtools:"
  samtools --version 2>&1 | head -n 1 || true
} > "$OUT_DIR/run_info.txt"

if [ -n "${RESFINDER_DB:-}" ] && [ -n "${POINTFINDER_DB:-}" ]; then
  mkdir -p "$OUT_DIR/pointfinder"
  for g in "${genomes[@]}"; do
    name=$(basename "$g" ."$GENOME_EXT")
    "$RESFINDER_PYTHON" -m resfinder \
      -ifa "$g" \
      -o "$OUT_DIR/pointfinder/$name" \
      -s "Campylobacter" -acq --point -u \
      -db_res "$RESFINDER_DB" -db_point "$POINTFINDER_DB"
  done
  {
    echo "ResFinder:"
    "$RESFINDER_PYTHON" -m resfinder --version 2>&1 || true
    echo "resfinder_db commit: $(git -C "$RESFINDER_DB" rev-parse HEAD 2>/dev/null || echo unknown)"
    echo "pointfinder_db commit: $(git -C "$POINTFINDER_DB" rev-parse HEAD 2>/dev/null || echo unknown)"
  } >> "$OUT_DIR/run_info.txt"
else
  echo "PointFinder skipped (RESFINDER_DB and POINTFINDER_DB not set)"
fi

echo "Finished: $(date)"
