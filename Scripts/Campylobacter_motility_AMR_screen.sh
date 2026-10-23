#!/bin/bash -e
#SBATCH --account       massey03212
#SBATCH --job-name      campy_screen
#SBATCH --time          04:00:00
#SBATCH --mem           16GB
#SBATCH --cpus-per-task 8
#SBATCH --error         logs/%x_%j.err
#SBATCH --output        logs/%x_%j.out

# campylobacter_motility_AMR_screen.sh
#
# Screens a directory of genomes (MAGs and/or public assemblies, one FASTA per genome) for
#   1. chemotaxis and flagellar-motor genes by Pfam domains (HMMER hmmsearch --cut_ga on
#      Prodigal proteins, Pfam v32.0): CheA, CheW, FliM, FliG, FliN, flagellin, MotA
#   2. antimicrobial resistance genes with AMRFinderPlus (--plus, nucleotide mode)
#   3. resistance and virulence genes with ABRicate against ncbi, card, resfinder,
#      argannot, megares and vfdb (default thresholds: >=80% identity, >=80% coverage)
#
# Gene calls (per protein, from Pfam domains on the same protein):
#   CheA       Hpt (PF01627) + H-kinase_dim (PF02895) + CheW (PF01584)
#   CheW       PF01584 without PF01627 / PF02895
#   FliM       PF02154            FliG   PF01706            MotA  PF01618
#   FliN       PF01052 without PF02154
#   Flagellin  Flagellin_N (PF00669) + Flagellin_C (PF00700); flaA and flaB are not separated
#   CheY-like  PF00072 (Response_reg) as the only hit, protein <= 160 aa, within 10 genes of a
#              CheA gene. PF00072 is shared by many response regulators and this rule did not
#              recover cheY in Campylobacter, so report cheY from VFDB instead.
#
# Required environment variables:
#   GENOME_DIR     directory with one FASTA per genome (<name>.fa)
#   OUT_DIR        output directory (created)
#   PFAM_DB        path to Pfam-A.hmm (Pfam v32.0)
#   AMRFINDER_DB   path to the AMRFinderPlus database directory
# Optional:
#   GENOME_EXT     genome file extension without dot (default: fa)
#   THREADS        threads when not run under Slurm (default: 4)
#   PYTHON         python3 executable (standard library only; default: python3)
#   MOD_PRODIGAL, MOD_HMMER, MOD_AMRFINDER, MOD_ABRICATE
#                  module names (defaults are NeSI names; ignored if "module" is not available,
#                  in which case the four tools must be on PATH)
#
# Usage (create logs/ first):
#   mkdir -p logs
#   GENOME_DIR=/path/genomes OUT_DIR=/path/out PFAM_DB=/path/Pfam-A.hmm \
#   AMRFINDER_DB=/path/amrfinderplus/latest sbatch campylobacter_motility_AMR_screen.sh
#
# Outputs in OUT_DIR:
#   run_info.txt                      tool and database versions
#   che_flag_by_genome.tsv            gene counts per genome
#   che_flag_by_gene.tsv              one row per called gene
#   amrfinderplus/<genome>.tsv, amrfinderplus_all.tsv
#   abricate/<db>/<genome>.<db>.tsv, abricate_all.tsv

set -eu

: "${GENOME_DIR:?set GENOME_DIR}"
: "${OUT_DIR:?set OUT_DIR}"
: "${PFAM_DB:?set PFAM_DB}"
: "${AMRFINDER_DB:?set AMRFINDER_DB}"
GENOME_EXT=${GENOME_EXT:-fa}
THREADS=${SLURM_CPUS_PER_TASK:-${THREADS:-4}}
PYTHON=${PYTHON:-python3}
MOD_PRODIGAL=${MOD_PRODIGAL:-prodigal/2.6.3-GCCcore-7.4.0}
MOD_HMMER=${MOD_HMMER:-HMMER/3.4-GCC-15.2.0}
MOD_AMRFINDER=${MOD_AMRFINDER:-AMRFinderPlus/4.2.7-foss-2023a}
MOD_ABRICATE=${MOD_ABRICATE:-ABRicate/1.0.0-GCC-11.3.0-Perl-5.34.1}
ACCS="PF01627|PF02895|PF01584|PF00072|PF02154|PF01706|PF01052|PF00669|PF00700|PF01618"

load() {
  if command -v module > /dev/null 2>&1; then
    module purge
    module load "$1"
  fi
}

for f in "$PFAM_DB" "$AMRFINDER_DB"; do
  if [ ! -e "$f" ]; then echo "ERROR: missing: $f" >&2; exit 1; fi
done
shopt -s nullglob
genomes=("$GENOME_DIR"/*."$GENOME_EXT")
shopt -u nullglob
if [ ${#genomes[@]} -eq 0 ]; then echo "ERROR: no .$GENOME_EXT files in $GENOME_DIR" >&2; exit 1; fi
mkdir -p "$OUT_DIR/prodigal" "$OUT_DIR/amrfinderplus" "$OUT_DIR/abricate"
echo "Started: $(date)"
echo "Genomes found: ${#genomes[@]}"

load "$MOD_PRODIGAL"
: > "$OUT_DIR/prodigal/all.faa"
for g in "${genomes[@]}"; do
  name=$(basename "$g" ."$GENOME_EXT")
  prodigal -i "$g" -a "$OUT_DIR/prodigal/${name}.faa" -p meta -q > /dev/null
  awk -v n="$name" '/^>/ { split($1, a, " "); print ">" n "|" substr(a[1], 2); next } { print }' \
      "$OUT_DIR/prodigal/${name}.faa" >> "$OUT_DIR/prodigal/all.faa"
done
echo "Proteins: $(grep -c '^>' "$OUT_DIR/prodigal/all.faa")"

awk -v pat="^(${ACCS})(\\\\.|\$)" '{ rec = rec $0 "\n" }
     /^ACC / { if ($2 ~ pat) keep = 1 }
     /^\/\// { if (keep) printf "%s", rec; rec = ""; keep = 0 }' \
    "$PFAM_DB" > "$OUT_DIR/che_flag.hmm"
echo "Profiles extracted (expected 10):"
grep -E '^(NAME|ACC)' "$OUT_DIR/che_flag.hmm" | paste - -
n_models=$(grep -c '^NAME' "$OUT_DIR/che_flag.hmm" || true)
if [ "$n_models" -ne 10 ]; then echo "ERROR: expected 10 profiles, got $n_models" >&2; exit 1; fi

load "$MOD_HMMER"
hmmsearch --cut_ga --noali -o /dev/null \
    --domtblout "$OUT_DIR/che_flag.domtbl" \
    --cpu "$THREADS" \
    "$OUT_DIR/che_flag.hmm" "$OUT_DIR/prodigal/all.faa"
echo "Domain hits: $(grep -vc '^#' "$OUT_DIR/che_flag.domtbl" || true)"

"$PYTHON" - "$OUT_DIR/che_flag.domtbl" "$OUT_DIR" "$GENOME_DIR" "$GENOME_EXT" << 'PYEOF'
import sys, os, collections, re
domtbl, out, gdir, ext = sys.argv[1:5]
suffix = "." + ext
genomes = sorted(f[:-len(suffix)] for f in os.listdir(gdir) if f.endswith(suffix))
hits = collections.defaultdict(set)
tlen = {}
for line in open(domtbl):
    if line.startswith("#"):
        continue
    p = line.split(None, 22)
    hits[p[0]].add(p[4].split(".")[0])
    tlen[p[0]] = int(p[2])

def classify(s):
    o = []
    if {"PF01627", "PF02895", "PF01584"} <= s:
        o.append("CheA")
    elif "PF01584" in s and not ({"PF01627", "PF02895"} & s):
        o.append("CheW")
    if "PF02154" in s:
        o.append("FliM")
    if "PF01706" in s:
        o.append("FliG")
    if "PF01052" in s and "PF02154" not in s:
        o.append("FliN")
    if {"PF00669", "PF00700"} <= s:
        o.append("Flagellin")
    if "PF01618" in s:
        o.append("MotA")
    return o

calls = {t: classify(s) for t, s in hits.items()}
calls = {t: c for t, c in calls.items() if c}
for t, s in hits.items():
    if s == {"PF00072"} and tlen[t] <= 160:
        calls.setdefault(t, []).append("Response_reg_short")

def parts(t):
    g, rest = t.split("|", 1)
    m = re.match(r"^(.*)_(\d+)$", rest)
    return g, (m.group(1), int(m.group(2))) if m else (rest, 0)

chea = collections.defaultdict(list)
for t, c in calls.items():
    if "CheA" in c:
        g, (ctg, i) = parts(t)
        chea[(g, ctg)].append(i)
order = ["CheA", "CheW", "CheY-like", "FliM", "FliG", "FliN", "Flagellin", "MotA"]
summ = {g: collections.Counter() for g in genomes}
with open(os.path.join(out, "che_flag_by_gene.tsv"), "w") as fh:
    fh.write("gene\tgenome\tclasses\tnear_CheA_within_10_genes\n")
    for t, c in sorted(calls.items()):
        g, (ctg, i) = parts(t)
        near = int(any(abs(i - j) <= 10 for j in chea.get((g, ctg), [])))
        cl = [x for x in c if x != "Response_reg_short"]
        if "Response_reg_short" in c:
            cl.append("CheY-like" if near else "short_response_regulator")
        fh.write(f"{t}\t{g}\t{','.join(cl)}\t{near}\n")
        for x in cl:
            summ[g][x] += 1
with open(os.path.join(out, "che_flag_by_genome.tsv"), "w") as fh:
    fh.write("genome\t" + "\t".join(order) + "\tshort_response_regulator\tflagellar_core\n")
    for g in genomes:
        c = summ[g]
        flag = int(c["FliM"] > 0 and c["FliG"] > 0 and c["Flagellin"] > 0 and c["MotA"] > 0)
        fh.write(g + "\t" + "\t".join(str(c[x]) for x in order) + f"\t{c['short_response_regulator']}\t{flag}\n")
print(f"{len(genomes)} genomes | CheA: {sum(1 for g in genomes if summ[g]['CheA'])} | "
      f"FliM: {sum(1 for g in genomes if summ[g]['FliM'])} | FliN: {sum(1 for g in genomes if summ[g]['FliN'])}")
PYEOF

load "$MOD_AMRFINDER"
first=1
for g in "${genomes[@]}"; do
  name=$(basename "$g" ."$GENOME_EXT")
  amrfinder -n "$g" --plus --database "$AMRFINDER_DB" --threads "$THREADS" \
      -o "$OUT_DIR/amrfinderplus/${name}.tsv"
  if [ $first -eq 1 ]; then
    { printf 'genome\t'; head -n 1 "$OUT_DIR/amrfinderplus/${name}.tsv"; } > "$OUT_DIR/amrfinderplus_all.tsv"
    first=0
  fi
  tail -n +2 "$OUT_DIR/amrfinderplus/${name}.tsv" | awk -v n="$name" '{ print n "\t" $0 }' >> "$OUT_DIR/amrfinderplus_all.tsv"
done
{
  echo "AMRFinderPlus:"
  amrfinder -V 2>&1 || true
} > "$OUT_DIR/run_info.txt"
echo "AMRFinderPlus done: $(date)"

load "$MOD_ABRICATE"
printf 'genome\tdatabase\tSEQUENCE\tSTART\tEND\tSTRAND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t%%COVERAGE\t%%IDENTITY\tDATABASE\tACCESSION\tPRODUCT\tRESISTANCE\n' > "$OUT_DIR/abricate_all.tsv"
for db in ncbi card resfinder argannot megares vfdb; do
  mkdir -p "$OUT_DIR/abricate/$db"
  for g in "${genomes[@]}"; do
    name=$(basename "$g" ."$GENOME_EXT")
    abricate --db "$db" --threads "$THREADS" "$g" > "$OUT_DIR/abricate/$db/${name}.${db}.tsv"
    tail -n +2 "$OUT_DIR/abricate/$db/${name}.${db}.tsv" | cut -f2- | awk -v n="$name" -v d="$db" '{ print n "\t" d "\t" $0 }' >> "$OUT_DIR/abricate_all.tsv"
  done
done
{
  echo "ABRicate:"
  abricate --version 2>&1 || true
  abricate --list 2>&1 || true
} >> "$OUT_DIR/run_info.txt"
{
  echo "Prodigal / HMMER / Pfam:"
  echo "Pfam database: $PFAM_DB"
  echo "Genomes screened: ${#genomes[@]}"
} >> "$OUT_DIR/run_info.txt"

echo "Finished: $(date)"