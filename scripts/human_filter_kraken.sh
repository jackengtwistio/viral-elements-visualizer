
#!/usr/bin/env bash
set -euo pipefail

# Usage:
# bash scripts/human_filter_kraken.sh \
#   FASTQ HUMAN_FASTA KRAKEN_DB OUTPUT_DIR [THREADS] [MEMORY_GB]

FASTQ=${1:?Provide FASTQ file}
HUMAN=${2:?Provide human FASTA}
DB=${3:?Provide Kraken2 database}
OUT=${4:?Provide output directory}
THREADS=${5:-4}
MEMORY_GB=${6:-8}

# Validate inputs
[[ -f "$FASTQ" ]] || { echo "FASTQ not found"; exit 1; }
[[ -f "$HUMAN" ]] || { echo "Human FASTA not found"; exit 1; }
[[ -d "$DB" ]] || { echo "Kraken database not found"; exit 1; }

[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || {
    echo "THREADS must be a positive integer"
    exit 1
}

[[ "$MEMORY_GB" =~ ^[1-9][0-9]*$ ]] || {
    echo "MEMORY_GB must be a positive integer"
    exit 1
}

command -v rnaspades.py >/dev/null 2>&1 || {
    echo "rnaSPAdes not found: install SPAdes or add rnaspades.py to PATH"
    exit 1
}

# Protect existing results
if [[ -e "$OUT/unmapped.fastq.gz" ||
      -e "$OUT/kraken.tsv" ||
      -e "$OUT/kraken.report" ]]; then
    echo "Output files already exist in $OUT"
    exit 1
fi

# 1. Build or reuse HISAT2 index
INDEX="${HUMAN%.*}"

index_complete() {
    local ext i ok

    for ext in ht2 ht2l; do
        ok=1
        for i in {1..8}; do
            [[ -s "${INDEX}.${i}.${ext}" ]] || {
                ok=0
                break
            }
        done
        (( ok )) && return 0
    done

    return 1
}

if index_complete; then
    echo "[1/3] Reusing HISAT2 index: $INDEX"
else
    if compgen -G "${INDEX}.*.ht2*" > /dev/null; then
        echo "Incomplete HISAT2 index found: $INDEX"
        echo "Stopping to avoid overwriting existing index files."
        exit 1
    fi

    echo "[1/4] Building HISAT2 index..."
    hisat2-build "$HUMAN" "$INDEX"
fi

mkdir -p "$OUT"

# 2. Align to human and retain unmapped reads
echo "[2/4] Filtering human reads..."

hisat2 \
    -x "$INDEX" \
    -U "$FASTQ" \
    -p "$THREADS" \
    --un-gz "$OUT/unmapped.fastq.gz" \
    -S /dev/null \
    2> "$OUT/hisat2.log"



# 3. De novo assembly of all human-unmapped RNA-seq reads
echo "[3/4] Assembling unmapped RNA-seq reads with rnaSPAdes..."

rnaspades.py \
    -s "$OUT/unmapped.fastq.gz" \
    -o "$OUT/assembly" \
    -t "$THREADS" \
    -m "$MEMORY_GB"

ASSEMBLED="$OUT/assembly/transcripts.fasta"

[[ -s "$ASSEMBLED" ]] || {
    echo "Assembly did not produce $ASSEMBLED" >&2
    exit 1
}

# 4. Classify unmapped reads against viral database
echo "[4/4] Running Kraken2..."

micromamba run -n kraken kraken2 \
    --db "$DB" \
    --threads "$THREADS" \
    --gzip-compressed \
    --report "$OUT/kraken.report" \
    --output "$OUT/kraken.tsv" \
    "$ASSEMBLED/transcripts.fasta"

echo "Completed: $OUT"
echo "Assembled unmapped transcripts: $ASSEMBLED"
