#!/usr/bin/env bash
set -euo pipefail

# Usage:
# bash scripts/human_filter_kraken.sh \
#   READ1 HUMAN_FASTA KRAKEN_DB OUTPUT_DIR [THREADS] [MEMORY_GB] [--read2 READ2]
#
# Without --read2: single-end RNA-seq.
# With --read2: paired-end RNA-seq (both mates processed together).
# All input paths and the output directory are provided by the caller.

usage() {
    printf 'Usage: bash %s READ1 HUMAN_FASTA KRAKEN_DB OUTPUT_DIR [THREADS] [MEMORY_GB] [--read2 READ2]\n' "$0"
}

if [[ ${1:-} == -h || ${1:-} == --help ]]; then
    usage
    exit 0
fi
if (( $# < 4 )); then
    usage >&2
    exit 1
fi

READ1=$1
HUMAN=$2
DB=$3
OUT=$4
shift 4

THREADS=4
MEMORY_GB=8
READ2=''

# Retain compatibility with the original positional THREADS/MEMORY_GB.
if (( $# > 0 )) && [[ $1 != --* ]]; then
    THREADS=$1
    shift
fi
if (( $# > 0 )) && [[ $1 != --* ]]; then
    MEMORY_GB=$1
    shift
fi
while (( $# > 0 )); do
    case "$1" in
        --read2)
            if (( $# < 2 )) || [[ -z $2 || $2 == --* ]]; then
                echo 'ERROR: --read2 requires a FASTQ path.' >&2
                exit 1
            fi
            [[ -z $READ2 ]] || { echo 'ERROR: --read2 specified twice.' >&2; exit 1; }
            READ2=$2
            shift 2
            ;;
        *)
            printf 'ERROR: Unknown argument: %s\n' "$1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

# Validate inputs and resources
[[ -s "$READ1" ]] || { echo "Read 1 FASTQ missing or empty: $READ1" >&2; exit 1; }
if [[ -n "$READ2" ]]; then
    [[ -s "$READ2" ]] || { echo "Read 2 FASTQ missing or empty: $READ2" >&2; exit 1; }
    [[ "$READ1" != "$READ2" ]] || { echo 'Read 1 and Read 2 must be different files.' >&2; exit 1; }
fi
[[ -s "$HUMAN" ]] || { echo "Human FASTA missing or empty: $HUMAN" >&2; exit 1; }
[[ -d "$DB" ]] || { echo "Kraken database not found: $DB" >&2; exit 1; }

[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || {
    echo "THREADS must be a positive integer"
    exit 1
}
[[ "$MEMORY_GB" =~ ^[1-9][0-9]*$ ]] || {
    echo "MEMORY_GB must be a positive integer"
    exit 1
}
for program in hisat2 hisat2-build rnaspades.py gzip awk; do
    command -v "$program" >/dev/null 2>&1 || {
        echo "Required program not found: $program" >&2
        exit 1
    }
done

if command -v kraken2 >/dev/null 2>&1; then
    KRAKEN_CMD=(kraken2)
elif command -v micromamba >/dev/null 2>&1; then
    KRAKEN_CMD=(micromamba run -n kraken kraken2)
else
    echo 'Kraken2 not available directly or via micromamba.' >&2
    exit 1
fi

# Protect existing results
if [[ -e "$OUT/unmapped.fastq.gz" ||
      -e "$OUT/unmapped_1.fastq.gz" ||
      -e "$OUT/unmapped_2.fastq.gz" ||
      -e "$OUT/hisat2.log" ||
      -e "$OUT/kraken.tsv" ||
      -e "$OUT/kraken.report" ||
      -e "$OUT/assembly" ]]; then
    echo "Output files already exist in $OUT; choose a new output directory." >&2
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
    echo "[1/4] Reusing HISAT2 index: $INDEX"
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

if [[ -n "$READ2" ]]; then
    echo '  Mode: paired-end'
    hisat2 \
        -x "$INDEX" \
        -1 "$READ1" -2 "$READ2" \
        -p "$THREADS" \
        --un-conc-gz "$OUT/unmapped_%.fastq.gz" \
        -S /dev/null \
        2> "$OUT/hisat2.log"
    UNMAPPED1="$OUT/unmapped_1.fastq.gz"
    UNMAPPED2="$OUT/unmapped_2.fastq.gz"
else
    echo '  Mode: single-end'
    hisat2 \
        -x "$INDEX" \
        -U "$READ1" \
        -p "$THREADS" \
        --un-gz "$OUT/unmapped.fastq.gz" \
        -S /dev/null \
        2> "$OUT/hisat2.log"
    UNMAPPED1="$OUT/unmapped.fastq.gz"
fi

# A gzip file may be nonempty even if it contains zero FASTQ records.
has_fastq_reads() {
    [[ -f "$1" ]] && gzip -cd "$1" | awk 'END { exit (NR < 4 || NR % 4 != 0) }'
}
if ! has_fastq_reads "$UNMAPPED1"; then
    echo "No valid unmapped FASTQ records in $UNMAPPED1; stopping before assembly." >&2
    exit 0
fi
if [[ -n "$READ2" ]] && ! has_fastq_reads "$UNMAPPED2"; then
    echo "No valid unmapped FASTQ records in $UNMAPPED2; stopping before assembly." >&2
    exit 0
fi

# 3. Assemble retained RNA-seq reads.
echo "[3/4] Assembling RNA-seq reads with rnaSPAdes..."

if [[ -n "$READ2" ]]; then
    rnaspades.py \
        -1 "$UNMAPPED1" -2 "$UNMAPPED2" \
        -o "$OUT/assembly" \
        -t "$THREADS" \
        -m "$MEMORY_GB"
else
    rnaspades.py \
        -s "$UNMAPPED1" \
        -o "$OUT/assembly" \
        -t "$THREADS" \
        -m "$MEMORY_GB"
fi

ASSEMBLED="$OUT/assembly/transcripts.fasta"
[[ -s "$ASSEMBLED" ]] || {
    echo "Assembly did not produce $ASSEMBLED" >&2
    exit 1
}

# 4. Classify assembled transcripts (not the raw FASTQ reads).
echo "[4/4] Classifying assembled transcripts with Kraken2..."
"${KRAKEN_CMD[@]}" \
    --db "$DB" \
    --threads "$THREADS" \
    --report "$OUT/kraken.report" \
    --output "$OUT/kraken.tsv" \
    "$ASSEMBLED"

echo "Completed: $OUT"
echo "Assembled unmapped transcripts: $ASSEMBLED"



