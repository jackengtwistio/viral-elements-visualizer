
#!/usr/bin/env bash
set -euo pipefail

# Use existing Kraken2 results; do not rerun Kraken2.
# Required: --fastq --kraken-tsv --reference --out --taxid
# Optional: --threads (default 4)

usage() {
    echo "Usage: bash $0 --fastq FILE --kraken-tsv FILE --reference FASTA --out DIR --taxid ID [--threads N]"
}

FASTQ='' KRAKEN_TSV='' REF='' OUT='' TAXID='' THREADS=4

while (( $# )); do
    case "$1" in
        --fastq) FASTQ=${2:?Missing FASTQ}; shift 2 ;;
        --kraken-tsv) KRAKEN_TSV=${2:?Missing Kraken TSV}; shift 2 ;;
        --reference) REF=${2:?Missing reference FASTA}; shift 2 ;;
        --out) OUT=${2:?Missing output directory}; shift 2 ;;
        --taxid) TAXID=${2:?Missing taxid}; shift 2 ;;
        --threads) THREADS=${2:?Missing thread count}; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
    esac
done

if [[ -z "$FASTQ" || -z "$KRAKEN_TSV" || -z "$REF" || -z "$OUT" || -z "$TAXID" ]]; then
    usage >&2
    exit 1
fi

[[ -f "$FASTQ" && -f "$KRAKEN_TSV" && -f "$REF" ]] || {
    echo 'Input FASTQ, Kraken TSV, or reference FASTA not found.' >&2
    exit 1
}

[[ "$TAXID" =~ ^[0-9]+$ && "$THREADS" =~ ^[1-9][0-9]*$ ]] || {
    echo 'Taxid must be numeric and threads must be a positive integer.' >&2
    exit 1
}

mkdir -p "$OUT"
IDS="$OUT/$TAXID.ids"
READS="$OUT/$TAXID.fastq"
BAM="$OUT/$TAXID.sorted.bam"

# 1. Extract reads assigned exactly to the selected taxid.
awk -F '\t' -v taxid="$TAXID" \
    '$1=="C" && $3==taxid {print $2}' \
    "$KRAKEN_TSV" > "$IDS"

[[ -s "$IDS" ]] || {
    echo "No reads assigned to taxid $TAXID"
    exit 0
}

echo "Extracting $(wc -l < "$IDS") reads assigned to taxid $TAXID..."

seqtk subseq "$FASTQ" "$IDS" > "$READS"

[[ -s "$READS" ]] || {
    echo 'No matching reads in the input FASTQ.' >&2
    exit 1
}

# 2. Align reads to the viral reference.
minimap2 -t "$THREADS" -ax sr "$REF" "$READS" |
    samtools sort -@ "$THREADS" -o "$BAM" -

samtools index "$BAM"
samtools faidx "$REF"

# 3. Report coverage.
samtools coverage "$BAM"

# 4. Automatically find the position with the highest read depth.
PEAK=$(samtools depth "$BAM" | awk '
    $3 > max_depth {max_depth=$3; location=$1 ":" $2}
    END {if (max_depth > 0) print location "\t" max_depth}
')

if [[ -z "$PEAK" ]]; then
    echo 'No mapped positions found; skipping visualization.'
    exit 0
fi

IFS=$'\t' read -r VIEW_POS MAX_DEPTH <<< "$PEAK"

printf 'Highest-depth position: %s (depth %s)\n' "$VIEW_POS" "$MAX_DEPTH"

# 5. Ask whether to visualize the alignment at that position.
TVIEW=(samtools tview -p "$VIEW_POS" "$BAM" "$REF")

if [[ -t 0 && -t 1 ]]; then
    read -r -p 'Visualize alignment with samtools tview? [y/N] ' answer || answer=''

    case "$answer" in
        y|Y|yes|YES|Yes) "${TVIEW[@]}" ;;
    esac
else
    printf 'To visualize later: '
    printf '%q ' "${TVIEW[@]}"
    printf '\n'
fi
