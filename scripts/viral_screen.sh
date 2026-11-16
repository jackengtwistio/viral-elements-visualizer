
#!/usr/bin/env bash
set -euo pipefail

FASTQ=$1
# SAMPLE=$2
# TAXID=${3:-}

DB="$HOME/EVE_project/refs/viral"
REF="refs/viral/isfahan.fa"
OUT="results"

mkdir -p "$OUT"

# 1. Detect viral candidates (run only once per sample)
if [[ ! -s "$OUT/kraken.tsv" ]]; then
    micromamba run -n kraken kraken2 \
        --db "$DB" --threads 4 \
        --report "$OUT/kraken.report" \
        --output "$OUT/kraken.tsv" \
        "$FASTQ"
fi

# Stop here if no virus was selected
if [[ -z "$TAXID" ]]; then
    echo "Check viral candidates: $OUT/kraken.report"
    exit 0
fi

# 2. Extract reads assigned to the selected virus
awk -F '\t' -v t="$TAXID" \
    '$1=="C" && $3==t {print $2}' \
    "$OUT/kraken.tsv" > "$OUT/$TAXID.ids"

if [[ ! -s "$OUT/$TAXID.ids" ]]; then
    echo "No reads assigned to taxid $TAXID"
    exit 0
fi

seqtk subseq "$FASTQ" "$OUT/$TAXID.ids" \
    > "$OUT/$TAXID.fastq"

# 3. Align reads to viral reference
minimap2 -ax sr "$REF" "$OUT/$TAXID.fastq" |
    samtools sort -o "$OUT/$TAXID.bam" -

samtools index "$OUT/$TAXID.bam"
samtools faidx "$REF"

# 4. Show alignment statistics
samtools coverage "$OUT/$TAXID.bam"
