#!/usr/bin/env bash
# ============================================================
# Run FastQC for all files belonging to a given biosample
# Usage: ./fastqc.sh <biosample ID>
# Input files must use the Illumina read naming convention:
# 1827-22_S1_L001_R1_001.fastq.gz
# ││││││││││││││││││││││││││
# └─── biosample ID (1827-22)
#        └── biosample number (S1)
#           └── lane (L002)
#                └── read direction (R1/R2)
#                   └── file index (_001)
# ============================================================

set -euo pipefail

INPUT_PATH="$1"
BIOSAMPLE=$(basename "$INPUT_PATH")
OUTPUT_DIR="fastqc/${BIOSAMPLE}"
REPORT_CSV="${OUTPUT_DIR}/fastqc_summary_${BIOSAMPLE}.csv"

MIN_QUALITY=20      # minimum mean per-base quality for PASS
GC_MIN=64.5         # lower GC threshold for M. tuberculosis
GC_MAX=66.5         # upper GC threshold for M. tuberculosis

mkdir -p "$OUTPUT_DIR"

>&2 echo "[RUN] Running FastQC for biosample: $BIOSAMPLE"
>&2 echo "[DIR] Searching files in: $(dirname "$INPUT_PATH")"
>&2 echo "[OUT] Output directory: $OUTPUT_DIR"
>&2 echo "---------------------------------------------"

# ======== FIND FILES ========
FILES=$(find "$(dirname "$INPUT_PATH")" -type f \
    \( -name "*${BIOSAMPLE}*.fastq" -o -name "*${BIOSAMPLE}*.fastq.gz" \) \
    | sort)

if [[ -z "$FILES" ]]; then
    >&2 echo "[ERROR] No files found containing '${BIOSAMPLE}' in their names."
    exit 1
fi

FILE_COUNT=$(echo "$FILES" | wc -l)

# ======== CHECK FILE COUNT IS EVEN ========
if (( FILE_COUNT % 2 != 0 )); then
    >&2 echo "[ERROR] Found an odd number of files (${FILE_COUNT}) for biosample '${BIOSAMPLE}'."
    >&2 echo "[HINT] Each biosample must have paired-end reads (e.g., R1 and R2 per lane)."
    exit 1
fi

>&2 echo "[INFO] Files found (${FILE_COUNT}):"
echo "$FILES" | sed 's/^/   • /'
>&2 echo "---------------------------------------------"

# ======== RUN FASTQC ========
for FILE in $FILES; do
    >&2 echo "[RUN] FastQC on: $FILE"
    fastqc -o "$OUTPUT_DIR" "$FILE"
done

>&2 echo "[OK] FastQC completed for all files."
>&2 echo "[INFO] Unzipping FastQC reports..."
for ZIP in "${OUTPUT_DIR}"/*_fastqc.zip; do
    unzip -qo "$ZIP" -d "$OUTPUT_DIR"
done

>&2 echo "[INFO] Generating summary report: ${REPORT_CSV}"
>&2 echo "---------------------------------------------"

echo "biosample,filename,total_sequences,poor_quality_sequences,sequence_length,percent_GC,mean_per_base_quality,status" > "$REPORT_CSV"

for QC_DIR in "${OUTPUT_DIR}"/*_fastqc; do
    DATA_FILE="${QC_DIR}/fastqc_data.txt"
    [[ ! -f "$DATA_FILE" ]] && continue

    filename=$(grep "^Filename" "$DATA_FILE" | cut -f2)
    total=$(grep "^Total Sequences" "$DATA_FILE" | cut -f2)
    poor=$(grep "^Sequences flagged as poor quality" "$DATA_FILE" | cut -f2)
    length=$(grep "^Sequence length" "$DATA_FILE" | cut -f2)
    gc=$(grep "^%GC" "$DATA_FILE" | cut -f2)

    # Mean "Per base sequence quality"
    mean_quality=$(awk '
        /^#Base/{getline; while ($1 ~ /^[0-9]/) {sum+=$2; n++; getline} }
        END{if (n>0) printf "%.2f", sum/n; else print "NA"}' "$DATA_FILE")

    # Determine PASS/FAIL
    status="PASS"

    if [[ "$mean_quality" == "NA" || "$gc" == "" ]]; then
        status="FAIL"
    else
        # Check mean quality
        if (( $(echo "$mean_quality < $MIN_QUALITY" | bc -l) )); then
            status="FAIL"
        fi
        # Check GC range
        if (( $(echo "$gc < $GC_MIN" | bc -l) )) || (( $(echo "$gc > $GC_MAX" | bc -l) )); then
            status="FAIL"
        fi
    fi

    echo "${BIOSAMPLE},${filename},${total},${poor},${length},${gc},${mean_quality},${status}" >> "$REPORT_CSV"
done

>&2 echo "[OK] Summary CSV generated: $REPORT_CSV"
>&2 echo "---------------------------------------------"
>&2 echo "[DONE] Analysis completed for biosample: $BIOSAMPLE"
>&2 echo "[OUT] Reports available in: $OUTPUT_DIR"
