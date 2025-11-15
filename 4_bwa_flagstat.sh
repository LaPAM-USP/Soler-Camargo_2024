#!/usr/bin/env bash
# ============================================================
# Mapping paired-end reads using BWA-MEM and generating BAMs
# Usage: ./bwa.sh <biosample>
# Input:  trimmomatic/<biosample>/*.fastq.gz
# Output: bwa/<biosample>/*.bam and summary CSV
# Reference: mtbRef/NC0009623.fasta
# ============================================================

set -euo pipefail

BIOSAMPLE="${1:-}"
TRIM_DIR="trimmomatic/${BIOSAMPLE}"
REF_DIR="mtbRef"
REF="${REF_DIR}/NC0009623.fasta"
OUTPUT_DIR="bwa/${BIOSAMPLE}"
SUMMARY_CSV="${OUTPUT_DIR}/${BIOSAMPLE}_bwa_summary.csv"

MIN_MAPPED=95       # Minimum % mapped reads to PASS
MIN_COVERAGE=90     # Minimum % genome coverage to PASS

# -------------------- CHECK INPUT --------------------
if [[ -z "$BIOSAMPLE" ]]; then
    echo "Usage: ./4_bwa.sh <biosample>"
    exit 1
fi

if [[ ! -d "$TRIM_DIR" ]]; then
    echo "[ERROR] Trimmed reads directory not found: ${TRIM_DIR}"
    exit 1
fi

if [[ ! -f "$REF" ]]; then
    echo "[ERROR] Reference genome not found: ${REF}"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# -------------------- DEPENDENCY CHECKS --------------------
for cmd in bwa samtools java; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "[ERROR] Required command not found: $cmd"
        exit 1
    fi
done

if [[ ! -f "picard.jar" ]]; then
    echo "[ERROR] picard.jar not found in current directory."
    exit 1
fi

echo "[RUN] Running BWA-MEM mapping for biosample: ${BIOSAMPLE}"
echo "[IN]  Trimmed reads: ${TRIM_DIR}"
echo "[REF] Reference genome: ${REF}"
echo "[OUT] Output directory: ${OUTPUT_DIR}"
echo "---------------------------------------------"

# -------------------- BUILD INDEX (only once) --------------------
INDEX_FILES=("${REF}.bwt" "${REF}.pac" "${REF}.ann" "${REF}.amb" "${REF}.sa")
if [[ ! -f "${INDEX_FILES[0]}" ]]; then
    echo "[INFO] BWA index not found. Building index for ${REF}..."
    bwa index "$REF"
else
    echo "[INFO] BWA index already present. Skipping index build."
fi

# -------------------- HEADER OF CSV --------------------
echo "biosample,filename,total_reads,mapped_reads,mapped_pct,duplicates_pct,properly_paired_pct,coverage_pct,status" > "$SUMMARY_CSV"

# -------------------- MAPPING LOOP --------------------
for R1_FILE in "${TRIM_DIR}"/*_R1_001.fastq.gz; do
    [[ ! -f "$R1_FILE" ]] && continue
    R2_FILE="${R1_FILE/_R1_/_R2_}"
    if [[ ! -f "$R2_FILE" ]]; then
        echo "[WARN] Missing R2 for ${R1_FILE}. Skipping."
        continue
    fi

    BASENAME=$(basename "$R1_FILE" | sed 's/_R1_001.fastq.gz//')
    SAMPLE_NAME="${BIOSAMPLE}_${BASENAME}"
    TMP_DIR="${OUTPUT_DIR}/tmp"
    mkdir -p "$TMP_DIR"

    echo "[RUN] Mapping ${SAMPLE_NAME}..."

    RAW_BAM="${TMP_DIR}/${SAMPLE_NAME}.bam"
    SORTED_BAM="${TMP_DIR}/${SAMPLE_NAME}_sorted.bam"
    NAMED_BAM="${TMP_DIR}/${SAMPLE_NAME}_named.bam"
    MARKED_BAM="${OUTPUT_DIR}/${SAMPLE_NAME}_dupMarked.bam"
    METRICS_FILE="${TMP_DIR}/${SAMPLE_NAME}_dupMetrics.txt"
    FLAGSTAT_FILE="${TMP_DIR}/${SAMPLE_NAME}_flagstat.txt"

    # ----- Mapping -----
    bwa mem -t 16 "$REF" "$R1_FILE" "$R2_FILE" | samtools view -bS - > "$RAW_BAM"

    # ----- Sorting -----
    samtools sort "$RAW_BAM" -o "$SORTED_BAM"
    rm -f "$RAW_BAM"

    # ----- Indexing sorted -----
    samtools index "$SORTED_BAM"

    # ----- Add read groups -----
    java -jar picard.jar AddOrReplaceReadGroups \
        -I "$SORTED_BAM" \
        -O "$NAMED_BAM" \
        -RGID 1 -RGLB "${BIOSAMPLE}" -RGPL ILLUMINA -RGPU unit1 -RGSM "$SAMPLE_NAME" \
        -VALIDATION_STRINGENCY LENIENT

    # ----- Mark duplicates -----
    java -jar picard.jar MarkDuplicates \
        -I "$NAMED_BAM" \
        -O "$MARKED_BAM" \
        -M "$METRICS_FILE" \
        -VALIDATION_STRINGENCY LENIENT

    # ----- Index final BAM -----
    samtools index "$MARKED_BAM"

    # ----- Flagstat summary -----
    samtools flagstat "$MARKED_BAM" > "$FLAGSTAT_FILE"

    # ----- Extract mapping stats -----
    total_reads=$(grep "in total" "$FLAGSTAT_FILE" | awk '{print $1}')
    mapped_reads=$(grep " mapped (" "$FLAGSTAT_FILE" | awk '{print $1}')
    mapped_pct=$(grep " mapped (" "$FLAGSTAT_FILE" | sed 's/.*(\(.*%\).*/\1/')
    paired_pct=$(grep "properly paired (" "$FLAGSTAT_FILE" | sed 's/.*(\(.*%\).*/\1/')
    dup_pct=$(grep " duplicates" "$FLAGSTAT_FILE" | sed 's/.*(\(.*%\).*/\1/')

    # ----- Compute genome coverage -----
    coverage_pct=$(samtools coverage "$MARKED_BAM" 2>/dev/null | awk 'NR==2 {print $6}')
    [[ -z "$coverage_pct" ]] && coverage_pct="0.00"

    # ----- QC decision -----
    mapped_value=$(echo "$mapped_pct" | tr -d '%')
    coverage_value=$(echo "$coverage_pct" | tr -d '%')

    status="PASS"

    if (( $(echo "$mapped_value < $MIN_MAPPED" | bc -l) )) || \
       (( $(echo "$coverage_value < $MIN_COVERAGE" | bc -l) )); then
        status="FAIL"
    fi

    # ----- Write line to CSV -----
    echo "${BIOSAMPLE},${SAMPLE_NAME},${total_reads},${mapped_reads},${mapped_pct},${dup_pct},${paired_pct},${coverage_pct},${status}" >> "$SUMMARY_CSV"

    echo "[OK] Mapping completed for ${SAMPLE_NAME} → ${status}"
    echo "---------------------------------------------"
done

# -------------------- CLEAN TMP FILES --------------------
echo "[INFO] Cleaning up temporary files..."
rm -rf "${OUTPUT_DIR}/tmp"
echo "[INFO] Temporary folder removed: ${OUTPUT_DIR}/tmp"

# -------------------- COMPLETION --------------------
echo "[DONE] BWA mapping finished for biosample: ${BIOSAMPLE}"
echo "[OUT] Summary CSV: ${SUMMARY_CSV}"
echo "[OUT] BAM files: ${OUTPUT_DIR}/"
