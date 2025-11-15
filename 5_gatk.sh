#!/usr/bin/env bash
# ============================================================
# Variant calling with GATK HaplotypeCaller (GVCF + VCF)
# Usage: ./gatk.sh <biosample>
# Input:  bwa/<biosample>/<biosample>.bam
# Output: gatk/<biosample>/<biosample>.g.vcf.gz and *_gatk.vcf.gz
# Reference: mtbRef/NC0009623.fasta
# ============================================================

set -euo pipefail

BIOSAMPLE="${1:-}"
BWA_DIR="bwa/${BIOSAMPLE}"
REF_DIR="mtbRef"
REF="${REF_DIR}/NC0009623.fasta"
DICT="${REF_DIR}/NC0009623.dict"
OUTPUT_DIR="gatk/${BIOSAMPLE}"
GATK_CMD="gatk"

# -------------------- CHECK INPUT --------------------
if [[ -z "$BIOSAMPLE" ]]; then
    echo "Usage: ./gatk.sh <biosample>"
    exit 1
fi

if [[ ! -d "$BWA_DIR" ]]; then
    echo "[ERROR] BAM directory not found: ${BWA_DIR}"
    exit 1
fi

if [[ ! -f "$REF" ]]; then
    echo "[ERROR] Reference genome not found: ${REF}"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# -------------------- DEPENDENCY CHECKS --------------------
for cmd in "$GATK_CMD" samtools; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "[ERROR] Required command not found: $cmd"
        exit 1
    fi
done

echo "[RUN] GATK variant calling for biosample: ${BIOSAMPLE}"
echo "[IN]  BAM directory: ${BWA_DIR}"
echo "[REF] Reference genome: ${REF}"
echo "[OUT] Output directory: ${OUTPUT_DIR}"
echo "---------------------------------------------"

# -------------------- PREPARE REFERENCE --------------------
if [[ ! -f "${REF}.fai" ]]; then
    echo "[INFO] Creating FASTA index..."
    samtools faidx "$REF"
else
    echo "[INFO] FASTA index already exists: ${REF}.fai"
fi

if [[ ! -f "$DICT" ]]; then
    echo "[INFO] Creating sequence dictionary..."
    $GATK_CMD CreateSequenceDictionary -R "$REF" -O "$DICT"
else
    echo "[INFO] Dictionary already exists: $DICT"
fi

# -------------------- FIND BAM --------------------
BAM_FILE="${BWA_DIR}/${BIOSAMPLE}.bam"

if [[ ! -f "$BAM_FILE" ]]; then
    echo "[ERROR] BAM file not found: ${BAM_FILE}"
    exit 1
fi

echo "[INFO] Using BAM file: ${BAM_FILE}"

# ============================================================
# 1. HAPLOTYPECALLER (GVCF MODE)
# ============================================================
GVCF="${OUTPUT_DIR}/${BIOSAMPLE}.g.vcf.gz"
LOG_GVCF="${OUTPUT_DIR}/${BIOSAMPLE}_gvcf.log"

echo "[RUN] Calling variants in GVCF mode (no MNP)..."
$GATK_CMD HaplotypeCaller \
    -R "$REF" \
    -I "$BAM_FILE" \
    -O "$GVCF" \
    -ERC GVCF \
    -A FisherStrand \
    -A StrandOddsRatio \
    -A QualByDepth \
    -A MappingQualityRankSumTest \
    -A ReadPosRankSumTest \
    -A DepthPerAlleleBySample \
    -A Coverage \
    2>&1 | tee "$LOG_GVCF"

echo "[OK] GVCF generated: ${GVCF}"

# ============================================================
# 2. HAPLOTYPECALLER (VCF MODE – with MNP)
# ============================================================
VCF="${OUTPUT_DIR}/${BIOSAMPLE}_gatk.vcf.gz"
LOG_VCF="${OUTPUT_DIR}/${BIOSAMPLE}_gatk_vcf.log"

echo "[RUN] Calling variants in standard VCF mode (--max-mnp-distance 1)..."
$GATK_CMD HaplotypeCaller \
    -R "$REF" \
    -I "$BAM_FILE" \
    -O "$VCF" \
    --max-mnp-distance 1 \
    -A FisherStrand \
    -A StrandOddsRatio \
    -A QualByDepth \
    -A MappingQualityRankSumTest \
    -A ReadPosRankSumTest \
    -A DepthPerAlleleBySample \
    -A Coverage \
    2>&1 | tee "$LOG_VCF"

echo "[OK] VCF generated: ${VCF}"

# -------------------- COMPLETION --------------------
echo "---------------------------------------------"
echo "[DONE] GATK variant calling completed for biosample: ${BIOSAMPLE}"
echo "[OUT] GVCF: ${GVCF}"
echo "[OUT] VCF:  ${VCF}"
echo "[OUT] Logs: ${OUTPUT_DIR}/"
