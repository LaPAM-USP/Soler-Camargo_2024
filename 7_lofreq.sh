#!/usr/bin/env bash
# ============================================================
# Variant calling with LoFreq (overwrite enabled)
# Usage: ./lofreq.sh <biosample>
# Input:  bwa/<biosample>/<biosample>.bam
# Output: lofreq/<biosample>/*.vcf.gz
# Reference: mtbRef/NC0009623.fasta
# ============================================================

set -euo pipefail

BIOSAMPLE="${1:-}"
BWA_DIR="bwa/${BIOSAMPLE}"
REF_DIR="mtbRef"
REF="${REF_DIR}/NC0009623.fasta"
OUTPUT_DIR="lofreq/${BIOSAMPLE}"

# -------------------- CHECK INPUT --------------------
if [[ -z "$BIOSAMPLE" ]]; then
    echo "Usage: ./lofreq.sh <biosample>"
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
for cmd in lofreq samtools bgzip bcftools; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "[ERROR] Required command not found: $cmd"
        exit 1
    fi
done

echo "[RUN] Running LoFreq for biosample: ${BIOSAMPLE}"
echo "[IN]  BAM directory: ${BWA_DIR}"
echo "[REF] Reference genome: ${REF}"
echo "[OUT] Output directory: ${OUTPUT_DIR}"
echo "---------------------------------------------"

# -------------------- FIND BAM FILE --------------------
# Updated: expect final BAM called <biosample>.bam
BAM_FILE="${BWA_DIR}/${BIOSAMPLE}.bam"

if [[ ! -f "$BAM_FILE" ]]; then
    echo "[ERROR] BAM file not found: ${BAM_FILE}"
    exit 1
fi

echo "[INFO] Using BAM file: ${BAM_FILE}"

# -------------------- RUN LoFreq (PIPE OPTIMIZED) --------------------
VCF_OUTPUT="${OUTPUT_DIR}/${BIOSAMPLE}_lofreq.vcf"

echo "[RUN] Adding indel + alignment qualities and calling variants..."
lofreq indelqual --dindel -f "$REF" -o - "$BAM_FILE" | \
lofreq alnqual -b -r - "$REF" | \
lofreq call --force-overwrite -f "$REF" -o "$VCF_OUTPUT" --call-indels -
echo "[OK] Variant calling complete: $(basename "$VCF_OUTPUT")"

# -------------------- COMPRESS & INDEX VCF --------------------
echo "[RUN] Compressing and indexing VCF..."
bgzip -f "$VCF_OUTPUT"
bcftools index -f "${VCF_OUTPUT}.gz"
echo "[OK] Compressed and indexed VCF: $(basename "${VCF_OUTPUT}.gz")"

# -------------------- CLEANUP --------------------
echo "[CLEANUP] Removing temporary BAMs (if any)..."
find "$OUTPUT_DIR" -type f \( -name "*_lofreq.bam" -o -name "*_aln_qual.bam" \) -delete
echo "[OK] Temporary files removed."

# -------------------- COMPLETION --------------------
echo "[OUT] Final VCF: ${VCF_OUTPUT}.gz"
echo "[OUT] Index: ${VCF_OUTPUT}.gz.csi"
echo "[OUT] Directory: ${OUTPUT_DIR}"
