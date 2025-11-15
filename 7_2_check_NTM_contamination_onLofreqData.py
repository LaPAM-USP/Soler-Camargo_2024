#!/usr/bin/env bash
# ============================================================
# NTM contamination filter based on LoFreq VCF
# Usage: ./ntmFilter.sh <biosample>
# Input:  lofreq/<biosample>/<biosample>_lofreq.vcf.gz
# Output: ntm/<biosample>/<biosample>_ntm_summary.csv
# ============================================================

set -euo pipefail


BIOSAMPLE="${1:-}"
LOFREQ_DIR="lofreq/${BIOSAMPLE}"
OUTPUT_DIR="ntmFilter/${BIOSAMPLE}"
SUMMARY_CSV="${OUTPUT_DIR}/${BIOSAMPLE}_ntm_summary.csv"

REF_CHR="NC_000962.3"
NTM_POS=1472307
CUTOFF=0.20  # 20% threshold

# -------------------- CHECK INPUT --------------------
if [[ -z "$BIOSAMPLE" ]]; then
    echo "Usage: ./ntmFilter.sh <biosample>"
    exit 1
fi

VCF_FILE="${LOFREQ_DIR}/${BIOSAMPLE}_lofreq.vcf.gz"
if [[ ! -f "$VCF_FILE" ]]; then
    echo "[ERROR] LoFreq VCF not found: ${VCF_FILE}"
    exit 1
fi

if [[ ! -f "${VCF_FILE}.csi" && ! -f "${VCF_FILE}.tbi" ]]; then
    echo "[ERROR] Index file (.csi or .tbi) not found for ${VCF_FILE}"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo "[RUN] Checking NTM contamination for biosample: ${BIOSAMPLE}"
echo "[IN]  Input VCF: ${VCF_FILE}"
echo "[OUT] Output directory: ${OUTPUT_DIR}"
echo "---------------------------------------------"

# -------------------- FIND VARIANT --------------------
AF=$(bcftools query -r "${REF_CHR}:${NTM_POS}-${NTM_POS}" -f '%INFO/AF\n' "$VCF_FILE" 2>/dev/null | head -n 1 || true)
GT=$(bcftools query -r "${REF_CHR}:${NTM_POS}-${NTM_POS}" -f '[%GT]\n' "$VCF_FILE" 2>/dev/null | head -n 1 || true)

if [[ -z "$AF" ]]; then
    echo "[INFO] No variant found at ${REF_CHR}:${NTM_POS}. Assuming AF=0 (no NTM)."
    AF="0"
    GT="0/0"
fi

# Remove any trailing characters and ensure numeric
AF_NUM=$(printf "%.4f" "$AF" 2>/dev/null || echo "0.0000")

# Determine PASS/FAIL
status="PASS"
if (( $(echo "$AF_NUM >= $CUTOFF" | bc -l) )); then
    status="FAIL"
fi

# -------------------- WRITE SUMMARY --------------------
echo "biosample,filename,genotype,position,allele_frequency,status" > "$SUMMARY_CSV"
echo "${BIOSAMPLE},$(basename "$VCF_FILE"),${GT},${NTM_POS},${AF_NUM},${status}" >> "$SUMMARY_CSV"

echo "[OK] Summary CSV generated: ${SUMMARY_CSV}"
echo "---------------------------------------------"
echo "[DONE] NTM filter completed for biosample: ${BIOSAMPLE}"
echo "[OUT] Status: ${status} (AF=${AF_NUM})"
echo "[OUT] Results saved in: ${OUTPUT_DIR}"
