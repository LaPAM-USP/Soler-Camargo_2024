#!/usr/bin/env bash
# ============================================================
# NTM contamination filter (LoFreq + GATK GVCF interval support)
# Usage: ./ntmFilter.sh <biosample>
# ============================================================

set -euo pipefail

BIOSAMPLE="${1:-}"
LOFREQ_DIR="lofreq/${BIOSAMPLE}"
GATK_DIR="gatk/${BIOSAMPLE}"
OUTPUT_DIR="ntmFilter/${BIOSAMPLE}"
SUMMARY_CSV="${OUTPUT_DIR}/${BIOSAMPLE}_ntm_summary.csv"

REF_CHR="NC_000962.3"
NTM_POS=1472307
CUTOFF=0.20   # Only used for LoFreq AF

mkdir -p "$OUTPUT_DIR"

# -------------------- CHECK INPUT --------------------
if [[ -z "$BIOSAMPLE" ]]; then
    echo "Usage: ./ntmFilter.sh <biosample>"
    exit 1
fi

VCF_LOFREQ="${LOFREQ_DIR}/${BIOSAMPLE}_lofreq.vcf.gz"
if [[ ! -f "$VCF_LOFREQ" ]]; then
    echo "[ERROR] LoFreq VCF not found: $VCF_LOFREQ"
    exit 1
fi

GVCF="${GATK_DIR}/${BIOSAMPLE}.g.vcf.gz"
if [[ ! -f "$GVCF" ]]; then
    echo "[ERROR] GATK GVCF not found: $GVCF"
    exit 1
fi

echo "[RUN] Checking NTM contamination for biosample: ${BIOSAMPLE}"
echo "[IN]  LoFreq VCF: ${VCF_LOFREQ}"
echo "[IN]  GATK GVCF:  ${GVCF}"
echo "[OUT] Output:     ${OUTPUT_DIR}"
echo "---------------------------------------------"

# ============================================================
# STEP 1 — CHECK LoFreq VARIANT
# ============================================================

AF=$(bcftools query -r "${REF_CHR}:${NTM_POS}-${NTM_POS}" -f '%INFO/AF\n' "$VCF_LOFREQ" 2>/dev/null | head -n 1 || true)
GT_LF=$(bcftools query -r "${REF_CHR}:${NTM_POS}-${NTM_POS}" -f '[%GT]\n' "$VCF_LOFREQ" 2>/dev/null | head -n 1 || true)

if [[ -n "${AF}" ]]; then
    # LoFreq detected a variant
    AF_NUM=$(printf "%.4f" "$AF")

    if (( $(echo "$AF_NUM >= $CUTOFF" | bc -l) )); then
        STATUS="FAIL"
    else
        STATUS="PASS"
    fi

    DP_NUM="NA"
    GT="$GT_LF"
else
    # ============================================================
    # STEP 2 — NO VARIANT IN LOFREQ → CHECK GVCF BLOCK
    # ============================================================

    GVCF_LINE=$(bcftools view \
        --include "POS <= ${NTM_POS} && INFO/END >= ${NTM_POS}" \
        -H "$GVCF" 2>/dev/null | head -n 1 || true)

    if [[ -z "$GVCF_LINE" ]]; then
        STATUS="NOCOV"
        AF_NUM="0.0000"
        GT="./."
        DP_NUM=0
    else
        # Extract FORMAT field
        FORMAT_FIELD=$(echo "$GVCF_LINE" | awk '{print $10}')

        # Extract GT and DP from the FORMAT
        GT=$(echo "$FORMAT_FIELD" | cut -d: -f1)
        DP_NUM=$(echo "$FORMAT_FIELD" | cut -d: -f2)
        DP_NUM=$(echo "${DP_NUM:-0}" | grep -Eo '^[0-9]+' || echo 0)

        # Interpret GT → AF → STATUS
        case "$GT" in
            "0/0") AF_NUM="0.0000"; STATUS="PASS" ;;
            "0/1") AF_NUM="0.5000"; STATUS="FAIL" ;;
            "1/1") AF_NUM="1.0000"; STATUS="FAIL" ;;
            "./.") AF_NUM="0.0000"; STATUS="NOCOV"; DP_NUM=0 ;;
            *)     AF_NUM="0.0000"; STATUS="NOCOV" ;;
        esac

        # Coverage rule
        if [[ "$DP_NUM" -eq 0 ]]; then
            STATUS="NOCOV"
        fi
    fi
fi

# ============================================================
# SAVE RESULTS
# ============================================================

echo "biosample,vcf,genotype,position,allele_frequency,depth,status" > "$SUMMARY_CSV"
echo "${BIOSAMPLE},$(basename "$VCF_LOFREQ"),${GT},${NTM_POS},${AF_NUM},${DP_NUM},${STATUS}" >> "$SUMMARY_CSV"

echo "[OK] Summary CSV generated: ${SUMMARY_CSV}"
echo "[DONE] NTM filter completed for biosample: ${BIOSAMPLE}"
echo "[OUT] Status: $STATUS"
echo "[OUT] AF: $AF_NUM"
echo "[OUT] DP: $DP_NUM"
echo "[OUT] Saved in: ${OUTPUT_DIR}"
