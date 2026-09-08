#!/usr/bin/env bash
# ==============================================================================
# Script: trim_16s.sh
# Purpose: Trim V3-V4 primers from paired-end 16S rRNA sequencing reads using Cutadapt
# Primers:
#   - Forward (341F): CCTACGGGNGGCWGCAG (or ^CCTAYGGGDBGCWGCAG)
#   - Reverse (806R): GGACTACNVGGGTWTCTAAT (or ^GACTACNVGGGTMTCTAATCC)
# ==============================================================================

set -euo pipefail

# Ensure cutadapt is available
if ! command -v cutadapt &> /dev/null; then
    echo "Error: cutadapt could not be found. Please activate your conda environment (e.g., conda activate cutadapt_env)." >&2
    exit 1
fi

echo "Starting 16S rRNA V3-V4 primer trimming..."

for R1 in *_R1*.fastq*; do
    # Skip files that have already been trimmed or don't match
    if [[ ! -f "$R1" ]] || [[ "$R1" == trim_* ]]; then
        continue
    fi

    # Determine base sample name and file extensions
    if [[ "$R1" == *_R1_001.fastq.gz ]]; then
        BASE="${R1%_R1_001.fastq.gz}"
        R2="${BASE}_R2_001.fastq.gz"
        OUT_R1="trim_${BASE}_R1_001.fastq.gz"
        OUT_R2="trim_${BASE}_R2_001.fastq.gz"
    elif [[ "$R1" == *_R1.fastq.gz ]]; then
        BASE="${R1%_R1.fastq.gz}"
        R2="${BASE}_R2.fastq.gz"
        OUT_R1="trim_${BASE}_R1.fastq.gz"
        OUT_R2="trim_${BASE}_R2.fastq.gz"
    elif [[ "$R1" == *_R1.fastq ]]; then
        BASE="${R1%_R1.fastq}"
        R2="${BASE}_R2.fastq"
        OUT_R1="trim_${BASE}_R1.fastq.gz"
        OUT_R2="trim_${BASE}_R2.fastq.gz"
    else
        echo "Skipping unrecognized naming format: $R1"
        continue
    fi

    # Verify matching reverse read exists
    if [[ ! -f "$R2" ]]; then
        echo "Warning: Matching reverse read not found for $R1 (expected $R2). Skipping." >&2
        continue
    fi

    echo "----------------------------------------"
    echo "Processing sample: $BASE"
    echo "  Forward: $R1"
    echo "  Reverse: $R2"

    cutadapt \
        -g CCTACGGGNGGCWGCAG \
        -G GGACTACNVGGGTWTCTAAT \
        --discard-untrimmed \
        -o "$OUT_R1" \
        -p "$OUT_R2" \
        "$R1" "$R2"

    echo "  -> Trimmed output saved to: $OUT_R1 and $OUT_R2"
done

echo "----------------------------------------"
echo "Trimming complete! All trimmed files prefixed with 'trim_'."
