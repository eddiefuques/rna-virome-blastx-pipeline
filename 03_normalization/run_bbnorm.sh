#!/bin/bash

################################################################################
# Script: run_bbnorm.sh
# Description: Read depth normalization using BBTools bbnorm.sh
# Usage: ./run_bbnorm.sh -i INPUT_DIR -o OUTPUT_DIR [-t THREADS] [--target TARGET]
#
# Inputs:
#   INPUT_DIR:  Directory containing *_R1_trimmed.fq.gz and *_R2_trimmed.fq.gz
#   OUTPUT_DIR: Directory for normalized FASTQ files
#   THREADS:    Number of threads (default: 20)
#   TARGET:     Target coverage depth (default: 40)
#
# Outputs:
#   {SAMPLE}_R1_norm.fq.gz
#   {SAMPLE}_R2_norm.fq.gz
#
# Dependencies: bbnorm.sh (BBTools)
################################################################################

set -euo pipefail

# ============================================================================
# CONFIGURATION SECTION
# ============================================================================
THREADS=20
TARGET=40
INPUT_DIR=""
OUTPUT_DIR=""

# ============================================================================
# FUNCTION: Print usage
# ============================================================================
usage() {
    echo "Usage: $0 -i INPUT_DIR -o OUTPUT_DIR [-t THREADS] [--target TARGET]"
    echo ""
    echo "  -i INPUT_DIR    Directory containing trimmed paired FASTQ files"
    echo "  -o OUTPUT_DIR   Directory for normalized output"
    echo "  -t THREADS      Number of threads (default: 20)"
    echo "  --target TARGET Target coverage depth (default: 40)"
    echo ""
    exit 1
}

# ============================================================================
# Parse command-line arguments
# ============================================================================
while [[ $# -gt 0 ]]; do
    case $1 in
        -i)
            INPUT_DIR="$2"
            shift 2
            ;;
        -o)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -t)
            THREADS="$2"
            shift 2
            ;;
        --target)
            TARGET="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Validate required arguments
if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" ]]; then
    echo "Error: -i and -o options are required."
    usage
fi

# Validate input directory
if [[ ! -d "$INPUT_DIR" ]]; then
    echo "Error: INPUT_DIR '$INPUT_DIR' does not exist."
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# ============================================================================
# Extract unique sample names and run bbnorm
# ============================================================================
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting bbnorm normalization..."

# Extract unique sample prefixes (remove trailing _R1_trimmed.fq.gz)
SAMPLES=$(ls -1 "$INPUT_DIR"/*_R1_trimmed.fq.gz 2>/dev/null | \
          xargs -I {} basename {} _R1_trimmed.fq.gz | sort -u)

if [[ -z "$SAMPLES" ]]; then
    echo "Error: No input files matching pattern *_R1_trimmed.fq.gz found in $INPUT_DIR"
    exit 1
fi

SAMPLE_COUNT=0
for SAMPLE in $SAMPLES; do
    SAMPLE_COUNT=$((SAMPLE_COUNT + 1))
    R1_INPUT="${INPUT_DIR}/${SAMPLE}_R1_trimmed.fq.gz"
    R2_INPUT="${INPUT_DIR}/${SAMPLE}_R2_trimmed.fq.gz"
    R1_OUTPUT="${OUTPUT_DIR}/${SAMPLE}_R1_norm.fq.gz"
    R2_OUTPUT="${OUTPUT_DIR}/${SAMPLE}_R2_norm.fq.gz"

    # Check input files exist
    if [[ ! -f "$R1_INPUT" ]]; then
        echo "Warning: R1 file not found: $R1_INPUT"
        continue
    fi
    if [[ ! -f "$R2_INPUT" ]]; then
        echo "Warning: R2 file not found: $R2_INPUT"
        continue
    fi

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing sample $SAMPLE_COUNT: $SAMPLE"

    bbnorm.sh \
        in1="$R1_INPUT" in2="$R2_INPUT" \
        out1="$R1_OUTPUT" out2="$R2_OUTPUT" \
        target="$TARGET" min=0 threads="$THREADS"
done

echo "[$(date '+%Y-%m-%d %H:%M:%S')] bbnorm normalization complete. Processed $SAMPLE_COUNT samples."
