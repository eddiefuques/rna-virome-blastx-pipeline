#!/bin/bash

################################################################################
# Script: run_fastp.sh
# Description: Quality and adapter trimming of paired-end reads using fastp
# Usage: ./run_fastp.sh -i INPUT_DIR -o OUTPUT_DIR [-t THREADS]
#
# Inputs:
#   INPUT_DIR:  Directory containing paired FASTQ files in format:
#               {SAMPLE}_L001_R1_001.fastq.gz and {SAMPLE}_L001_R2_001.fastq.gz
#   OUTPUT_DIR: Directory for trimmed FASTQ files and reports
#   THREADS:    Number of threads (default: 12)
#
# Outputs:
#   {SAMPLE}_R1_trimmed.fq.gz
#   {SAMPLE}_R2_trimmed.fq.gz
#   {SAMPLE}.fastp.html
#   {SAMPLE}.fastp.json
#
# Trimming parameters:
#   - Min average quality: 25
#   - Min base quality: 20
#   - Min read length: 50
#   - Min complexity: 30%
#   - Remove poly-A/G: yes
#   - Max Ns: 2
#   - Auto-detect adapters: yes
#
# Dependencies: fastp
################################################################################

set -euo pipefail

# ============================================================================
# CONFIGURATION SECTION
# ============================================================================
THREADS=12
INPUT_DIR=""
OUTPUT_DIR=""

# Trimming parameters (change as needed)
AVERAGE_QUAL=25
MIN_BASE_QUAL=20
MIN_LENGTH=50
MIN_COMPLEXITY=30
MAX_N=2

# ============================================================================
# FUNCTION: Print usage
# ============================================================================
usage() {
    echo "Usage: $0 -i INPUT_DIR -o OUTPUT_DIR [-t THREADS]"
    echo ""
    echo "  -i INPUT_DIR   Directory containing paired FASTQ files"
    echo "  -o OUTPUT_DIR  Directory for trimmed output"
    echo "  -t THREADS     Number of threads (default: 12)"
    echo ""
    exit 1
}

# ============================================================================
# Parse command-line arguments
# ============================================================================
while getopts "i:o:t:h" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
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
# Extract unique sample names and run fastp (fastq file name structure: {sample_ID}_L001_R1_001.fastq.gz)
# ============================================================================
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting fastp trimming..."

# Extract unique sample prefixes (remove trailing _L001_R1_001.fastq.gz)
SAMPLES=$(ls -1 "$INPUT_DIR"/*_L001_R1_001.fastq.gz 2>/dev/null | \
          xargs -I {} basename {} _L001_R1_001.fastq.gz | sort -u)

if [[ -z "$SAMPLES" ]]; then
    echo "Error: No input files matching pattern *_L001_R1_001.fastq.gz found in $INPUT_DIR"
    exit 1
fi

SAMPLE_COUNT=0
for SAMPLE in $SAMPLES; do
    SAMPLE_COUNT=$((SAMPLE_COUNT + 1))
    R1_INPUT="${INPUT_DIR}/${SAMPLE}_L001_R1_001.fastq.gz"
    R2_INPUT="${INPUT_DIR}/${SAMPLE}_L001_R2_001.fastq.gz"
    R1_OUTPUT="${OUTPUT_DIR}/${SAMPLE}_R1_trimmed.fq.gz"
    R2_OUTPUT="${OUTPUT_DIR}/${SAMPLE}_R2_trimmed.fq.gz"

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

    fastp \
        -i "$R1_INPUT" -I "$R2_INPUT" \
        -o "$R1_OUTPUT" -O "$R2_OUTPUT" \
        --detect_adapter_for_pe \
        --average_qual "$AVERAGE_QUAL" \
        -q "$MIN_BASE_QUAL" \
        -l "$MIN_LENGTH" \
        -y -Y "$MIN_COMPLEXITY" \
        -g -x \
        -n "$MAX_N" \
        -c \
        --overrepresentation_analysis \
        --html "${OUTPUT_DIR}/${SAMPLE}.fastp.html" \
        --json "${OUTPUT_DIR}/${SAMPLE}.fastp.json" \
        --thread "$THREADS" \
        --report_title "$SAMPLE"
done

echo "[$(date '+%Y-%m-%d %H:%M:%S')] fastp trimming complete. Processed $SAMPLE_COUNT samples."
