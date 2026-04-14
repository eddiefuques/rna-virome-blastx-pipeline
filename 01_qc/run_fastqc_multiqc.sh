#!/bin/bash

################################################################################
# Script: run_fastqc_multiqc.sh
# Description: Run FastQC on FASTQ files, then aggregate results with MultiQC
# Usage: ./run_fastqc_multiqc.sh -i INPUT_DIR -o OUTPUT_DIR [-t THREADS]
#
# Inputs:
#   INPUT_DIR:  Directory containing .fastq.gz files
#   OUTPUT_DIR: Directory for FastQC and MultiQC results
#   THREADS:    Number of threads for FastQC (default: 10)
#
# Outputs:
#   OUTPUT_DIR/fastqc_results/      - Individual FastQC HTML/ZIP files
#   OUTPUT_DIR/multiqc_report.html  - Aggregated MultiQC report
#
# Dependencies: fastqc, multiqc
################################################################################

set -euo pipefail

# ============================================================================
# CONFIGURATION SECTION
# ============================================================================
THREADS=10
INPUT_DIR=""
OUTPUT_DIR=""

# ============================================================================
# FUNCTION: Print usage
# ============================================================================
usage() {
    echo "Usage: $0 -i INPUT_DIR -o OUTPUT_DIR [-t THREADS]"
    echo ""
    echo "  -i INPUT_DIR   Directory containing .fastq.gz files"
    echo "  -o OUTPUT_DIR  Directory for FastQC and MultiQC output"
    echo "  -t THREADS     Number of threads for FastQC (default: 10)"
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
# Run FastQC
# ============================================================================
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting FastQC on input files..."
fastqc -t "$THREADS" "$INPUT_DIR"/*.fastq.gz -o "$OUTPUT_DIR"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] FastQC completed."

# ============================================================================
# Run MultiQC
# ============================================================================
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting MultiQC aggregation..."
multiqc "$OUTPUT_DIR" -o "$OUTPUT_DIR" -n multiqc_report
echo "[$(date '+%Y-%m-%d %H:%M:%S')] MultiQC completed."

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Pipeline complete. Results in: $OUTPUT_DIR"
