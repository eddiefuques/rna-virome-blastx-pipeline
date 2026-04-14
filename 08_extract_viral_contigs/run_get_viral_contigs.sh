#!/bin/bash

################################################################################
# Script: run_get_viral_contigs.sh
# Description: Bash wrapper to extract viral contigs for all samples
# Usage: ./run_get_viral_contigs.sh -b BLAST_DIR -f FASTA_DIR -o OUTPUT_DIR [--min_id MIN_ID] [-p PYTHON_SCRIPT]
#
# Inputs:
#   BLAST_DIR:      Directory containing *_diamond_best_hits.txt files
#   FASTA_DIR:      Directory containing *_rnaspades_min500bp_transcripts.fasta files
#   OUTPUT_DIR:     Directory for viral contigs output
#   MIN_ID:         Minimum percent identity (default: 70)
#   PYTHON_SCRIPT:  Path to get_viral_contigs.py (default: same directory as this script)
#
# Outputs:
#   {SAMPLE}_viral_contigs_only.fasta
#
# Dependencies: Python 3, get_viral_contigs.py
#
################################################################################

set -euo pipefail

# ============================================================================
# CONFIGURATION SECTION
# ============================================================================
BLAST_DIR=""
FASTA_DIR=""
OUTPUT_DIR=""
MIN_ID=70
PYTHON_SCRIPT=""

# Script directory (for locating Python script)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ============================================================================
# FUNCTION: Print usage
# ============================================================================
usage() {
    echo "Usage: $0 -b BLAST_DIR -f FASTA_DIR -o OUTPUT_DIR [--min_id MIN_ID] [-p PYTHON_SCRIPT]"
    echo ""
    echo "  -b BLAST_DIR       Directory containing *_diamond_best_hits.txt files"
    echo "  -f FASTA_DIR       Directory containing *_rnaspades_min500bp_transcripts.fasta files"
    echo "  -o OUTPUT_DIR      Directory for viral contigs output"
    echo "  --min_id MIN_ID    Minimum percent identity threshold (default: 70)"
    echo "  -p PYTHON_SCRIPT   Path to get_viral_contigs.py (default: $SCRIPT_DIR/get_viral_contigs.py)"
    echo ""
    exit 1
}

# ============================================================================
# Parse command-line arguments
# ============================================================================
while [[ $# -gt 0 ]]; do
    case $1 in
        -b)
            BLAST_DIR="$2"
            shift 2
            ;;
        -f)
            FASTA_DIR="$2"
            shift 2
            ;;
        -o)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --min_id)
            MIN_ID="$2"
            shift 2
            ;;
        -p)
            PYTHON_SCRIPT="$2"
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
if [[ -z "$BLAST_DIR" || -z "$FASTA_DIR" || -z "$OUTPUT_DIR" ]]; then
    echo "Error: -b, -f, and -o options are required."
    usage
fi

# Set default Python script path if not provided
if [[ -z "$PYTHON_SCRIPT" ]]; then
    PYTHON_SCRIPT="$SCRIPT_DIR/get_viral_contigs.py"
fi

# Validate directories
if [[ ! -d "$BLAST_DIR" ]]; then
    echo "Error: BLAST_DIR '$BLAST_DIR' does not exist."
    exit 1
fi

if [[ ! -d "$FASTA_DIR" ]]; then
    echo "Error: FASTA_DIR '$FASTA_DIR' does not exist."
    exit 1
fi

# Validate Python script
if [[ ! -f "$PYTHON_SCRIPT" ]]; then
    echo "Error: Python script '$PYTHON_SCRIPT' not found."
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# ============================================================================
# Process all samples
# ============================================================================
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting viral contigs extraction (min_id >= $MIN_ID%)..."

SAMPLE_COUNT=0
for BLAST_FILE in "$BLAST_DIR"/*_diamond_best_hits.txt; do
    # Check if files exist (in case glob expands to no matches)
    if [[ ! -f "$BLAST_FILE" ]]; then
        echo "Warning: No BLAST files found matching pattern"
        continue
    fi

    SAMPLE_COUNT=$((SAMPLE_COUNT + 1))
    SAMPLE_NAME=$(basename "$BLAST_FILE" _diamond_best_hits.txt)
    FASTA_FILE="${FASTA_DIR}/${SAMPLE_NAME}_rnaspades_min500bp_transcripts.fasta"
    OUTPUT_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}_viral_contigs_only.fasta"

    # Check if corresponding FASTA file exists
    if [[ ! -f "$FASTA_FILE" ]]; then
        echo "Warning: FASTA file not found for sample $SAMPLE_NAME: $FASTA_FILE"
        continue
    fi

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing sample $SAMPLE_COUNT: $SAMPLE_NAME"

    # Run Python script to extract viral contigs
    python3 "$PYTHON_SCRIPT" \
        --blast "$BLAST_FILE" \
        --fasta "$FASTA_FILE" \
        --output "$OUTPUT_FILE" \
        --min_id "$MIN_ID"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')]   Completed: $OUTPUT_FILE"
done

if [[ $SAMPLE_COUNT -eq 0 ]]; then
    echo "Warning: No BLAST files found matching pattern *_diamond_best_hits.txt"
    exit 1
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Viral contigs extraction complete. Processed $SAMPLE_COUNT samples."
