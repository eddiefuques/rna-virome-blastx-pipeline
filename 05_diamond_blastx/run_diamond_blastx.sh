#!/bin/bash

################################################################################
# Script: run_diamond_blastx.sh
# Description: Protein-level sequence alignment using DIAMOND blastx
# Usage: ./run_diamond_blastx.sh -c CONTIGS_DIR -o OUTPUT_DIR -d DB_PATH [-t THREADS]
#
# Inputs:
#   CONTIGS_DIR:  Directory containing *_rnaspades_min500bp_transcripts.fasta files
#   OUTPUT_DIR:   Directory for DIAMOND blastx output
#   DB_PATH:      Path to DIAMOND database (.dmnd file)
#   THREADS:      Number of threads (default: 20)
#
# Outputs:
#   {SAMPLE}_diamond_definitivo_tabular_out.txt - DIAMOND blastx results (tabular format)
#
# DIAMOND parameters (from original pipeline):
#   - Top hits: 3
#   - Sensitivity: sensitive
#   - Minimum bitscore: 300
#   - Minimum identity: 60%
#   - Output format: 6 (tabular)
#
# Dependencies: diamond
################################################################################

set -euo pipefail

# ============================================================================
# CONFIGURATION SECTION
# ============================================================================
THREADS=20
CONTIGS_DIR=""
OUTPUT_DIR=""
DB_PATH=""

# DIAMOND parameters
TOP_HITS=3
SENSITIVITY="sensitive"
MIN_BITSCORE=300
MIN_IDENTITY=60

# ============================================================================
# FUNCTION: Print usage
# ============================================================================
usage() {
    echo "Usage: $0 -c CONTIGS_DIR -o OUTPUT_DIR -d DB_PATH [-t THREADS]"
    echo ""
    echo "  -c CONTIGS_DIR  Directory containing *_rnaspades_min500bp_transcripts.fasta"
    echo "  -o OUTPUT_DIR   Directory for DIAMOND blastx output"
    echo "  -d DB_PATH      Path to DIAMOND database file (.dmnd)"
    echo "  -t THREADS      Number of threads (default: 20)"
    echo ""
    exit 1
}

# ============================================================================
# Parse command-line arguments
# ============================================================================
while getopts "c:o:d:t:h" opt; do
    case $opt in
        c) CONTIGS_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        d) DB_PATH="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Validate required arguments
if [[ -z "$CONTIGS_DIR" || -z "$OUTPUT_DIR" || -z "$DB_PATH" ]]; then
    echo "Error: -c, -o, and -d options are required."
    usage
fi

# Validate input directory
if [[ ! -d "$CONTIGS_DIR" ]]; then
    echo "Error: CONTIGS_DIR '$CONTIGS_DIR' does not exist."
    exit 1
fi

# Validate database
if [[ ! -f "$DB_PATH" ]]; then
    echo "Error: Database file '$DB_PATH' not found."
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# ============================================================================
# Run DIAMOND blastx on all contig files
# ============================================================================
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting DIAMOND blastx..."

SAMPLE_COUNT=0
for CONTIG_FILE in "$CONTIGS_DIR"/*_rnaspades_min500bp_transcripts.fasta; do
    # Check if files exist (in case glob expands to no matches)
    if [[ ! -f "$CONTIG_FILE" ]]; then
        echo "Warning: No contig files found matching pattern"
        continue
    fi

    SAMPLE_COUNT=$((SAMPLE_COUNT + 1))
    SAMPLE_NAME=$(basename "$CONTIG_FILE" _rnaspades_min500bp_transcripts.fasta)
    OUTPUT_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}_diamond_definitivo_tabular_out.txt"
    TEMP_DIR="${OUTPUT_DIR}/.diamond_temp_${SAMPLE_NAME}"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing sample $SAMPLE_COUNT: $SAMPLE_NAME"

    # Create temporary directory for DIAMOND
    mkdir -p "$TEMP_DIR"

    # Run DIAMOND blastx
    diamond blastx \
        -d "$DB_PATH" \
        -q "$CONTIG_FILE" \
        --out "$OUTPUT_FILE" \
        -t "$TEMP_DIR" \
        --top "$TOP_HITS" \
        --"$SENSITIVITY" \
        --min-score "$MIN_BITSCORE" \
        --id "$MIN_IDENTITY" \
        --outfmt 6 \
        -p "$THREADS"

    # Clean up temporary directory
    rm -rf "$TEMP_DIR"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')]   Completed: $OUTPUT_FILE"
done

if [[ $SAMPLE_COUNT -eq 0 ]]; then
    echo "Warning: No input files found matching pattern *_rnaspades_min500bp_transcripts.fasta"
    exit 1
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] DIAMOND blastx complete. Processed $SAMPLE_COUNT samples."
