#!/bin/bash

################################################################################
# Script: run_get_best_hits.sh
# Description: Bash wrapper to extract best BLAST hits for all samples
# Usage: ./run_get_best_hits.sh -i INPUT_DIR -o OUTPUT_DIR [-p PYTHON_SCRIPT]
#
# Inputs:
#   INPUT_DIR:      Directory containing *_diamond_tabular_out.txt files
#   OUTPUT_DIR:     Directory for best hits output
#   PYTHON_SCRIPT:  Path to get_best_hits.py (default: same directory as this script)
#
# Outputs:
#   {SAMPLE}_diamond_best_hits.txt
#
# Dependencies: Python 3, get_best_hits.py
#
################################################################################

set -euo pipefail

# ============================================================================
# CONFIGURATION SECTION
# ============================================================================
INPUT_DIR=""
OUTPUT_DIR=""
PYTHON_SCRIPT=""

# Script directory (for locating Python script)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ============================================================================
# FUNCTION: Print usage
# ============================================================================
usage() {
    echo "Usage: $0 -i INPUT_DIR -o OUTPUT_DIR [-p PYTHON_SCRIPT]"
    echo ""
    echo "  -i INPUT_DIR       Directory containing *_diamond_tabular_out.txt files"
    echo "  -o OUTPUT_DIR      Directory for best hits output"
    echo "  -p PYTHON_SCRIPT   Path to get_best_hits.py (default: $SCRIPT_DIR/get_best_hits.py)"
    echo ""
    exit 1
}

# ============================================================================
# Parse command-line arguments
# ============================================================================
while getopts "i:o:p:h" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        p) PYTHON_SCRIPT="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Validate required arguments
if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" ]]; then
    echo "Error: -i and -o options are required."
    usage
fi

# Set default Python script path if not provided
if [[ -z "$PYTHON_SCRIPT" ]]; then
    PYTHON_SCRIPT="$SCRIPT_DIR/get_best_hits.py"
fi

# Validate input directory
if [[ ! -d "$INPUT_DIR" ]]; then
    echo "Error: INPUT_DIR '$INPUT_DIR' does not exist."
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
# Process all DIAMOND tabular output files
# ============================================================================
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting best hits extraction..."

SAMPLE_COUNT=0
for INPUT_FILE in "$INPUT_DIR"/*_diamond_tabular_out.txt; do
    # Check if files exist (in case glob expands to no matches)
    if [[ ! -f "$INPUT_FILE" ]]; then
        echo "Warning: No input files found matching pattern"
        continue
    fi

    SAMPLE_COUNT=$((SAMPLE_COUNT + 1))
    SAMPLE_NAME=$(basename "$INPUT_FILE" _diamond_tabular_out.txt)
    OUTPUT_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}_diamond_best_hits.txt"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing sample $SAMPLE_COUNT: $SAMPLE_NAME"

    # Run Python script to extract best hits
    python3 "$PYTHON_SCRIPT" "$INPUT_FILE" "$OUTPUT_FILE"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')]   Completed: $OUTPUT_FILE"
done

if [[ $SAMPLE_COUNT -eq 0 ]]; then
    echo "Warning: No input files found matching pattern *_diamond_tabular_out.txt"
    exit 1
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Best hits extraction complete. Processed $SAMPLE_COUNT samples."
