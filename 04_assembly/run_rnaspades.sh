#!/bin/bash

################################################################################
# Script: run_rnaspades.sh
# Description: De novo assembly of normalized RNA-seq reads using rnaSPAdes
#              Automatically filters transcripts by minimum length (500bp)
# Usage: ./run_rnaspades.sh -i INPUT_DIR -o OUTPUT_DIR [-t THREADS] [-k KMERS]
#
# Inputs:
#   INPUT_DIR:  Directory containing *_R1_norm.fq.gz and *_R2_norm.fq.gz files
#   OUTPUT_DIR: Directory for assembly output
#   THREADS:    Number of threads (default: 20)
#   KMERS:      Comma-separated k-mer values (default: 77,99,127)
#
# Outputs:
#   {SAMPLE}/{SAMPLE}_rnaspades_alltranscripts.fasta    - All transcripts
#   {SAMPLE}/{SAMPLE}_rnaspades_min500bp_transcripts.fasta - Filtered (>=500bp)
#
# Dependencies: rnaspades.py, perl (for filter_contigs_by_length.pl)
# Notes:
#   - Requires filter_contigs_by_length.pl in same directory or accessible via PATH
#   - Assembly can be computationally intensive; allocate sufficient memory
#
################################################################################

set -euo pipefail

# ============================================================================
# CONFIGURATION SECTION
# ============================================================================
THREADS=20
KMERS="77,99,127"
MIN_LENGTH=500
INPUT_DIR=""
OUTPUT_DIR=""

# Script directory (for locating filter script)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ============================================================================
# FUNCTION: Print usage
# ============================================================================
usage() {
    echo "Usage: $0 -i INPUT_DIR -o OUTPUT_DIR [-t THREADS] [-k KMERS]"
    echo ""
    echo "  -i INPUT_DIR   Directory containing normalized paired FASTQ files"
    echo "  -o OUTPUT_DIR  Directory for assembly output"
    echo "  -t THREADS     Number of threads (default: 20)"
    echo "  -k KMERS       Comma-separated k-mer values (default: 77,99,127)"
    echo ""
    exit 1
}

# ============================================================================
# Parse command-line arguments
# ============================================================================
while getopts "i:o:t:k:h" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        k) KMERS="$OPTARG" ;;
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

# Check for filter script
if [[ ! -f "$SCRIPT_DIR/filter_contigs_by_length.pl" ]]; then
    echo "Error: filter_contigs_by_length.pl not found in $SCRIPT_DIR"
    exit 1
fi

# Make filter script executable
chmod +x "$SCRIPT_DIR/filter_contigs_by_length.pl"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# ============================================================================
# Extract unique sample names and run rnaSPAdes
# ============================================================================
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting rnaSPAdes assembly..."

# Extract unique sample prefixes (remove trailing _R1_norm.fq.gz)
SAMPLES=$(ls -1 "$INPUT_DIR"/*_R1_norm.fq.gz 2>/dev/null | \
          xargs -I {} basename {} _R1_norm.fq.gz | sort -u)

if [[ -z "$SAMPLES" ]]; then
    echo "Error: No input files matching pattern *_R1_norm.fq.gz found in $INPUT_DIR"
    exit 1
fi

SAMPLE_COUNT=0
for SAMPLE in $SAMPLES; do
    SAMPLE_COUNT=$((SAMPLE_COUNT + 1))
    R1_INPUT="${INPUT_DIR}/${SAMPLE}_R1_norm.fq.gz"
    R2_INPUT="${INPUT_DIR}/${SAMPLE}_R2_norm.fq.gz"
    SAMPLE_OUTPUT="${OUTPUT_DIR}/${SAMPLE}_spades_output"

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

    # Run rnaSPAdes
    rnaspades.py \
        -1 "$R1_INPUT" \
        -2 "$R2_INPUT" \
        -k "$KMERS" \
        -t "$THREADS" \
        -o "$SAMPLE_OUTPUT"

    # Rename transcripts.fasta to include sample name
    if [[ -f "$SAMPLE_OUTPUT/transcripts.fasta" ]]; then
        mv "$SAMPLE_OUTPUT/transcripts.fasta" \
           "$SAMPLE_OUTPUT/${SAMPLE}_rnaspades_alltranscripts.fasta"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')]   Assembly complete: ${SAMPLE}_rnaspades_alltranscripts.fasta"

        # Filter by minimum length (500bp)
        echo "[$(date '+%Y-%m-%d %H:%M:%S')]   Filtering contigs (min length: ${MIN_LENGTH}bp)..."
        "$SCRIPT_DIR/filter_contigs_by_length.pl" "$MIN_LENGTH" \
            "$SAMPLE_OUTPUT/${SAMPLE}_rnaspades_alltranscripts.fasta" > \
            "$SAMPLE_OUTPUT/${SAMPLE}_rnaspades_min${MIN_LENGTH}bp_transcripts.fasta"

        echo "[$(date '+%Y-%m-%d %H:%M:%S')]   Filtered output: ${SAMPLE}_rnaspades_min${MIN_LENGTH}bp_transcripts.fasta"
    else
        echo "Warning: transcripts.fasta not found in $SAMPLE_OUTPUT"
    fi
done

echo "[$(date '+%Y-%m-%d %H:%M:%S')] rnaSPAdes assembly and filtering complete. Processed $SAMPLE_COUNT samples."
