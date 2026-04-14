#!/bin/bash

################################################################################
# Script: run_kallisto_pipeline.sh
# Description: Quantify transcript abundance using Kallisto
# Usage: ./run_kallisto_pipeline.sh -c CONTIGS_DIR -r READS_DIR -o OUTPUT_DIR [-t THREADS] [-b BOOTSTRAPS]
#
# Inputs:
#   CONTIGS_DIR:  Directory containing *_viral_contigs_only.fasta files
#   READS_DIR:    Directory containing *_R1_trimmed.fq.gz and *_R2_trimmed.fq.gz files
#   OUTPUT_DIR:   Directory for Kallisto output
#   THREADS:      Number of threads (default: 33)
#   BOOTSTRAPS:   Number of bootstrap samples (default: 100)
#
# Outputs:
#   {SAMPLE}_kallisto.idx                    - Kallisto index
#   {SAMPLE}_kallisto_quant_output/
#     ├── abundance.tsv                       - Abundance estimates (TPM, counts)
#     └── abundance.h5                        - HDF5 format abundance file
#   {SAMPLE}_abundance.tsv                   - Copy of abundance.tsv for analysis
#
# Pipeline steps:
#   1. Build Kallisto index from viral contigs
#   2. Quantify abundances using trimmed paired-end reads
#   3. Copy abundance.tsv to output directory with sample name
#
# Dependencies: kallisto
#
################################################################################

set -euo pipefail

# ============================================================================
# CONFIGURATION SECTION
# ============================================================================
THREADS=33
BOOTSTRAPS=100
CONTIGS_DIR=""
READS_DIR=""
OUTPUT_DIR=""

# ============================================================================
# FUNCTION: Print usage
# ============================================================================
usage() {
    echo "Usage: $0 -c CONTIGS_DIR -r READS_DIR -o OUTPUT_DIR [-t THREADS] [-b BOOTSTRAPS]"
    echo ""
    echo "  -c CONTIGS_DIR  Directory containing *_viral_contigs_only.fasta files"
    echo "  -r READS_DIR    Directory containing trimmed paired FASTQ files"
    echo "  -o OUTPUT_DIR   Directory for Kallisto output"
    echo "  -t THREADS      Number of threads (default: 33)"
    echo "  -b BOOTSTRAPS   Number of bootstrap samples (default: 100)"
    echo ""
    exit 1
}

# ============================================================================
# Parse command-line arguments
# ============================================================================
while getopts "c:r:o:t:b:h" opt; do
    case $opt in
        c) CONTIGS_DIR="$OPTARG" ;;
        r) READS_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        b) BOOTSTRAPS="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Validate required arguments
if [[ -z "$CONTIGS_DIR" || -z "$READS_DIR" || -z "$OUTPUT_DIR" ]]; then
    echo "Error: -c, -r, and -o options are required."
    usage
fi

# Validate input directories
if [[ ! -d "$CONTIGS_DIR" ]]; then
    echo "Error: CONTIGS_DIR '$CONTIGS_DIR' does not exist."
    exit 1
fi

if [[ ! -d "$READS_DIR" ]]; then
    echo "Error: READS_DIR '$READS_DIR' does not exist."
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# ============================================================================
# Extract unique sample names and run Kallisto pipeline
# ============================================================================
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting Kallisto abundance quantification..."

# Extract unique sample prefixes (remove trailing _viral_contigs_only.fasta)
SAMPLES=$(ls -1 "$CONTIGS_DIR"/*_viral_contigs_only.fasta 2>/dev/null | \
          xargs -I {} basename {} _viral_contigs_only.fasta | sort -u)

if [[ -z "$SAMPLES" ]]; then
    echo "Error: No input files matching pattern *_viral_contigs_only.fasta found in $CONTIGS_DIR"
    exit 1
fi

SAMPLE_COUNT=0
for SAMPLE in $SAMPLES; do
    SAMPLE_COUNT=$((SAMPLE_COUNT + 1))
    CONTIGS_FILE="${CONTIGS_DIR}/${SAMPLE}_viral_contigs_only.fasta"
    R1_READS="${READS_DIR}/${SAMPLE}_R1_trimmed.fq.gz"
    R2_READS="${READS_DIR}/${SAMPLE}_R2_trimmed.fq.gz"
    KALLISTO_INDEX="${OUTPUT_DIR}/${SAMPLE}_kallisto.idx"
    QUANT_OUTPUT="${OUTPUT_DIR}/${SAMPLE}_kallisto_quant_output"
    ABUNDANCE_COPY="${OUTPUT_DIR}/${SAMPLE}_abundance.tsv"

    # Check input files exist
    if [[ ! -f "$CONTIGS_FILE" ]]; then
        echo "Warning: Contigs file not found: $CONTIGS_FILE"
        continue
    fi
    if [[ ! -f "$R1_READS" ]]; then
        echo "Warning: R1 reads file not found: $R1_READS"
        continue
    fi
    if [[ ! -f "$R2_READS" ]]; then
        echo "Warning: R2 reads file not found: $R2_READS"
        continue
    fi

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing sample $SAMPLE_COUNT: $SAMPLE"

    # ========================================================================
    # Step 1: Build Kallisto index
    # ========================================================================
    echo "[$(date '+%Y-%m-%d %H:%M:%S')]   Building Kallisto index..."
    kallisto index \
        --make-unique \
        -i "$KALLISTO_INDEX" \
        "$CONTIGS_FILE"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')]   Index built: $KALLISTO_INDEX"

    # ========================================================================
    # Step 2: Run Kallisto quantification
    # ========================================================================
    echo "[$(date '+%Y-%m-%d %H:%M:%S')]   Running Kallisto quantification..."
    kallisto quant \
        --plaintext \
        -i "$KALLISTO_INDEX" \
        "$R1_READS" "$R2_READS" \
        -o "$QUANT_OUTPUT" \
        -t "$THREADS" \
        -b "$BOOTSTRAPS" \
        --bias

    echo "[$(date '+%Y-%m-%d %H:%M:%S')]   Quantification complete: $QUANT_OUTPUT"

    # ========================================================================
    # Step 3: Copy abundance.tsv to output directory with sample name
    # ========================================================================
    if [[ -f "$QUANT_OUTPUT/abundance.tsv" ]]; then
        cp "$QUANT_OUTPUT/abundance.tsv" "$ABUNDANCE_COPY"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')]   Copied abundance file: $ABUNDANCE_COPY"
    else
        echo "Warning: abundance.tsv not found in $QUANT_OUTPUT"
    fi
done

if [[ $SAMPLE_COUNT -eq 0 ]]; then
    echo "Error: No samples processed successfully"
    exit 1
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Kallisto quantification complete. Processed $SAMPLE_COUNT samples."
