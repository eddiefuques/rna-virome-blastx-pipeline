#!/usr/bin/env python3

################################################################################
# Script: add_ppm_values.py
# Description: Add PPM (Parts Per Million relative to sample total reads) to Kallisto abundance
# Usage: python3 add_ppm_values.py KALLISTO_TSV READS_MAPPING_TSV SAMPLE_NAME OUTPUT_TSV
#
# Inputs:
#   KALLISTO_TSV:     Kallisto abundance.tsv file (or abundance with CPM)
#   READS_MAPPING_TSV: TSV file with columns: sample_id, total_reads
#   SAMPLE_NAME:      Sample identifier to look up in reads mapping file
#   OUTPUT_TSV:       Output file with added relative_abundance (PPM) column
#
# Outputs:
#   OUTPUT_TSV   TSV file with all original columns + new 'relative_abundance' column
#
# Calculation:
#   relative_abundance (PPM) = (est_counts / total_reads_for_sample) * 1,000,000
#
# Notes:
#   - Preserves all original columns
#   - Adds 'relative_abundance' as new rightmost column
#   - Reads mapping file must have 'sample_id' and 'total_reads' columns
#   - Handles headers automatically
#   - Skips malformed rows with warning
#
# Dependencies: Python 3 (only standard library)
#
################################################################################

import argparse
import sys

# ============================================================================
# FUNCTION: Read sample-to-total-reads mapping
# ============================================================================
def read_sample_mapping(mapping_file, sample_name):
    """
    Read sample mapping file and extract total reads for given sample.

    Args:
        mapping_file: Path to TSV with columns: sample_id, total_reads
        sample_name: Sample identifier to look up

    Returns:
        Total reads for the sample
    """
    try:
        with open(mapping_file, 'r') as f:
            header = None
            sample_id_index = None
            total_reads_index = None

            for line_num, line in enumerate(f, start=1):
                line = line.rstrip('\n')
                parts = line.split('\t')

                if line_num == 1:
                    # Header line
                    header = parts
                    try:
                        sample_id_index = header.index('sample_id')
                        total_reads_index = header.index('total_reads')
                    except ValueError:
                        print(f"Error: Required columns 'sample_id' and/or 'total_reads' not found in header.",
                              file=sys.stderr)
                        sys.exit(1)
                else:
                    # Data line
                    if len(parts) != len(header):
                        continue

                    try:
                        current_sample = parts[sample_id_index].strip()
                        if current_sample == sample_name:
                            total_reads = float(parts[total_reads_index])
                            return total_reads
                    except (ValueError, IndexError):
                        continue

        print(f"Error: Sample '{sample_name}' not found in mapping file.",
              file=sys.stderr)
        sys.exit(1)

    except IOError as e:
        print(f"Error reading mapping file: {e}", file=sys.stderr)
        sys.exit(1)

# ============================================================================
# FUNCTION: Add PPM column to abundance TSV
# ============================================================================
def add_ppm_column(kallisto_file, output_file, total_reads):
    """
    Add relative abundance (PPM) column to Kallisto abundance TSV.

    Args:
        kallisto_file: Path to input TSV file (Kallisto abundance)
        output_file: Path to output TSV file
        total_reads: Total number of reads for this sample
    """
    rows = []
    header = None
    est_counts_index = None
    valid_row_count = 0

    try:
        # Read file and find est_counts column
        with open(kallisto_file, 'r') as f:
            for line_num, line in enumerate(f, start=1):
                line = line.rstrip('\n')

                if line_num == 1:
                    # Header line
                    header = line.split('\t')
                    try:
                        est_counts_index = header.index('est_counts')
                    except ValueError:
                        print(f"Error: 'est_counts' column not found in header.",
                              file=sys.stderr)
                        sys.exit(1)
                    rows.append(header)
                else:
                    # Data line
                    parts = line.split('\t')

                    if len(parts) != len(header):
                        print(f"Warning: Row {line_num} column count mismatch. Skipping.",
                              file=sys.stderr)
                        continue

                    try:
                        est_count = float(parts[est_counts_index])
                        valid_row_count += 1
                        rows.append(parts)
                    except (ValueError, IndexError):
                        print(f"Warning: Row {line_num} has invalid est_counts value. Skipping.",
                              file=sys.stderr)
                        continue

    except IOError as e:
        print(f"Error reading Kallisto file: {e}", file=sys.stderr)
        sys.exit(1)

    if total_reads <= 0:
        print(f"Error: Total reads must be positive, got {total_reads}.",
              file=sys.stderr)
        sys.exit(1)

    # Write output with PPM values
    try:
        with open(output_file, 'w') as out_f:
            for i, row in enumerate(rows):
                if i == 0:
                    # Header: add 'relative_abundance' column
                    output_line = '\t'.join(row) + '\trelative_abundance\n'
                else:
                    # Data: calculate and add PPM
                    try:
                        est_count = float(row[est_counts_index])
                        ppm = (est_count / total_reads) * 1_000_000
                        output_line = '\t'.join(row) + f'\t{ppm:.6f}\n'
                    except (ValueError, IndexError):
                        continue

                out_f.write(output_line)

    except IOError as e:
        print(f"Error writing output file: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Added PPM (relative_abundance) values: {valid_row_count} rows processed",
          file=sys.stderr)
    print(f"Sample total reads: {total_reads:,.0f}", file=sys.stderr)
    print(f"Output saved to: {output_file}")

# ============================================================================
# FUNCTION: Main entry point
# ============================================================================
def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Add PPM (relative abundance) to Kallisto abundance data based on sample total reads.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
PPM Calculation:
  relative_abundance (PPM) = (est_counts / total_reads_for_sample) * 1,000,000

Reads mapping file format (tab-separated):
  sample_id    total_reads
  sample1      50000000
  sample2      75000000

Examples:
  python3 add_ppm_values.py abundance.tsv reads.tsv sample1 abundance_with_ppm.tsv
  python3 add_ppm_values.py kallisto.tsv samples_reads.tsv my_sample output.tsv
        """
    )

    parser.add_argument("kallisto_tsv",
                        help="Input Kallisto abundance TSV file")
    parser.add_argument("reads_mapping_tsv",
                        help="TSV file with sample_id and total_reads columns")
    parser.add_argument("sample_name",
                        help="Sample identifier to look up in reads mapping file")
    parser.add_argument("output_tsv",
                        help="Output TSV file with added relative_abundance column")

    args = parser.parse_args()

    # Validate input files
    try:
        with open(args.kallisto_tsv, 'r') as f:
            pass
    except FileNotFoundError:
        print(f"Error: Kallisto TSV file '{args.kallisto_tsv}' not found.",
              file=sys.stderr)
        sys.exit(1)

    try:
        with open(args.reads_mapping_tsv, 'r') as f:
            pass
    except FileNotFoundError:
        print(f"Error: Reads mapping file '{args.reads_mapping_tsv}' not found.",
              file=sys.stderr)
        sys.exit(1)

    # Get total reads for sample
    total_reads = read_sample_mapping(args.reads_mapping_tsv, args.sample_name)
    print(f"Found sample '{args.sample_name}' with {total_reads:,.0f} total reads.",
          file=sys.stderr)

    # Add PPM column
    add_ppm_column(args.kallisto_tsv, args.output_tsv, total_reads)

if __name__ == "__main__":
    main()
