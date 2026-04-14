#!/usr/bin/env python3

################################################################################
# Script: add_cpm_values.py
# Description: Add CPM (Counts Per Million) normalization to Kallisto abundance TSV
# Usage: python3 add_cpm_values.py INPUT_TSV OUTPUT_TSV
#
# Inputs:
#   INPUT_TSV:   Kallisto abundance.tsv file (with est_counts column)
#   OUTPUT_TSV:  Output file with added cpm column
#
# Outputs:
#   OUTPUT_TSV   TSV file with all original columns + new 'cpm' column
#
# Calculation:
#   CPM = (est_counts / sum(all est_counts)) * 1,000,000
#
# Notes:
#   - Preserves all original columns
#   - Adds 'cpm' as new rightmost column
#   - Handles headers automatically
#   - Skips malformed rows with warning
#
# Dependencies: Python 3 (only standard library)
#
################################################################################

import argparse
import sys

# ============================================================================
# FUNCTION: Read TSV file and add CPM column
# ============================================================================
def add_cpm_column(input_file, output_file):
    """
    Read Kallisto abundance TSV and add CPM column.

    Args:
        input_file: Path to input TSV file (Kallisto abundance.tsv)
        output_file: Path to output TSV file
    """
    rows = []
    header = None
    est_counts_index = None
    total_counts = 0.0
    valid_row_count = 0

    try:
        # First pass: read file, find est_counts column, and sum counts
        with open(input_file, 'r') as f:
            for line_num, line in enumerate(f, start=1):
                line = line.rstrip('\n')

                if line_num == 1:
                    # Header line
                    header = line.split('\t')
                    try:
                        est_counts_index = header.index('est_counts')
                    except ValueError:
                        print(f"Error: 'est_counts' column not found in header.", file=sys.stderr)
                        sys.exit(1)
                    rows.append(header)
                else:
                    # Data line
                    parts = line.split('\t')

                    if len(parts) != len(header):
                        print(f"Warning: Row {line_num} has {len(parts)} columns, expected {len(header)}. Skipping.",
                              file=sys.stderr)
                        continue

                    try:
                        est_count = float(parts[est_counts_index])
                        total_counts += est_count
                        valid_row_count += 1
                        rows.append(parts)
                    except (ValueError, IndexError) as e:
                        print(f"Warning: Row {line_num} has invalid est_counts value. Skipping.",
                              file=sys.stderr)
                        continue

    except IOError as e:
        print(f"Error reading input file: {e}", file=sys.stderr)
        sys.exit(1)

    if total_counts == 0:
        print("Warning: No valid counts found or sum is zero.", file=sys.stderr)
        total_counts = 1.0  # Prevent division by zero

    # Second pass: write output with CPM values
    try:
        with open(output_file, 'w') as out_f:
            for i, row in enumerate(rows):
                if i == 0:
                    # Header: add 'cpm' column
                    output_line = '\t'.join(row) + '\tcpm\n'
                else:
                    # Data: calculate and add CPM
                    try:
                        est_count = float(row[est_counts_index])
                        cpm = (est_count / total_counts) * 1_000_000
                        output_line = '\t'.join(row) + f'\t{cpm:.6f}\n'
                    except (ValueError, IndexError):
                        continue

                out_f.write(output_line)

    except IOError as e:
        print(f"Error writing output file: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Added CPM values: {valid_row_count} rows processed", file=sys.stderr)
    print(f"Total est_counts sum: {total_counts:,.0f}", file=sys.stderr)
    print(f"Output saved to: {output_file}")

# ============================================================================
# FUNCTION: Main entry point
# ============================================================================
def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Add CPM (Counts Per Million) normalization to Kallisto abundance data.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
CPM Calculation:
  CPM = (est_counts / sum(all est_counts)) * 1,000,000

Examples:
  python3 add_cpm_values.py abundance.tsv abundance_with_cpm.tsv
  python3 add_cpm_values.py input.tsv output.tsv
        """
    )

    parser.add_argument("input_tsv",
                        help="Input Kallisto abundance TSV file")
    parser.add_argument("output_tsv",
                        help="Output TSV file with added CPM column")

    args = parser.parse_args()

    # Validate input file
    try:
        with open(args.input_tsv, 'r') as f:
            pass
    except FileNotFoundError:
        print(f"Error: Input file '{args.input_tsv}' not found.", file=sys.stderr)
        sys.exit(1)

    # Process file
    add_cpm_column(args.input_tsv, args.output_tsv)

if __name__ == "__main__":
    main()
