#!/usr/bin/env python3

################################################################################
# Script: add_tpm_to_blast.py
# Description: Join BLAST/taxonomy results with Kallisto TPM abundance estimates
# Usage: python3 add_tpm_to_blast.py ABUNDANCE_FILE BLAST_TAXONOMY_FILE OUTPUT_FILE
#
# Inputs:
#   ABUNDANCE_FILE:      Kallisto abundance.tsv file (contains target_id, tpm columns)
#   BLAST_TAXONOMY_FILE: BLAST results with taxonomy (from add_taxonomy_phyloR.R)
#   OUTPUT_FILE:         Output file with merged results
#
# Outputs:
#   OUTPUT_FILE  Combined table with columns:
#     query_acc_ver, subject_acc_ver, pident, length, mismatch, gapopen,
#     q_start, q_end, s_start, s_end, evalue, bit_score, positives,
#     taxid, kingdom, family, genus, species, tpm
#
# Logic:
#   - Reads Kallisto abundance (target_id -> tpm)
#   - Reads BLAST+taxonomy table (query_acc_ver -> various columns)
#   - Joins on query_acc_ver (BLAST) = target_id (Kallisto)
#   - Adds tpm column to each BLAST hit
#
# Notes:
#   - BLAST file should come from add_taxonomy_phyloR.R (has taxonomy columns)
#   - target_id in abundance must match query_acc_ver in BLAST
#   - Unmatched BLAST hits are kept with tpm = NA
#   - Preserves all original BLAST columns
#
# Dependencies: Python 3 (only standard library)
#
################################################################################

import argparse
import sys

# ============================================================================
# FUNCTION: Read Kallisto abundance file
# ============================================================================
def read_abundance_file(abundance_file):
    """
    Read Kallisto abundance file and create target_id -> tpm mapping.

    Args:
        abundance_file: Path to Kallisto abundance.tsv

    Returns:
        Dictionary mapping target_id to tpm value
    """
    abundance_map = {}

    try:
        with open(abundance_file, 'r') as f:
            header = None
            target_id_index = None
            tpm_index = None

            for line_num, line in enumerate(f, start=1):
                line = line.rstrip('\n')
                parts = line.split('\t')

                if line_num == 1:
                    # Header line
                    header = parts
                    try:
                        target_id_index = header.index('target_id')
                        tpm_index = header.index('tpm')
                    except ValueError:
                        print(f"Error: Required columns 'target_id' and/or 'tpm' not found in header.",
                              file=sys.stderr)
                        sys.exit(1)
                else:
                    # Data line
                    if len(parts) != len(header):
                        print(f"Warning: Row {line_num} column count mismatch. Skipping.",
                              file=sys.stderr)
                        continue

                    try:
                        target_id = parts[target_id_index].strip()
                        tpm = float(parts[tpm_index])
                        abundance_map[target_id] = tpm
                    except (ValueError, IndexError):
                        print(f"Warning: Row {line_num} has invalid data. Skipping.",
                              file=sys.stderr)
                        continue

    except IOError as e:
        print(f"Error reading abundance file: {e}", file=sys.stderr)
        sys.exit(1)

    return abundance_map

# ============================================================================
# FUNCTION: Join BLAST and abundance data
# ============================================================================
def join_blast_and_abundance(blast_file, abundance_map, output_file):
    """
    Read BLAST file, add TPM values, and write output.

    Args:
        blast_file: Path to BLAST+taxonomy file
        abundance_map: Dictionary of target_id -> tpm
        output_file: Path to output file
    """
    rows = []
    header = None
    query_acc_ver_index = None
    rows_with_tpm = 0
    rows_without_tpm = 0

    try:
        # Read BLAST file
        with open(blast_file, 'r') as f:
            for line_num, line in enumerate(f, start=1):
                line = line.rstrip('\n')
                parts = line.split('\t')

                if line_num == 1:
                    # Header line
                    header = parts
                    try:
                        query_acc_ver_index = header.index('query_acc_ver')
                    except ValueError:
                        print(f"Error: 'query_acc_ver' column not found in BLAST file header.",
                              file=sys.stderr)
                        sys.exit(1)
                    rows.append(header)
                else:
                    # Data line
                    if len(parts) != len(header):
                        print(f"Warning: Row {line_num} column count mismatch. Skipping.",
                              file=sys.stderr)
                        continue

                    try:
                        query_id = parts[query_acc_ver_index].strip()
                        tpm = abundance_map.get(query_id, "NA")

                        if tpm != "NA":
                            rows_with_tpm += 1
                        else:
                            rows_without_tpm += 1

                        # Store row with tpm value
                        rows.append((parts, tpm))
                    except (ValueError, IndexError):
                        print(f"Warning: Row {line_num} has invalid data. Skipping.",
                              file=sys.stderr)
                        continue

    except IOError as e:
        print(f"Error reading BLAST file: {e}", file=sys.stderr)
        sys.exit(1)

    # Write output
    try:
        with open(output_file, 'w') as out_f:
            for i, row in enumerate(rows):
                if i == 0:
                    # Header: add 'tpm' column
                    output_line = '\t'.join(row) + '\ttpm\n'
                else:
                    # Data: format with tpm value
                    parts, tpm = row
                    if tpm == "NA":
                        tpm_str = "NA"
                    else:
                        tpm_str = f"{tpm:.6f}"
                    output_line = '\t'.join(parts) + f'\t{tpm_str}\n'

                out_f.write(output_line)

    except IOError as e:
        print(f"Error writing output file: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Merged BLAST and TPM data:", file=sys.stderr)
    print(f"  Rows with TPM values: {rows_with_tpm}", file=sys.stderr)
    print(f"  Rows without TPM values: {rows_without_tpm}", file=sys.stderr)
    print(f"Output saved to: {output_file}")

# ============================================================================
# FUNCTION: Main entry point
# ============================================================================
def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Add TPM (abundance) values to BLAST+taxonomy results.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python3 add_tpm_to_blast.py abundance.tsv blast_taxonomy.txt output_merged.txt

Input files:
  - Abundance file: Kallisto abundance.tsv (must have target_id and tpm columns)
  - BLAST file: Output from add_taxonomy_phyloR.R (must have query_acc_ver column)

Output:
  - Merged file with all BLAST columns + new tpm column (rightmost)
        """
    )

    parser.add_argument("abundance_file",
                        help="Path to Kallisto abundance.tsv file")
    parser.add_argument("blast_taxonomy_file",
                        help="Path to BLAST+taxonomy file (from add_taxonomy_phyloR.R)")
    parser.add_argument("output_file",
                        help="Path to output merged file")

    args = parser.parse_args()

    # Validate input files
    try:
        with open(args.abundance_file, 'r') as f:
            pass
    except FileNotFoundError:
        print(f"Error: Abundance file '{args.abundance_file}' not found.",
              file=sys.stderr)
        sys.exit(1)

    try:
        with open(args.blast_taxonomy_file, 'r') as f:
            pass
    except FileNotFoundError:
        print(f"Error: BLAST file '{args.blast_taxonomy_file}' not found.",
              file=sys.stderr)
        sys.exit(1)

    # Read abundance data
    print(f"Reading abundance data...", file=sys.stderr)
    abundance_map = read_abundance_file(args.abundance_file)
    print(f"Loaded {len(abundance_map)} unique targets with TPM values.",
          file=sys.stderr)

    # Join with BLAST data
    print(f"Reading BLAST data and merging...", file=sys.stderr)
    join_blast_and_abundance(args.blast_taxonomy_file, abundance_map, args.output_file)

if __name__ == "__main__":
    main()
