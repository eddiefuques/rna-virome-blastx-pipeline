#!/usr/bin/env python3

################################################################################
# Script: merge_blast_kallisto.py
# Description: Merge BLAST+taxonomy results with Kallisto master abundance table
# Usage: python3 merge_blast_kallisto.py --blast BLAST_FILE --kallisto KALLISTO_FILE --output OUTPUT_FILE
#
# Inputs:
#   BLAST_FILE:     BLAST results with taxonomy and TPM (from add_tpm_to_blast.py)
#   KALLISTO_FILE:  Master Kallisto abundance table with all samples
#   OUTPUT_FILE:    Output file with merged results
#
# Outputs:
#   OUTPUT_FILE  Combined table with columns:
#     query_acc_ver, subject_acc_ver, pident, ..., kingdom, family, genus, species,
#     tpm, est_counts, length, effective_length
#
# Logic:
#   - Reads BLAST+taxonomy+TPM file (query_acc_ver = contig_id)
#   - Reads Kallisto master table (target_id = contig_id)
#   - Joins on query_acc_ver (BLAST) = target_id (Kallisto)
#   - Combines all columns from both files
#
# Notes:
#   - BLAST file comes from add_tpm_to_blast.py
#   - Kallisto file should be combined sample abundance table
#   - Preserves all original columns from both files
#   - Unmatched rows are kept with NA for missing columns
#
# Dependencies: Python 3 (only standard library)
#
################################################################################

import argparse
import sys

# ============================================================================
# FUNCTION: Read Kallisto master abundance file
# ============================================================================
def read_kallisto_file(kallisto_file):
    """
    Read Kallisto master abundance file and create target_id mapping.

    Args:
        kallisto_file: Path to Kallisto master abundance file

    Returns:
        Tuple of (header list, dictionary mapping target_id to data)
    """
    kallisto_data = {}
    header = None

    try:
        with open(kallisto_file, 'r') as f:
            for line_num, line in enumerate(f, start=1):
                line = line.rstrip('\n')
                parts = line.split('\t')

                if line_num == 1:
                    # Header line
                    header = parts
                    try:
                        target_id_index = header.index('target_id')
                    except ValueError:
                        print(f"Error: 'target_id' column not found in Kallisto header.",
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
                        kallisto_data[target_id] = parts
                    except (ValueError, IndexError):
                        print(f"Warning: Row {line_num} has invalid data. Skipping.",
                              file=sys.stderr)
                        continue

    except IOError as e:
        print(f"Error reading Kallisto file: {e}", file=sys.stderr)
        sys.exit(1)

    return header, kallisto_data

# ============================================================================
# FUNCTION: Merge BLAST and Kallisto data
# ============================================================================
def merge_blast_and_kallisto(blast_file, kallisto_header, kallisto_data, output_file):
    """
    Read BLAST file, add Kallisto data, and write merged output.

    Args:
        blast_file: Path to BLAST+taxonomy+TPM file
        kallisto_header: Header from Kallisto file
        kallisto_data: Dictionary of target_id -> kallisto data
        output_file: Path to output file
    """
    blast_header = None
    query_acc_ver_index = None
    rows_with_kallisto = 0
    rows_without_kallisto = 0

    try:
        # Read BLAST file and write merged output
        with open(blast_file, 'r') as in_f, open(output_file, 'w') as out_f:
            for line_num, line in enumerate(in_f, start=1):
                line = line.rstrip('\n')
                parts = line.split('\t')

                if line_num == 1:
                    # Header line
                    blast_header = parts
                    try:
                        query_acc_ver_index = blast_header.index('query_acc_ver')
                    except ValueError:
                        print(f"Error: 'query_acc_ver' column not found in BLAST header.",
                              file=sys.stderr)
                        sys.exit(1)

                    # Write combined header
                    # Exclude 'target_id' from Kallisto header to avoid duplication
                    kallisto_cols = [col for col in kallisto_header
                                     if col != 'target_id']
                    combined_header = blast_header + kallisto_cols
                    out_f.write('\t'.join(combined_header) + '\n')
                else:
                    # Data line
                    if len(parts) != len(blast_header):
                        print(f"Warning: Row {line_num} column count mismatch. Skipping.",
                              file=sys.stderr)
                        continue

                    try:
                        query_id = parts[query_acc_ver_index].strip()
                        kallisto_row = kallisto_data.get(query_id)

                        if kallisto_row is not None:
                            rows_with_kallisto += 1
                            # Exclude target_id column from Kallisto data
                            target_id_index = kallisto_header.index('target_id')
                            kallisto_cols = [kallisto_row[i] for i in range(len(kallisto_row))
                                             if i != target_id_index]
                        else:
                            rows_without_kallisto += 1
                            # Create NA values for missing Kallisto columns
                            target_id_index = kallisto_header.index('target_id')
                            kallisto_cols = ['NA' for col in kallisto_header
                                             if col != 'target_id']

                        # Write merged row
                        merged_row = parts + kallisto_cols
                        out_f.write('\t'.join(merged_row) + '\n')

                    except (ValueError, IndexError) as e:
                        print(f"Warning: Row {line_num} has invalid data: {e}. Skipping.",
                              file=sys.stderr)
                        continue

    except IOError as e:
        print(f"Error during merge: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Merged BLAST and Kallisto data:", file=sys.stderr)
    print(f"  Rows with Kallisto data: {rows_with_kallisto}", file=sys.stderr)
    print(f"  Rows without Kallisto data: {rows_without_kallisto}", file=sys.stderr)
    print(f"Output saved to: {output_file}")

# ============================================================================
# FUNCTION: Main entry point
# ============================================================================
def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Merge BLAST+taxonomy results with Kallisto master abundance table.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python3 merge_blast_kallisto.py --blast blast_tpm.txt --kallisto master_abundance.tsv --output merged_final.tsv

Input files:
  - BLAST file: Output from add_tpm_to_blast.py (must have query_acc_ver)
  - Kallisto file: Master abundance table (must have target_id column)

Output:
  - Merged file with all columns from both files, joined on contig ID
        """
    )

    parser.add_argument("--blast",
                        required=True,
                        help="Path to BLAST+taxonomy+TPM file")
    parser.add_argument("--kallisto",
                        required=True,
                        help="Path to Kallisto master abundance table")
    parser.add_argument("--output",
                        required=True,
                        help="Path to output merged file")

    args = parser.parse_args()

    # Validate input files
    try:
        with open(args.blast, 'r') as f:
            pass
    except FileNotFoundError:
        print(f"Error: BLAST file '{args.blast}' not found.",
              file=sys.stderr)
        sys.exit(1)

    try:
        with open(args.kallisto, 'r') as f:
            pass
    except FileNotFoundError:
        print(f"Error: Kallisto file '{args.kallisto}' not found.",
              file=sys.stderr)
        sys.exit(1)

    # Read Kallisto data
    print(f"Reading Kallisto master abundance table...", file=sys.stderr)
    kallisto_header, kallisto_data = read_kallisto_file(args.kallisto)
    print(f"Loaded {len(kallisto_data)} targets from Kallisto file.",
          file=sys.stderr)

    # Merge with BLAST data
    print(f"Reading BLAST data and merging...", file=sys.stderr)
    merge_blast_and_kallisto(args.blast, kallisto_header, kallisto_data, args.output)

if __name__ == "__main__":
    main()
