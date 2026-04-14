#!/usr/bin/env python3

################################################################################
# Script: get_viral_contigs.py
# Description: Extract viral contigs from assembly based on BLAST taxonomy results
# Usage: python3 get_viral_contigs.py --blast BLAST_FILE --fasta CONTIGS_FILE --output OUTPUT_FILE [--min_id MIN_ID]
#
# Inputs:
#   BLAST_FILE:    Path to BLAST tabular output file
#   CONTIGS_FILE:  Path to contig FASTA file
#   OUTPUT_FILE:   Path to output FASTA file with filtered viral contigs
#   MIN_ID:        Minimum percent identity threshold (default: 70)
#
# Outputs:
#   OUTPUT_FILE    FASTA file containing only contigs with BLAST hits at MIN_ID threshold
#
# Logic:
#   - Reads BLAST tabular output
#   - Filters hits by minimum percent identity (column 2)
#   - Extracts corresponding sequence headers
#   - Outputs matching sequences from FASTA file
#   - Preserves sequence order from original FASTA
#
# BLAST tabular columns (format 6):
#   0: query_id, 1: subject_id, 2: pident, 3: length, 4: mismatch,
#   5: gapopen, 6: qstart, 7: qend, 8: sstart, 9: send, 10: evalue, 11: bitscore
#
# Dependencies: Python 3 (only standard library)
#
################################################################################

import argparse
import sys

# ============================================================================
# FUNCTION: Read BLAST tabular output and extract query headers
# ============================================================================
def read_blast_tabular_to_list(blast_file, min_identity):
    """
    Read BLAST tabular output and extract query sequence headers passing filter.

    Args:
        blast_file: Path to BLAST tabular output file
        min_identity: Minimum percent identity threshold

    Returns:
        List of unique query sequence headers passing filter
    """
    query_sequence_headers = []

    try:
        with open(blast_file, 'r') as f:
            for line in f:
                parts = line.strip().split("\t")

                # Skip malformed lines
                if len(parts) < 3:
                    continue

                try:
                    percent_identity = float(parts[2])
                except (ValueError, IndexError):
                    continue

                # Filter by minimum identity
                if percent_identity >= min_identity:
                    query_seq = parts[0]
                    if query_seq not in query_sequence_headers:
                        query_sequence_headers.append(query_seq)

    except IOError as e:
        print(f"Error reading BLAST file: {e}", file=sys.stderr)
        sys.exit(1)

    return query_sequence_headers

# ============================================================================
# FUNCTION: Read FASTA file into dictionary
# ============================================================================
def read_fasta(fasta_file):
    """
    Read FASTA file into dictionary.

    Args:
        fasta_file: Path to input FASTA file

    Returns:
        Dictionary mapping sequence headers to sequences
    """
    sequences = {}
    current_seq_name = None

    try:
        with open(fasta_file, 'r') as f:
            for line in f:
                line = line.strip()

                if line.startswith(">"):
                    # New sequence header (remove leading ">")
                    current_seq_name = line[1:]
                    sequences[current_seq_name] = []

                else:
                    # Sequence data
                    if current_seq_name is not None:
                        sequences[current_seq_name].append(line)

        # Join sequence fragments
        for seq_name, seq_list in sequences.items():
            sequences[seq_name] = ''.join(seq_list)

    except IOError as e:
        print(f"Error reading FASTA file: {e}", file=sys.stderr)
        sys.exit(1)

    return sequences

# ============================================================================
# FUNCTION: Write filtered sequences to output FASTA
# ============================================================================
def write_filtered_fasta(query_headers, sequences, output_file):
    """
    Write filtered sequences to output FASTA file.

    Args:
        query_headers: List of sequence headers to extract
        sequences: Dictionary of sequences
        output_file: Path to output file
    """
    missing_count = 0

    try:
        with open(output_file, 'w') as output_handle:
            for query in query_headers:
                if query in sequences:
                    output_handle.write(f">{query}\n{sequences[query]}\n")
                else:
                    print(f"Warning: Query '{query}' not found in FASTA file.", file=sys.stderr)
                    missing_count += 1

    except IOError as e:
        print(f"Error writing output file: {e}", file=sys.stderr)
        sys.exit(1)

    if missing_count > 0:
        print(f"Warning: {missing_count} queries not found in FASTA file.", file=sys.stderr)

# ============================================================================
# FUNCTION: Main entry point
# ============================================================================
def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Extract viral contigs based on BLAST taxonomy results.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 get_viral_contigs.py --blast results.txt --fasta contigs.fasta --output viral.fasta
  python3 get_viral_contigs.py --blast results.txt --fasta contigs.fasta --output viral.fasta --min_id 60
        """
    )

    parser.add_argument("--blast",
                        required=True,
                        help="Path to BLAST tabular output file")
    parser.add_argument("--fasta",
                        required=True,
                        help="Path to contig FASTA file")
    parser.add_argument("--output",
                        required=True,
                        help="Path to output FASTA file")
    parser.add_argument("--min_id",
                        type=float,
                        default=70,
                        help="Minimum percent identity threshold (default: 70)")

    args = parser.parse_args()

    # Validate input files
    try:
        with open(args.blast, 'r') as f:
            pass
    except FileNotFoundError:
        print(f"Error: BLAST file '{args.blast}' not found.", file=sys.stderr)
        sys.exit(1)

    try:
        with open(args.fasta, 'r') as f:
            pass
    except FileNotFoundError:
        print(f"Error: FASTA file '{args.fasta}' not found.", file=sys.stderr)
        sys.exit(1)

    # Read BLAST output and extract query headers
    print(f"Reading BLAST results (min_id >= {args.min_id}%)...", file=sys.stderr)
    query_sequence_headers = read_blast_tabular_to_list(args.blast, args.min_id)
    print(f"Found {len(query_sequence_headers)} queries passing filter.", file=sys.stderr)

    # Read FASTA file
    print(f"Reading FASTA file...", file=sys.stderr)
    sequences = read_fasta(args.fasta)
    print(f"Loaded {len(sequences)} sequences.", file=sys.stderr)

    # Write filtered output
    print(f"Writing filtered contigs...", file=sys.stderr)
    write_filtered_fasta(query_sequence_headers, sequences, args.output)
    print(f"Viral contigs extracted: {len(query_sequence_headers)} sequences")
    print(f"Output saved to: {args.output}")

if __name__ == "__main__":
    main()
