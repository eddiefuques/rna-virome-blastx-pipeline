#!/usr/bin/env python3

################################################################################
# Script: get_best_hits.py
# Description: Extract best BLAST hit per query based on bitscore
# Usage: python3 get_best_hits.py INPUT_BLAST OUTPUT_FILE
#
# Inputs:
#   INPUT_BLAST:  Path to BLAST tabular output file (format 6)
#   OUTPUT_FILE:  Path to output file for best hits
#
# Outputs:
#   OUTPUT_FILE   Filtered BLAST results (one best hit per query, highest bitscore)
#
# Logic:
#   - For each unique query sequence ID, retains only the hit with highest bitscore
#   - Column 11 (bitscore) is used for ranking
#   - All other columns are preserved from original BLAST output
#
# Dependencies: Python 3 (only standard library)
#
################################################################################

import argparse
import sys

# ============================================================================
# FUNCTION: Parse BLAST output and extract best hits
# ============================================================================
def parse_blast_output(blast_output_file):
    """
    Parse BLAST tabular output and identify best hits per query.

    Args:
        blast_output_file: Path to BLAST tabular output (format 6)

    Returns:
        Dictionary mapping query_id to best hit information (line and bitscore)
    """
    best_hits = {}

    try:
        with open(blast_output_file, 'r') as file:
            for line in file:
                parts = line.strip().split('\t')

                # Skip malformed lines
                if len(parts) < 12:
                    continue

                query_id = parts[0]
                try:
                    bitscore = float(parts[11])
                except (ValueError, IndexError):
                    continue

                # Keep only the hit with highest bitscore for this query
                if query_id not in best_hits or best_hits[query_id]['bitscore'] < bitscore:
                    best_hits[query_id] = {'line': line.strip(), 'bitscore': bitscore}

    except IOError as e:
        print(f"Error reading file: {e}", file=sys.stderr)
        sys.exit(1)

    return best_hits

# ============================================================================
# FUNCTION: Write best hits to output file
# ============================================================================
def write_best_hits(best_hits, output_file):
    """
    Write best hits to output file.

    Args:
        best_hits: Dictionary of best hit information
        output_file: Path to output file
    """
    try:
        with open(output_file, 'w') as file:
            for hit_info in best_hits.values():
                file.write(f"{hit_info['line']}\n")

    except IOError as e:
        print(f"Error writing to file: {e}", file=sys.stderr)
        sys.exit(1)

# ============================================================================
# FUNCTION: Main entry point
# ============================================================================
def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Extract best BLAST hit based on bitscore.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 get_best_hits.py input_blast.txt output_best_hits.txt
  python3 get_best_hits.py results.txt results_best_only.txt
        """
    )

    parser.add_argument("blast_output",
                        help="Path to the BLAST tabular output file (format 6)")
    parser.add_argument("output_file",
                        help="Path to the output file for best hits")

    args = parser.parse_args()

    # Validate input file
    try:
        with open(args.blast_output, 'r') as f:
            pass
    except FileNotFoundError:
        print(f"Error: Input file '{args.blast_output}' not found.", file=sys.stderr)
        sys.exit(1)

    # Extract best hits
    best_hits = parse_blast_output(args.blast_output)

    if not best_hits:
        print("Warning: No valid BLAST hits found in input file.", file=sys.stderr)

    # Write best hits
    write_best_hits(best_hits, args.output_file)
    print(f"Best hits extracted: {len(best_hits)} queries with single best hits")
    print(f"Results saved to: {args.output_file}")

if __name__ == "__main__":
    main()
