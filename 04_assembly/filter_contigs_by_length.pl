#!/usr/bin/perl

################################################################################
# Script: filter_contigs_by_length.pl
# Description: Filter FASTA sequences by minimum length threshold
# Usage: filter_contigs_by_length.pl MIN_LENGTH input_file.fasta > output_file.fasta
#
# Inputs:
#   MIN_LENGTH      Minimum sequence length threshold (in base pairs)
#   input_file.fasta  Input FASTA file
#
# Outputs:
#   stdout          Filtered FASTA sequences (redirect with > output_file.fasta)
#
# Notes:
#   - Filters out sequences shorter than MIN_LENGTH
#   - Preserves sequence headers exactly as written in input
#   - Output is written to stdout (redirect as needed)
#
# Dependencies: Perl core modules (none beyond standard library)
################################################################################

use strict;
use warnings;

# ============================================================================
# Parse command-line arguments
# ============================================================================
my $minlen = shift or die "Error: MIN_LENGTH parameter not provided\n";

die "Error: MIN_LENGTH must be a positive integer\n" unless $minlen =~ /^\d+$/ && $minlen > 0;

if (!@ARGV) {
    die "Error: Input file not provided\n";
}

# ============================================================================
# Filter FASTA sequences by length
# ============================================================================
{
    local $/=">";  # Set record separator to ">"
    while (<>) {
        chomp;
        next unless /\w/;  # Skip empty records
        s/>$//gs;          # Remove trailing ">"

        my @chunk = split /\n/;
        my $header = shift @chunk;
        my $seqlen = length join "", @chunk;

        # Output sequence if it meets minimum length threshold
        print ">$_" if ($seqlen >= $minlen);
    }
    local $/="\n";  # Reset record separator
}

exit 0;
