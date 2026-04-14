#!/usr/bin/env Rscript

################################################################################
# Script: add_taxonomy_phyloR.R
# Description: Annotate BLAST tabular output with taxonomic information using phyloR
# Usage: Rscript add_taxonomy_phyloR.R INPUT_DIR OUTPUT_DIR
#
# Inputs:
#   INPUT_DIR:   Directory containing *_diamond_tabular_out.txt files
#   OUTPUT_DIR:  Directory for taxonomic annotated output
#
# Outputs:
#   {SAMPLE}_blast_and_taxonomy_diamond_results.txt  - Full results with taxonomy for the BLAST hits
#   {SAMPLE}_taxonomy_diamond_results.txt            - Condensed taxonomy table
#
# Taxonomy levels added:
#   - kingdom, family, genus, species
#
# Output columns (full):
#   All BLAST columns + kingdom + family + genus + species
#
# Output columns (condensed):
#   query_acc_ver, subject_acc_ver, kingdom, family, genus, species
#
# Dependencies: phyloR, dplyr, tidyverse, tibble
#
# Notes:
#   - Requires phyloR package to be installed
#   - NCBI taxonomy database must be accessible to phyloR
#   - Run times may be long depending on number of hits
#
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(phyloR, quietly = TRUE)
  library(dplyr, quietly = TRUE)
  library(tidyverse, quietly = TRUE)
  library(tibble, quietly = TRUE)
})

# ============================================================================
# FUNCTION: Print usage
# ============================================================================
print_usage <- function() {
  cat("Usage: Rscript add_taxonomy_phyloR.R INPUT_DIR OUTPUT_DIR\n\n")
  cat("  INPUT_DIR   Directory containing *_diamond_tabular_out.txt files\n")
  cat("  OUTPUT_DIR  Directory for annotated output\n\n")
  quit(status = 1)
}

# ============================================================================
# Parse command-line arguments
# ============================================================================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Error: INPUT_DIR and OUTPUT_DIR are required.\n\n")
  print_usage()
}

input_dir <- args[1]
output_dir <- args[2]

# Validate directories
if (!dir.exists(input_dir)) {
  cat("Error: INPUT_DIR '", input_dir, "' does not exist.\n", sep = "")
  quit(status = 1)
}

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Get list of BLAST output files
# ============================================================================
files <- list.files(path = input_dir,
                     pattern = "_diamond_tabular_out\\.txt$",
                     full.names = TRUE)

if (length(files) == 0) {
  cat("Warning: No files matching pattern *_diamond_tabular_out.txt found in ",
      input_dir, "\n", sep = "")
  quit(status = 1)
}

cat("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "] Starting taxonomy annotation for ", length(files), " file(s)...\n", sep = "")

# ============================================================================
# Process each BLAST output file
# ============================================================================
for (i in 1:length(files)) {
  current_file <- files[i]
  cat("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      "] Processing file ", i, "/", length(files), ": ",
      basename(current_file), "\n", sep = "")

  # Read BLAST output table
  tryCatch({
    blast_out_tbl <- read.table(current_file, header = FALSE,
                                 stringsAsFactors = FALSE)

    # Add empty column for taxonomy
    blast_out_tbl <- blast_out_tbl %>% add_column(new_column = NA)

    # Convert to tibble
    blast_out_tbl <- as_tibble(blast_out_tbl)

    # Add column names using phyloR function
    colnames(blast_out_tbl) <- phyloR::get_blast_outformat_7_colnames()

    # Add taxonomy progressively: kingdom -> family -> genus -> species
    cat("    Adding kingdom annotation...\n")
    with_kingdom <- blast_out_tbl %>%
      phyloR::add_taxonomy_columns(ncbi_accession_colname = "subject_acc_ver",
                                    taxonomy_level = "kingdom")

    cat("    Adding family annotation...\n")
    with_kingdom_family <- with_kingdom %>%
      phyloR::add_taxonomy_columns(ncbi_accession_colname = "subject_acc_ver",
                                    taxonomy_level = "family")

    cat("    Adding genus annotation...\n")
    with_kingdom_family_genus <- with_kingdom_family %>%
      phyloR::add_taxonomy_columns(ncbi_accession_colname = "subject_acc_ver",
                                    taxonomy_level = "genus")

    cat("    Adding species annotation...\n")
    with_kingdom_family_genus_species <- with_kingdom_family_genus %>%
      phyloR::add_taxonomy_columns(ncbi_accession_colname = "subject_acc_ver",
                                    taxonomy_level = "species")

    # Create condensed version (select key columns)
    condensed <- with_kingdom_family_genus_species %>%
      dplyr::select(query_acc_ver, subject_acc_ver, kingdom, family, genus, species)

    # Extract sample name from filename
    file_base <- basename(current_file)
    sample_name <- head(strsplit(file_base, "_")[[1]], n = 1)

    # Write full annotation table
    full_output <- file.path(output_dir,
                             paste(sample_name, "_blast_and_taxonomy_diamond_results.txt", sep = ""))
    write.table(with_kingdom_family_genus_species, file = full_output,
                sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    cat("    Wrote full results: ", basename(full_output), "\n", sep = "")

    # Write condensed taxonomy table
    condensed_output <- file.path(output_dir,
                                  paste(sample_name, "_taxonomy_diamond_results.txt", sep = ""))
    write.table(condensed, file = condensed_output,
                sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    cat("    Wrote condensed results: ", basename(condensed_output), "\n", sep = "")

  }, error = function(e) {
    cat("    Error processing file: ", conditionMessage(e), "\n", sep = "")
  })
}

cat("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "] Taxonomy annotation complete.\n", sep = "")
