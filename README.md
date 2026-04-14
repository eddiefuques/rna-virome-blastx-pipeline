# RNA Virome Discovery Pipeline — BLAST-based Approach

A bioinformatics pipeline for the discovery and characterization of RNA viruses from Illumina metagenomics data. This pipeline was used to characterize the RNA virome of *Culex* and *Anopheles* mosquitoes across 37 sites in the Amazonian ecoregion.

> **Associated publication:** Fuques-Vásquez E. *et al.* (2026). Characterization of the RNA virome of *Culex* and *Anopheles* mosquitoes from 37 sites across the Amazonian ecoregion. *PeerJ*.

---

## Overview

This pipeline takes paired-end Illumina FASTQ files from RNA metagenomics experiments and produces a final annotated table of viral contigs with taxonomic information and normalized read abundance. Viral identification is based on protein-level homology to the NCBI RefSeq viral protein database using DIAMOND blastx.

**This pipeline is designed for RNA viromes.** The de novo assembly step uses rnaSPAdes, an assembler optimized for RNA-seq and RNA metagenomics data. For DNA viruses or a combined RNA/DNA approach using multi-tool viral classifiers, see the companion repository: [viral-contig-classification](https://github.com/YOUR_USERNAME/viral-contig-classification).

---

## Pipeline Overview

```
Raw paired-end FASTQ (RNA metagenomics)
    │
    ▼
01  FastQC + MultiQC        Quality assessment of raw reads
    │
    ▼
02  fastp                   Adapter trimming + quality filtering
    │                         avg quality ≥ 25 | base quality ≥ 20 | length ≥ 50 bp
    │                         low-complexity filter ≥ 30% | poly-G/A removal | max 2 Ns
    ▼
    FastQC + MultiQC        Quality assessment of trimmed reads (same script, re-run)
    │
    ▼
03  bbnorm                  k-mer normalization (target depth = 40×)
    │                         Reduces sequencing depth bias before assembly
    ▼
04  rnaSPAdes               De novo assembly (k-mers: 77, 99, 127)
    │                         + length filter ≥ 500 bp
    ▼
05  DIAMOND blastx          Protein-level homology to RefSeq viral proteins
    │                         --sensitive | min-score 300 | %id ≥ 60 | top 3 hits
    ▼
06  phyloR (R)              Retrieve taxonomy for each BLAST hit
    │                         Looks up kingdom, family, genus, and species from NCBI
    │                         for the matched protein in the database
    ▼
07  get_best_hits.py        Keep the single best hit per contig (by bitscore)
    │
    ▼
08  get_viral_contigs.py    Extract viral contigs to FASTA (%id ≥ 70% filter)
    │
    ▼
09  Kallisto                Map trimmed reads back to viral contigs
    │                         Produces TPM, CPM, and relative abundance (PPM)
    ▼
    merge_blast_kallisto.py Final table: BLAST hits + taxonomy + abundance per contig
```

---

## Output

The final output is a tab-separated table with one row per viral contig, containing:

| Column | Description |
|---|---|
| contig_id | Contig identifier from rnaSPAdes |
| subject_acc_ver | Best BLAST hit accession (RefSeq) |
| identity | % amino acid identity to database hit |
| alignment_length | Alignment length (aa) |
| evalue | BLAST E-value |
| bit_score | BLAST bitscore |
| kingdom | Taxonomy of BLAST hit: kingdom |
| family | Taxonomy of BLAST hit: family |
| genus | Taxonomy of BLAST hit: genus |
| species | Taxonomy of BLAST hit: species |
| tpm | Transcripts per million (Kallisto) |
| cpm | Counts per million |
| relative_abundance | Parts per million (reads mapped / total reads × 10⁶) |

**Note on taxonomy:** The `family`, `genus`, and `species` columns describe the taxonomy of the *database match* — retrieved via the phyloR R package from NCBI. This tells you what virus the best BLAST hit belongs to. Formal classification of novel or divergent viruses requires additional phylogenetic analyses (alignment, tree building) not included in this pipeline.

---

## Usage

### 1. Set up the environment

```bash
conda env create -f environment.yml
conda activate viral-metagenomics-p1
```

### 2. Run each step

Each script accepts command-line arguments (see individual script headers for full usage).

```bash
# QC on raw reads
bash 01_qc/run_fastqc_multiqc.sh --input_dir /path/to/raw_fastq --output_dir ./01_qc_raw

# Trimming
bash 02_trimming/run_fastp.sh --input_dir /path/to/raw_fastq --output_dir ./02_trimmed

# QC on trimmed reads (same script, different input)
bash 01_qc/run_fastqc_multiqc.sh --input_dir ./02_trimmed --output_dir ./01_qc_trimmed

# Normalization
bash 03_normalization/run_bbnorm.sh --input_dir ./02_trimmed --output_dir ./03_normalized

# Assembly + length filter
bash 04_assembly/run_rnaspades.sh --input_dir ./03_normalized --output_dir ./04_assembly

# DIAMOND blastx (requires RefSeq viral .dmnd database)
bash 05_diamond_blastx/run_diamond_blastx.sh \
    --contigs_dir ./04_assembly \
    --output_dir  ./05_diamond \
    --db          /path/to/viral.1.protein.faa.dmnd

# Retrieve taxonomy for BLAST hits
Rscript 06_taxonomy/add_taxonomy_phyloR.R ./05_diamond ./06_taxonomy_out

# Best hit per contig
bash 07_filter_blast_results/run_get_best_hits.sh \
    --input_dir ./05_diamond --output_dir ./07_best_hits

# Extract viral contigs to FASTA
bash 08_extract_viral_contigs/run_get_viral_contigs.sh \
    --blast_dir ./07_best_hits --contigs_dir ./04_assembly --output_dir ./08_viral_contigs

# Quantify with Kallisto + calculate CPM + merge final table
bash 09_abundance_kallisto/run_kallisto_pipeline.sh \
    --contigs_dir ./08_viral_contigs --reads_dir ./02_trimmed --output_dir ./09_abundance
python 09_abundance_kallisto/merge_blast_kallisto.py \
    --blast /path/to/blast_taxonomy.csv \
    --kallisto /path/to/kallisto_master.tsv \
    --output final_virome_table.csv
```

---

## Dependencies

| Tool | Version | Purpose |
|---|---|---|
| FastQC | ≥0.11 | Read quality assessment |
| MultiQC | ≥1.14 | Aggregate QC reports |
| fastp | ≥0.23 | Trimming and quality filtering |
| BBTools (bbnorm) | ≥39.0 | k-mer normalization |
| SPAdes (rnaSPAdes) | ≥3.15.3 | RNA de novo assembly |
| DIAMOND | ≥2.1.9 | Protein homology search |
| R + phyloR | ≥4.0 | Taxonomy retrieval for BLAST hits |
| Kallisto | ≥0.48 | Read quantification |
| Python + pandas | ≥3.8 | Data parsing and table merging |

**External database required:** RefSeq viral protein database in DIAMOND format.
```bash
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz
diamond makedb --in viral.1.protein.faa.gz -d viral.1.protein.faa
```

---

## Citation

If you use this pipeline, please cite:

> Fuques-Vásquez E. *et al.* (2026). Characterization of the RNA virome of *Culex* and *Anopheles* mosquitoes from 37 sites across the Amazonian ecoregion. *PeerJ*.

---

## Contact

Eduardo Fuques-Vásquez — eddiefuques@gmail.com
