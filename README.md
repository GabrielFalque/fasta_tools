# FASTA Tools

A collection of Python scripts for common bioinformatics tasks involving FASTA files.

## Requirements

* Python 3.8+
* Biopython
* Numpy

You can install the required libraries using pip:
```bash
pip install biopython numpy
```

## Scripts

### 1. `remove_gaps.py`

Removes gap characters (`-`) from sequences in FASTA or plain text files.

**Purpose:** Cleans alignment files or sequences containing gaps, producing sequences with only standard characters remaining. Reads FASTA (`.fasta`, `.fa`, etc.) or plain text (one sequence per line) files, identifies format by extension, removes all `-` characters, and writes the cleaned sequences to a new file. Empty sequences resulting from gap removal are skipped.

**Usage:**

```bash
python remove_gaps.py -i <input_file> -o <output_file>
```

**Arguments:**

* `-i`, `--input` (required): Path to the input file (FASTA or TXT format).
* `-o`, `--output` (required): Path to the output file where cleaned sequences will be saved.

**Example:**

```bash
python remove_gaps.py --input aligned_sequences.fasta --output sequences_no_gaps.fasta
```

---

### 2. `dna_to_prot.py`

Translates DNA sequences from a FASTA file into protein sequences.

**Purpose:** Converts coding DNA sequences into their corresponding protein sequences using the standard genetic code. The translation starts at the specified reading frame and stops at the first encountered in-frame stop codon.

**Usage:**

```bash
python dna_to_prot.py -i <input_dna.fasta> -o <output_protein.fasta> [-f <frame>]
```

**Arguments:**

* `-i`, `--input` (required): Path to the input FASTA file containing DNA sequences.
* `-o`, `--output` (required): Path to the output FASTA file for the resulting protein sequences.
* `-f`, `--frame` (optional): Reading frame for translation (1, 2, or 3). Default is `1`.

**Example:**

```bash
# Translate using frame 1 (default)
python dna_to_prot.py -i coding_regions.fasta -o proteins.fasta

# Translate using frame 2
python dna_to_prot.py --input coding_regions.fasta --output proteins_frame2.fasta --frame 2
```

---

### 3. `fastats.py`

Calculates and summarizes nucleotide/gap/ambiguous base statistics for sequences in a FASTA file.

**Purpose:** Provides insights into the composition of sequences within a FASTA file. For each sequence, it calculates the percentage of A, C, G, T, U, ambiguous bases (any other letter), and gaps (`-`). It then provides summary statistics (mean, min, max) for sequence length and percentage composition across all sequences in the file.

**Usage:**

```bash
python fastats.py <input_fasta> [-s]
```

**Arguments:**

* `input_fasta` (required): Path to the input FASTA file.
* `-s`, `--summary_only` (optional): If specified, only display the summary statistics, omitting the per-sequence percentage table.

**Output:**

Prints a tab-separated table (unless `-s` is used) with percentage composition for each sequence, followed by summary statistics (mean, min, max for length and percentages) printed to the console.

**Example:**

```bash
# Calculate and display per-sequence and summary stats
python fastats.py my_sequences.fasta

# Calculate and display only summary stats
python fastats.py my_sequences.fasta --summary_only
```

---

### 4. `extract_genes_aln.py`

Extracts individual gene alignments from a whole genome multiple sequence alignment (MSA) based on reference annotations.

**Purpose:** Useful for phylogenetic analysis or detailed study of specific genes across multiple aligned genomes. It takes a whole genome alignment, a reference annotation file (FASTA format with location tags), and the ID of the reference sequence within the alignment. It maps the gene coordinates from the annotations onto the gapped reference, extracts the corresponding alignment columns for all genomes, and saves each gene's alignment to a separate FASTA file.

**Annotation File Format:** The script expects a multi-FASTA file where headers contain gene identification and location information, similar to GenBank downloads. Specifically, it looks for:
* `[gene=GENE_NAME]` or `[locus_tag=LOCUS_TAG]` for the gene identifier.
* `[location=START..END]` or `[location=complement(START..END)]` for the coordinates (1-based, inclusive). It handles basic `<` and `>` boundary markers and takes the outermost coordinates for `join` or `order` locations (this may be inaccurate for complex features like spliced genes).

**Usage:**

```bash
python extract_genes_aln.py --annotations <ref_annotations.fasta> \
                            --alignment <whole_genome_aln.fasta> \
                            --ref_id <reference_sequence_id> \
                            --output_dir <directory_for_gene_alignments>
```

**Arguments:**

* `-a`, `--annotations` (required): Path to the reference gene annotation file (multi-FASTA format with specific header tags).
* `-g`, `--alignment` (required): Path to the whole genome alignment file (FASTA format).
* `-r`, `--ref_id` (required): Sequence ID (header ID) of the reference genome within the alignment file. This *must* match exactly the ID used in the alignment file (usually the part before the first space).
* `-o`, `--output_dir` (required): Path to the directory where the output FASTA files for each extracted gene alignment will be saved.

**Example:**

```bash
python extract_genes_aln.py \
    --annotations NC_045512.2_genes.fasta \
    --alignment covid_variants_aligned.fasta \
    --ref_id NC_045512.2 \
    --output_dir extracted_gene_alignments
```

This will create files like `extracted_gene_alignments/gene_S.fasta`, `extracted_gene_alignments/gene_ORF1ab.fasta`, etc., each containing the alignment slice corresponding to that gene across all sequences from `covid_variants_aligned.fasta`.
```

---
