# fasta_tools

This repo contains some python scripts in order to manipulate FASTA files.

## extract_genes_aln.py

Script to parse command-line arguments, process an alignment file, and generate separate gene alignment files based on a reference annotation file.

Arguments are provided via the command line and include paths to the reference file, alignment file, and output directory, as well as options for verbosity.

Usage :
```bash
python ./extract_genes_aln.py --ref_file complete_NC_045512.2_annotated_genes.txt \
--aligned_file complete_genome_mafft_aln.fasta \
--seq_id NC_045512.2 \
--output_directory /output/directory
```

