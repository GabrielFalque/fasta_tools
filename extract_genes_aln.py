#!/usr/bin/env python3

"""
Script to parse command-line arguments, process an alignment file, and
generate separate gene alignment files based on a reference annotation file.

Arguments are provided via the command line and include paths to the reference file,
alignment file, and output directory, as well as options for verbosity.

Usage :
./extract_genes_aln.py --ref_file complete_NC_045512.2_annotated_genes.txt \
--aligned_file complete_genome_mafft_aln.fasta \
--seq_id NC_045512.2 \
--output_directory /output/directory
"""

import warnings
warnings.filterwarnings("ignore")
import argparse
import io
import os
import sys
from Bio import SeqIO, Align,SearchIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Align import MultipleSeqAlignment
from Bio.Align import PairwiseAligner
from Bio.SeqRecord import SeqRecord
import subprocess
import re
from Bio import SeqIO

def main():

    parser = argparse.ArgumentParser(
        description=" ".join([
            "Script for separating an alignement file into multiple gene alignement files based on",
            "reference annotated file from NCBI dataset."
        ]
        )
    )
    # Arguments definition
    parser.add_argument(
        "-r", "--ref_file",
        type=str,
        required=True,
        help=" ".join([
            "Text (.txt) file containing annotated reference sequence with genes. This type of file can be",
            "downladed on NCBI database."
        ])
    )

    parser.add_argument(
        "-a", "--aligned_file",
        type=str,
        required=True,
        help="FASTA file containing aligned sequences."
    )

    parser.add_argument(
        "-i", "--seq_id",
        type=str,
        required=True,
        help="Reference sequence in alignement file. Please keep exact same name as in alignement file."
    )

    parser.add_argument(
        "-o", "--output_directory",
        type=str,
        required=True,
        help="Output directory to save results file."
    )

    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Verbose mode."
    )

    # Arguments analyse
    args = parser.parse_args()

    # Verify and create if directory does not exist
    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)
        if args.verbose:
            print(f"Directory '{args.output_directory}' created.")
    else:
        if args.verbose:
            print(f"Output directory '{args.output_directory}' already exists.")

    # Display options
    if args.verbose:
        print("Reference file:", args.ref_file)
        print("Alignement file:", args.aligned_file)
        print("Reference sequence ID:", args.seq_id)


    # Load sequences alignment file
    alignment = list(SeqIO.parse(args.aligned_file, "fasta"))

    # Search for reference sequence in alignement
    ref_seq = next((record for record in alignment if record.id == args.seq_id), None)

    if ref_seq is None:
        print(f"Error : sequence ID '{args.seq_id}' was not found in alignment file.")
        return

    print(f"{len(alignment)} aligned sequences loaded.")


    def parse_gene_info(ref_record):
        """
        Parses gene information from a reference annotation file.

        This function extracts gene names, their start and end positions, 
        and their sequences from the provided reference records. It assumes
        that the reference records contain annotation fields like `[gene=...]`
        and `[location=...]`.

        Args:
            ref_record (list of SeqRecord): A list of `SeqRecord` objects representing 
                the reference sequences with annotations.

        Returns:
            list of dict: A list of dictionaries where each dictionary contains:
                - "gene_name" (str): The name of the gene.
                - "start" (int): The start position of the gene (1-based).
                - "end" (int): The end position of the gene (1-based).
                - "sequence" (str): The nucleotide sequence of the gene.
        """

        genes_info = []

        for record in ref_record:
            description = record.description
            sequence = str(record.seq)

            # Extract name gene
            gene_name_match = re.search(r'\[gene=(.*?)\]', description)
            gene_name = gene_name_match.group(1) if gene_name_match else None

            # Extract gene location
            location_match = re.search(r'\[location=(\d+)\.\.(\d+)\]', description)
            if location_match:
                start = int(location_match.group(1))
                end = int(location_match.group(2))
            else:
                start, end = None, None

            # Add info to a list
            genes_info.append({
                "gene_name": gene_name,
                "start": start,
                "end": end,
                "sequence": sequence
            })

        return genes_info

    def get_coordinates_in_alignment(genes_info, ref_sequence_record):
        """
        Locate gene sequences within the aligned reference sequence and retrieve their coordinates.

        Args:
            genes_info (list): A list of dictionaries containing gene information (name, location, sequence).
            ref_sequence_record (SeqRecord): The reference sequence record from the alignment file.

        Returns:
            list: A list of dictionaries containing gene name and alignment coordinates (start and end positions).
        """
        seq_ref = str(ref_sequence_record.seq)
        aligned_genes = []

        for gene in genes_info:
            gene_seq = gene["sequence"].lower()
            gene_name = gene["gene_name"]
            target_index = 0  # Index in target_seq
            start = None
            end = None

            # Loop for each position in seq_ref
            for i, base in enumerate(seq_ref):
                # Verify if present nucleotide in seq_ref corresponds to present nucleotide in gene_seq.
                if base == gene_seq[target_index]:
                    if start is None:
                        if seq_ref[i:].replace("-","")[:len(gene_seq)]==gene_seq:
                            start = i + 1  # Position in 1-base
                            target_index += 1
                    else:
                        # Increment target index to move along gene_seq
                        target_index += 1
            
                        # If all target sequence is found
                        if target_index == len(gene_seq):
                            end = i + 1  # Position in 1-base
                            break
            
            print(f"Best ungapped alignment for {gene_name} found at position {start}-{end} ")

            # Saving gene information
            aligned_genes.append({
                "gene_name": gene_name,
                "start_in_ref": start,
                "end_in_ref": end
            })

        return aligned_genes

    ref_record = [record for record in SeqIO.parse(args.ref_file, "fasta")]
    genes_info = parse_gene_info(ref_record)

    # Extract aligned reference sequence.
    ref_sequence_record = None
    ref_sequence = None
    for record in alignment:
        if record.id == args.seq_id:  
            ref_sequence_record = record 
            ref_sequence = record.seq
            break

    aligned_genes_info = get_coordinates_in_alignment(genes_info, ref_sequence_record)

    # Using MAFFT to align every gene sequences
    for gene_info in aligned_genes_info:
        gene = gene_info["gene_name"]
        start = gene_info["start_in_ref"]
        end = gene_info["end_in_ref"]
        # Create file for gene to align
        aligned_gene_file = args.output_directory + f"/gene_{gene}.fasta"
        with open(aligned_gene_file, "w") as temp_file:

            # Add aligned sequences from gene to the file
            for record in alignment:
                if record.id != "seq_ref":
                    gene_seq_aligned = record.seq[start-1:end]
                    temp_file.write(f">{record.id}\n{gene_seq_aligned}\n")


    print("Genes alignement completed and FASTA files generated at ", args.output_directory)


if __name__ == "__main__":
    main()