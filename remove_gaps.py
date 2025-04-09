#!/usr/bin/env python3
"""
Removes gap characters ('-') from sequences in FASTA or plain text files.

Reads sequences from an input file, removes all occurrences of '-',
and writes the cleaned sequences to an output file.
Input file format is determined by the file extension (.fasta or .fa).
Files with other extensions are treated as plain text (one sequence per line).
"""

import argparse
import os
from pathlib import Path  # Use pathlib for robust path handling
from Bio import SeqIO
import sys  # For stderr

def remove_gaps_fasta(input_path: Path, output_path: Path) -> None:
    """
    Removes gaps from sequences in a FASTA file.

    Reads sequences using BioPython, removes '-' characters from each sequence,
    and writes the modified records to the output FASTA file.

    Args:
        input_path: Path object for the input FASTA file.
        output_path: Path object for the output FASTA file.

    Raises:
        FileNotFoundError: If the input file does not exist.
        Exception: If there's an error parsing the FASTA file.
    """
    if not input_path.is_file():
        raise FileNotFoundError(f"Error: Input FASTA file not found: {input_path}")

    try:
        with output_path.open("w") as output_handle:
            # Parse the input FASTA file
            for record in SeqIO.parse(str(input_path), "fasta"):
                # Remove gap characters from the sequence
                cleaned_seq = record.seq.replace("-", "")
                if cleaned_seq: # Only write if the sequence is not empty after cleaning
                    record.seq = cleaned_seq
                    SeqIO.write(record, output_handle, "fasta")
                else:
                    print(f"Warning: Sequence '{record.id}' became empty after removing gaps. Skipping.", file=sys.stderr)
    except Exception as e:
        raise Exception(f"Error processing FASTA file {input_path}: {e}") from e

def remove_gaps_txt(input_path: Path, output_path: Path) -> None:
    """
    Removes gaps from sequences in a plain text file.

    Assumes one sequence per line. Reads each line, removes '-' characters
    and leading/trailing whitespace, and writes the cleaned sequence
    to the output file if it's not empty.

    Args:
        input_path: Path object for the input TXT file.
        output_path: Path object for the output TXT file.

    Raises:
        FileNotFoundError: If the input file does not exist.
        IOError: If there's an error reading or writing files.
    """
    if not input_path.is_file():
        raise FileNotFoundError(f"Error: Input TXT file not found: {input_path}")

    try:
        with input_path.open("r") as input_handle, output_path.open("w") as output_handle:
            for line in input_handle:
                # Remove gaps and strip whitespace
                cleaned_line = line.strip().replace("-", "")
                if cleaned_line:  # Avoid writing empty lines
                    output_handle.write(cleaned_line + "\n")
    except IOError as e:
        raise IOError(f"Error processing text file {input_path}: {e}") from e

def main() -> None:
    """
    Main function to parse command-line arguments and process the file.
    """
    parser = argparse.ArgumentParser(
        description="Remove gap characters ('-') from sequences in a FASTA or TXT file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter # Show default values
    )
    parser.add_argument(
        "-i", "--input", required=True, type=Path,
        help="Input file path (FASTA or TXT format)."
    )
    parser.add_argument(
        "-o", "--output", required=True, type=Path,
        help="Output file path."
    )
    args = parser.parse_args()

    input_file: Path = args.input
    output_file: Path = args.output

    # Determine file format based on extension
    is_fasta = input_file.suffix.lower() in [".fasta", ".fa", ".fna", ".faa"] # Common FASTA extensions

    try:
        if is_fasta:
            print(f"Processing '{input_file}' as FASTA format...")
            remove_gaps_fasta(input_file, output_file)
        else:
            print(f"Processing '{input_file}' as plain text format (one sequence per line)...")
            remove_gaps_txt(input_file, output_file)

        print(f"Successfully removed gaps. Output saved to '{output_file}'.")

    except (FileNotFoundError, IOError, Exception) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1) # Exit with a non-zero status to indicate failure

if __name__ == "__main__":
    main()