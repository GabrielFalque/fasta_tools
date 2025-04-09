#!/usr/bin/env python3
"""
Translates DNA sequences from a FASTA file into protein sequences.

Reads DNA sequences, translates them using the specified reading frame,
and writes the resulting protein sequences to an output FASTA file.
Translation stops at the first in-frame stop codon.
"""

import argparse
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def translate_dna_fasta(input_path: Path, output_path: Path, frame: int = 1) -> int:
    """
    Translates DNA sequences from an input FASTA file to protein sequences.

    Args:
        input_path: Path object for the input DNA FASTA file.
        output_path: Path object for the output protein FASTA file.
        frame: Reading frame for translation (1, 2, or 3).

    Returns:
        The number of sequences successfully processed and written.

    Raises:
        FileNotFoundError: If the input file does not exist.
        ValueError: If the frame number is invalid or sequence length is
                    insufficient for the chosen frame.
        Exception: If there's an error parsing the FASTA file or during translation.
    """
    if not input_path.is_file():
        raise FileNotFoundError(f"Error: Input FASTA file not found: {input_path}")
    if frame not in [1, 2, 3]:
        raise ValueError(f"Error: Invalid frame '{frame}'. Must be 1, 2, or 3.")

    processed_count = 0
    try:
        with output_path.open("w") as output_handle:
            for record in SeqIO.parse(str(input_path), "fasta"):
                dna_seq = record.seq
                seq_len = len(dna_seq)

                # Adjust sequence start based on frame
                if frame > seq_len:
                     print(f"Warning: Sequence '{record.id}' (length {seq_len}) is too short for frame {frame}. Skipping.", file=sys.stderr)
                     continue

                # Select the subsequence for the correct frame (0-based index)
                framed_dna_seq = dna_seq[frame - 1:]

                # Check if the remaining sequence is a multiple of 3 (needed for translation)
                # Biopython's translate handles this, but warn if len % 3 != 0
                if len(framed_dna_seq) % 3 != 0:
                     print(f"Warning: Length of sequence '{record.id}' starting from frame {frame} ({len(framed_dna_seq)}) is not a multiple of 3. Translation might be partial.", file=sys.stderr)

                try:
                    # Translate to protein, stopping at the first stop codon
                    protein_seq: Seq = framed_dna_seq.translate(to_stop=True)
                except Exception as trans_err: # Catch potential translation errors
                    print(f"Warning: Could not translate sequence '{record.id}' (frame {frame}): {trans_err}. Skipping.", file=sys.stderr)
                    continue

                if not protein_seq:
                     print(f"Warning: Translation of sequence '{record.id}' (frame {frame}) resulted in an empty protein sequence (e.g., immediate stop codon or non-coding). Skipping output.", file=sys.stderr)
                     continue

                # Create a new SeqRecord for the protein sequence
                protein_record = SeqRecord(
                    protein_seq,
                    id=record.id,
                    # Keep original description and add translation info
                    description=f"{record.description} | Translated Frame {frame}"
                )
                SeqIO.write(protein_record, output_handle, "fasta")
                processed_count += 1

    except Exception as e:
        # Catch file parsing errors or other unexpected issues
        raise Exception(f"Error processing file {input_path}: {e}") from e

    if processed_count == 0:
         print(f"Warning: No sequences were successfully translated and written to '{output_path}'.", file=sys.stderr)


    return processed_count

def main() -> None:
    """
    Main function to parse command-line arguments and run the translation.
    """
    parser = argparse.ArgumentParser(
        description="Translate DNA sequences in a FASTA file into protein sequences, stopping at the first stop codon.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument(
        "-i", "--input", required=True, type=Path,
        help="Input FASTA file containing DNA sequences."
    )
    parser.add_argument(
        "-o", "--output", required=True, type=Path,
        help="Output FASTA file for the resulting protein sequences."
    )
    parser.add_argument(
        "-f", "--frame", type=int, choices=[1, 2, 3], default=1,
        help="Reading frame for translation (1, 2, or 3)."
    )
    args = parser.parse_args()

    try:
        num_translated = translate_dna_fasta(args.input, args.output, args.frame)
        print(f"Translation completed. {num_translated} protein sequences saved to '{args.output}'.")

    except (FileNotFoundError, ValueError, Exception) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()