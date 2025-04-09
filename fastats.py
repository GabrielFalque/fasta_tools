#!/usr/bin/env python3
"""
Calculates and summarizes nucleotide statistics for sequences in a FASTA file.

Computes the count and percentage of standard nucleotides (A, C, G, T, U),
ambiguous bases (any other letter), and gaps ('-') for each sequence.
Provides summary statistics (mean, min, max) across all sequences.
"""

import argparse
from pathlib import Path
import sys
from Bio import SeqIO # type: ignore # Biopython stubs might be missing
from collections import Counter, defaultdict
import numpy as np # type: ignore # Numpy stubs might be missing
from typing import List, Tuple, Dict, Any # For type hinting

# Define categories for statistics
NUCLEOTIDE_CATS = ["A", "C", "G", "T", "U"]
OTHER_CATS = ["Ambiguous", "Gaps"]
ALL_CATS = NUCLEOTIDE_CATS + OTHER_CATS

# Type alias for statistics dictionary
SeqStats = Dict[str, Any]

def analyze_sequence(record: SeqIO.SeqRecord) -> SeqStats | None:
    """
    Analyzes a single sequence record for nucleotide composition.

    Args:
        record: A Bio.SeqRecord object.

    Returns:
        A dictionary containing counts and percentages for each category,
        and total length, or None if the sequence is empty.
        Keys: 'id', 'length', 'counts', 'percentages'.
        'counts' and 'percentages' are dictionaries mapping category names
        (e.g., 'A', 'Gaps') to their values.
    """
    seq_id = record.id
    seq = str(record.seq).upper()  # Work with uppercase sequence
    seq_length = len(seq)

    if seq_length == 0:
        print(f"Warning: Sequence '{seq_id}' is empty. Skipping.", file=sys.stderr)
        return None

    # Count occurrences of each character
    counts = Counter(seq)

    # Initialize results dictionaries
    stat_counts: Dict[str, int] = defaultdict(int)
    stat_percentages: Dict[str, float] = defaultdict(float)

    # Calculate counts for defined categories
    for nucleotide in NUCLEOTIDE_CATS:
        stat_counts[nucleotide] = counts.get(nucleotide, 0)

    stat_counts["Gaps"] = counts.get("-", 0)

    # Calculate ambiguous bases count (anything not standard nucleotide or gap)
    standard_chars = set(NUCLEOTIDE_CATS + ["-"])
    stat_counts["Ambiguous"] = sum(count for char, count in counts.items() if char not in standard_chars)

    # Calculate percentages
    for category in ALL_CATS:
        stat_percentages[category] = (stat_counts[category] / seq_length) * 100

    return {
        "id": seq_id,
        "length": seq_length,
        "counts": stat_counts,
        "percentages": stat_percentages,
    }


def analyze_fasta_file(fasta_path: Path) -> List[SeqStats]:
    """
    Analyzes all sequences in a FASTA file.

    Args:
        fasta_path: Path object for the input FASTA file.

    Returns:
        A list of dictionaries, where each dictionary contains the statistics
        for one sequence (output from analyze_sequence).

    Raises:
        FileNotFoundError: If the input file does not exist.
        ValueError: If no valid sequences are found in the file.
        Exception: If there's an error parsing the FASTA file.
    """
    if not fasta_path.is_file():
        raise FileNotFoundError(f"Error: FASTA file not found: {fasta_path}")

    all_stats: List[SeqStats] = []
    try:
        for record in SeqIO.parse(str(fasta_path), "fasta"):
            seq_analysis = analyze_sequence(record)
            if seq_analysis:
                all_stats.append(seq_analysis)

    except Exception as e:
        raise Exception(f"Error parsing FASTA file {fasta_path}: {e}") from e

    if not all_stats:
        raise ValueError(f"Error: No valid sequences found or processed in '{fasta_path}'.")

    return all_stats

def print_per_sequence_stats(stats_list: List[SeqStats]) -> None:
    """Prints the statistics for each sequence in a table format (TSV)."""
    # Print header
    header = ["sequence_id", "length"] + [f"perc_{cat}" for cat in ALL_CATS]
    print("\t".join(header))

    # Print data for each sequence
    for stats in stats_list:
        row = [
            str(stats["id"]),
            str(stats["length"]),
        ]
        row.extend([f"{stats['percentages'][cat]:.2f}" for cat in ALL_CATS])
        print("\t".join(row))

def print_summary_statistics(stats_list: List[SeqStats]) -> None:
    """
    Calculates and prints summary statistics (mean, min, max) for percentages.

    Args:
        stats_list: List of sequence statistics dictionaries.
    """
    num_sequences = len(stats_list)
    if num_sequences == 0:
        print("\nNo sequences analyzed, cannot compute summary statistics.")
        return

    # Extract percentage data for numpy calculation
    percentage_data: Dict[str, List[float]] = defaultdict(list)
    length_data: List[int] = []

    for stats in stats_list:
        length_data.append(stats["length"])
        for category in ALL_CATS:
            percentage_data[category].append(stats['percentages'][category])

    print(f"\n--- Summary Statistics ({num_sequences} sequences) ---")

    # Calculate and print summary for length
    lengths = np.array(length_data)
    print(f"Sequence Length: Mean={np.mean(lengths):.2f}, Min={np.min(lengths)}, Max={np.max(lengths)}")

    # Calculate and print summary for percentages
    print("Composition (%):")
    for category in ALL_CATS:
        percentages = np.array(percentage_data[category])
        mean_perc = np.mean(percentages)
        min_perc = np.min(percentages)
        max_perc = np.max(percentages)
        print(f"  {category:>10}: Mean={mean_perc:>6.2f}, Min={min_perc:>6.2f}, Max={max_perc:>6.2f}")

def main() -> None:
    """
    Main function to handle command-line arguments and execute the analysis.
    """
    parser = argparse.ArgumentParser(
        description="Analyze nucleotide composition in a FASTA file. "
                    "Computes percentages of nucleotides (A, C, G, T, U), ambiguous bases, and gaps ('-').",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument(
        "input_fasta", type=Path,
        help="Path to the input FASTA file."
    )
    parser.add_argument(
        "-s", "--summary_only", action="store_true",
        help="Only display the summary statistics, not per-sequence results."
    )

    args = parser.parse_args()

    try:
        # Analyze the FASTA file
        all_sequence_stats = analyze_fasta_file(args.input_fasta)

        # Print results
        if not args.summary_only:
            print_per_sequence_stats(all_sequence_stats)

        print_summary_statistics(all_sequence_stats)

    except (FileNotFoundError, ValueError, Exception) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()