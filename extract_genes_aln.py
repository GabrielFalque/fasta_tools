#!/usr/bin/env python3
"""
Extracts individual gene alignments from a whole genome multiple sequence alignment (MSA).

Parses gene coordinates (start, end, strand) from a reference annotation file
(multi-FASTA format with GenBank-style location tags in headers).
Maps these coordinates to the corresponding positions in a gapped reference sequence
within the provided whole genome alignment (FASTA format).
Extracts the alignment columns corresponding to each gene for all sequences in the MSA.
Writes a separate FASTA alignment file for each successfully extracted gene.

Usage:
  python extract_genes_aln.py \\
      --annotations <reference_genes.fasta> \\
      --alignment <aligned_genomes.fasta> \\
      --ref_id <ID_of_reference_in_alignment_fasta> \\
      --output_dir <output_directory>
"""

import argparse
import os
import re
import sys
from pathlib import Path
from typing import List, Dict, Tuple, Optional # For type hinting

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# --- Annotation Parsing ---

def parse_genbank_location(location_str: str) -> Tuple[Optional[int], Optional[int], str]:
    """
    Parses GenBank feature location strings.

    Handles formats like 'start..end', 'complement(start..end)',
    and basic boundary markers ('<', '>'). It extracts the outermost
    coordinates from 'join' or 'order' statements.

    Args:
        location_str: The location string (e.g., "100..250", "complement(500..300)").

    Returns:
        A tuple (start, end, strand). Start and end are 1-based inclusive
        coordinates. Returns (None, None, strand) if parsing fails.
        Ensures start <= end.
    """
    strand = '+'
    original_input = location_str # Keep for potential error messages

    # Handle complement
    match_complement = re.match(r'complement\((.*)\)', location_str)
    if match_complement:
        strand = '-'
        location_str = match_complement.group(1)

    # Basic handling for join/order - extracts outermost range
    # WARNING: This simplification might be inaccurate for complex features (e.g., introns)
    location_str = location_str.replace('join(', '').replace('order(', '').replace(')', '')
    # Find all potential start..end pairs
    coords = re.findall(r'<?(\d+)\.\.>?(\d+)', location_str)

    if not coords:
        # Try single point location if no range found
        single_point = re.search(r'<?(\d+)>?', location_str)
        if single_point:
             pos = int(single_point.group(1))
             return pos, pos, strand # Treat single point as start=end
        else:
            print(f"Warning: Could not parse coordinates from location string: '{original_input}'", file=sys.stderr)
            return None, None, strand

    # Extract min start and max end from all found pairs
    try:
        all_pos = []
        for start_str, end_str in coords:
            all_pos.append(int(start_str))
            all_pos.append(int(end_str))

        start = min(all_pos)
        end = max(all_pos)

        # Ensure start <= end, even if original notation was end..start (common for complement)
        # We always store start <= end and use the strand info separately.
        if start > end:
            start, end = end, start # Swap

        return start, end, strand
    except ValueError:
        print(f"Warning: Error converting coordinates to numbers in location: '{original_input}'", file=sys.stderr)
        return None, None, strand


def parse_annotation_fasta(annotation_path: Path) -> List[Dict[str, any]]:
    """
    Reads a multi-FASTA reference gene file and extracts gene information.

    Parses headers like '>lcl|ID [gene=NAME] ... [location=LOC]' or using '[locus_tag=...]'.

    Args:
        annotation_path: Path to the annotation FASTA file.

    Returns:
        A list of dictionaries, each containing:
            'GeneName' (str),
            'Start' (int, 1-based ungapped),
            'End' (int, 1-based ungapped),
            'Strand' (str, '+' or '-'),
            'OriginalLocationStr' (str).
        Returns an empty list if the file is not found or no valid annotations are parsed.
    """
    gene_annotations = []
    required_fields_found = 0

    if not annotation_path.is_file():
         print(f"Error: Annotation file not found: {annotation_path}", file=sys.stderr)
         return []

    print(f"Parsing annotations from: {annotation_path}...")
    try:
        for record in SeqIO.parse(str(annotation_path), "fasta"):
            header = record.description
            gene_name: Optional[str] = None
            location_str: Optional[str] = None

            # Extract gene name ([gene=...] or [locus_tag=...])
            gene_match = re.search(r'\[(?:gene|locus_tag)=([^\]]+)\]', header)
            if gene_match:
                gene_name = gene_match.group(1)

            # Extract location string ([location=...])
            location_match = re.search(r'\[location=([^\]]+)\]', header)
            if location_match:
                location_str = location_match.group(1)

            if gene_name and location_str:
                start, end, strand = parse_genbank_location(location_str)
                if start is not None and end is not None:
                    gene_annotations.append({
                        'GeneName': gene_name,
                        'Start': start, # 1-based ungapped start
                        'End': end,     # 1-based ungapped end
                        'Strand': strand,
                        'OriginalLocationStr': location_str
                    })
                    required_fields_found += 1
                else:
                    print(f"Warning: Skipping record '{record.id}': Could not parse location '{location_str}'.", file=sys.stderr)
            # Optionally warn if gene/location missing but not skipping the record explicitly here
            # else:
            #     print(f"Warning: Gene name or location tag missing for record '{record.id}'. Cannot use for extraction.", file=sys.stderr)

    except Exception as e:
        print(f"Error parsing annotation file {annotation_path}: {e}", file=sys.stderr)
        return [] # Return empty list on critical error

    if required_fields_found == 0:
        print(f"Warning: No annotations with both a gene name/locus tag and a parsable location found in {annotation_path}.", file=sys.stderr)
    else:
        print(f"Found {required_fields_found} potential gene annotations.")

    return gene_annotations


# --- Coordinate Mapping ---

def map_coordinates_to_alignment(genes_info_list: List[Dict[str, any]],
                                 ref_aligned_record: SeqRecord
                                 ) -> List[Dict[str, any]]:
    """
    Maps ungapped 1-based gene coordinates to 0-based indices in the aligned reference sequence.

    Args:
        genes_info_list: List of gene info dicts from parse_annotation_fasta.
        ref_aligned_record: The reference SeqRecord *from the alignment*.

    Returns:
        Updated list of gene info dictionaries including:
            'Aln_Start_0based' (int): 0-based start index in the alignment.
            'Aln_End_Exclusive' (int): 0-based exclusive end index for slicing.
        Skips genes whose coordinates fall outside the reference sequence length
        or cannot be mapped.
    """
    if not ref_aligned_record:
        print("Error: Reference sequence record is missing, cannot map coordinates.", file=sys.stderr)
        return []

    ref_aligned_seq = str(ref_aligned_record.seq)
    mapped_genes_info = []

    # --- Pre-calculate mapping: ungapped_pos (1-based) -> aligned_index (0-based) ---
    ungapped_to_aligned_map: Dict[int, int] = {}
    ungapped_pos = 0
    for i, char in enumerate(ref_aligned_seq):
        if char != '-':
            ungapped_pos += 1
            ungapped_to_aligned_map[ungapped_pos] = i # Map 1-based ungapped to 0-based aligned index
    max_ref_ungapped_len = ungapped_pos
    # ---

    print(f"\nMapping coordinates relative to aligned reference '{ref_aligned_record.id}' "
          f"(ungapped length: {max_ref_ungapped_len})...")

    skipped_outside = 0
    skipped_mapping = 0
    for gene_info in genes_info_list:
        start_orig = gene_info['Start'] # 1-based ungapped
        end_orig = gene_info['End']     # 1-based ungapped

        # Validate original coordinates against the ungapped length of the reference
        if not (0 < start_orig <= max_ref_ungapped_len and 0 < end_orig <= max_ref_ungapped_len):
            # print(f"  Warning: Original coordinates {start_orig}..{end_orig} for gene '{gene_info['GeneName']}' "
            #       f"fall outside the reference's ungapped length ({max_ref_ungapped_len}). Skipping gene.", file=sys.stderr)
            skipped_outside += 1
            continue

        # Map ungapped start/end to aligned 0-based indices using the pre-calculated map
        try:
            aln_start_idx_0based = ungapped_to_aligned_map[start_orig]
            aln_end_idx_0based = ungapped_to_aligned_map[end_orig]

            # Store 0-based start index and exclusive end index (for Python slicing)
            # Ensure start <= end index even after mapping (should hold if start_orig <= end_orig)
            if aln_start_idx_0based <= aln_end_idx_0based:
                 gene_info['Aln_Start_0based'] = aln_start_idx_0based
                 gene_info['Aln_End_Exclusive'] = aln_end_idx_0based + 1 # Make end exclusive
                 mapped_genes_info.append(gene_info)
            else:
                 # This case should ideally not happen if start_orig <= end_orig
                 print(f"  Warning: Mapped alignment indices are reversed for gene '{gene_info['GeneName']}' "
                       f"({aln_start_idx_0based}..{aln_end_idx_0based}). Skipping gene.", file=sys.stderr)
                 skipped_mapping += 1
                 continue

        except KeyError as e:
             # This occurs if start_orig or end_orig doesn't correspond to an ungapped base
             # Should be caught by the length check above, but kept as safeguard.
             print(f"  Warning: Could not map coordinate {e} for gene '{gene_info['GeneName']}'. "
                   f"Original coords {start_orig}..{end_orig}. Skipping gene.", file=sys.stderr)
             skipped_mapping += 1
             continue

    print(f"Coordinate mapping done:")
    print(f"  Successfully mapped: {len(mapped_genes_info)} genes.")
    if skipped_outside > 0:
         print(f"  Skipped (outside ref len): {skipped_outside} genes.")
    if skipped_mapping > 0:
        print(f"  Skipped (mapping issue): {skipped_mapping} genes.")

    return mapped_genes_info

# --- File Handling ---

def sanitize_filename(name: str) -> str:
    """Removes or replaces characters potentially problematic for filenames."""
    # Remove potentially confusing characters like brackets, parentheses, slashes, colons
    name = re.sub(r'[\[\]()/:\\\'"]', '', name)
    # Replace whitespace with underscore
    name = name.replace(' ', '_')
    # Keep only alphanumeric, underscore, hyphen, dot
    # Allow dot for potential extensions if needed, but usually not for gene names
    name = re.sub(r'[^\w\-.]', '', name)
    # Avoid names starting with '.' or '-' which can be hidden or misinterpreted
    if name.startswith('.') or name.startswith('-'):
        name = '_' + name
    # Avoid empty filename
    if not name:
        name = "unnamed_gene"
    return name


# --- Main Execution ---
def main() -> None:
    """Parses arguments, orchestrates parsing, mapping, extraction, and writing."""
    parser = argparse.ArgumentParser(
        description="Extracts gene alignments from a whole genome alignment based on reference annotations.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-a", "--annotations",
        type=Path, required=True,
        help="Path to the reference gene annotation file (multi-FASTA with [gene=/locus_tag=] and [location=] tags)."
    )
    parser.add_argument(
        "-g", "--alignment",
        type=Path, required=True,
        help="Path to the whole genome multiple sequence alignment file (FASTA format)."
    )
    parser.add_argument(
        "-r", "--ref_id",
        type=str, required=True,
        help="Sequence ID (must exactly match a FASTA header ID) of the reference genome within the --alignment file."
    )
    parser.add_argument(
        "-o", "--output_dir",
        type=Path, required=True,
        help="Path to the output directory where individual gene alignment files will be saved."
    )
    args = parser.parse_args()

    # --- Validate inputs ---
    if not args.alignment.is_file():
        print(f"Error: Genome alignment file not found: {args.alignment}", file=sys.stderr)
        sys.exit(1)
    # Annotation file existence checked within parse_annotation_fasta

    # --- Create output directory ---
    try:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        print(f"Output directory: {args.output_dir}")
    except OSError as e:
        print(f"Error creating output directory '{args.output_dir}': {e}", file=sys.stderr)
        sys.exit(1)

    # --- Parse Annotations ---
    gene_annotations = parse_annotation_fasta(args.annotations)
    if not gene_annotations: # Checks for empty list
        print("Exiting: No valid annotations could be parsed.", file=sys.stderr)
        sys.exit(1)

    # --- Read Genome Alignment and Find Reference ---
    print(f"\nReading genome alignment: {args.alignment}...")
    ref_aligned_record: Optional[SeqRecord] = None
    genome_alignment_records: List[SeqRecord] = []
    try:
        # Read all records into memory. Warning: Can be memory intensive for large alignments.
        # Consider SeqIO.index for very large files if memory becomes an issue.
        temp_records = list(SeqIO.parse(str(args.alignment), "fasta"))

        if not temp_records:
             print(f"Error: No sequences found in genome alignment file: {args.alignment}", file=sys.stderr)
             sys.exit(1)

        for record in temp_records:
            genome_alignment_records.append(record)
            # Use record.id for matching, which is typically the part before the first space
            if record.id == args.ref_id:
                ref_aligned_record = record

        if not ref_aligned_record:
             print(f"Error: Reference sequence ID '{args.ref_id}' not found in the alignment file headers.", file=sys.stderr)
             print(f"       Available IDs: {[rec.id for rec in genome_alignment_records[:5]]}...", file=sys.stderr) # Show first few IDs
             sys.exit(1)

        print(f"Read {len(genome_alignment_records)} genome sequences. Found reference '{args.ref_id}'.")

    except Exception as e:
        print(f"Error reading or parsing genome alignment file {args.alignment}: {e}", file=sys.stderr)
        sys.exit(1)

    # --- Map Coordinates ---
    # Pass only the found reference record for mapping
    aligned_genes_info = map_coordinates_to_alignment(gene_annotations, ref_aligned_record)
    if not aligned_genes_info:
        print("Exiting: Could not map coordinates for any gene.", file=sys.stderr)
        sys.exit(1)


    # --- Extract and Write Gene Alignments ---
    print("\nExtracting and writing gene alignments...")
    genes_written = 0
    genes_failed = 0
    genes_skipped_extraction = 0

    for gene_info in aligned_genes_info:
        gene_name = gene_info['GeneName']
        safe_gene_name = sanitize_filename(gene_name)
        output_filename = args.output_dir / f"gene_{safe_gene_name}.fasta"
        # Short status message per gene
        # print(f"  Processing {gene_name} -> {output_filename} ...")

        gene_specific_records: List[SeqRecord] = []
        start_aln = gene_info['Aln_Start_0based']      # 0-based start index
        end_aln = gene_info['Aln_End_Exclusive']    # Exclusive end index for slicing
        strand = gene_info['Strand']

        # Double check mapped coords validity (should be okay if mapping worked)
        if start_aln < 0 or end_aln <= start_aln:
             print(f"    Warning: Invalid mapped coordinates ({start_aln}..{end_aln}) for gene {gene_name}. Skipping.", file=sys.stderr)
             genes_failed += 1
             continue

        sequences_extracted_count = 0
        extraction_skipped_for_genome = False
        for genome_record in genome_alignment_records:
            genome_seq = genome_record.seq
            genome_id = genome_record.id # Use the simple ID
            genome_len = len(genome_seq)

            # Ensure alignment coordinates are within the bounds of this specific genome sequence
            if end_aln > genome_len:
                # This might happen if alignments have different lengths, though less common in WGAs
                # print(f"    Warning: Mapped coordinates for gene {gene_name} ({start_aln}..{end_aln}) "
                #       f"exceed length of genome {genome_id} ({genome_len}). Skipping this genome for this gene.", file=sys.stderr)
                extraction_skipped_for_genome = True # Mark that at least one genome was skipped
                continue # Skip this genome for this gene

            # Extract subsequence using mapped alignment coordinates
            sub_sequence: Seq = genome_seq[start_aln:end_aln]

            # Handle reverse strand if necessary AFTER extraction
            if strand == '-':
                # Use Biopython's reverse complement which handles sequence types
                try:
                    sub_sequence = sub_sequence.reverse_complement()
                except Exception as rc_err:
                     print(f"    Warning: Could not reverse complement sequence for {genome_id} gene {gene_name}: {rc_err}. Using forward strand sequence.", file=sys.stderr)
                     # Decide whether to skip or use forward strand - using forward strand here

            # Create a new record for the extracted gene segment
            extracted_record = SeqRecord(
                sub_sequence,
                id=genome_id, # Keep the original genome ID
                description=f"gene={gene_name} | source_location={gene_info['OriginalLocationStr']} | source_strand={strand}"
            )
            gene_specific_records.append(extracted_record)
            sequences_extracted_count += 1

        # Write the alignment for this gene if any sequences were successfully extracted
        if gene_specific_records:
            try:
                with open(output_filename, "w") as outfile:
                    SeqIO.write(gene_specific_records, outfile, "fasta")
                genes_written += 1
            except IOError as e:
                print(f"    Error writing output file {output_filename}: {e}", file=sys.stderr)
                genes_failed += 1
            except Exception as e:
                 print(f"    Unexpected error writing file {output_filename}: {e}", file=sys.stderr)
                 genes_failed += 1
        else:
             # This gene had mapped coordinates, but no sequence could be extracted
             # (e.g., coordinates exceeded length for *all* genomes)
             genes_skipped_extraction += 1
             # print(f"    Warning: No sequences could be extracted for gene {gene_name}. No output file generated.", file=sys.stderr)


    print(f"\n--- Processing Summary ---")
    print(f"  Total annotations parsed:         {len(gene_annotations)}")
    print(f"  Annotations successfully mapped:  {len(aligned_genes_info)}")
    print(f"  Gene alignments written:          {genes_written}")
    if genes_skipped_extraction > 0:
        print(f"  Genes skipped (no seq extracted): {genes_skipped_extraction}")
    if genes_failed > 0:
        print(f"  Genes failed (error):             {genes_failed}")

    # Final check
    if genes_written == 0 and len(aligned_genes_info) > 0 :
        print("\nWarning: No gene alignment files were written, despite successful mapping. Check warnings above regarding sequence extraction issues (e.g., coordinate bounds).", file=sys.stderr)


if __name__ == "__main__":
    main()