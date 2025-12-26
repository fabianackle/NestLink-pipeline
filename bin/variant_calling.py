#!/usr/bin/env python3
import argparse
import re
import dnaio
import polars as pl
from Bio import SeqIO
from Bio.Seq import Seq


def extract_sequences(file_path, barcode, orf, translate_barcode):
    """
    Parse a FASTA file to extract barcode and ORF sequences.

    Args:
        file_path (str): Path to the FASTA file.
        barcode (tuple): A (start_pattern, end_pattern) tuple for the barcode.
        orf (tuple): A (start_pattern, end_pattern) tuple for the ORF.
        translate_barcode (bool): Translate DNA barcodes, used for flycodes.

    Returns:
        dict: Keys are cluster IDs; values are tuples of extracted sequences (barcode_seq, orf_seq) in-frame only.
    """
    def is_in_frame(start_match, end_match):
        return (end_match.start() - start_match.end()) % 3 == 0

    def extract_subsequence(sequence, start_match, end_match):
        return sequence[start_match.end():end_match.start()]

    sequences = {}

    bc_start = re.compile(barcode[0], re.IGNORECASE)
    bc_end = re.compile(barcode[1], re.IGNORECASE)

    orf_start = re.compile(orf[0], re.IGNORECASE)
    orf_end = re.compile(orf[1], re.IGNORECASE)

    with dnaio.open(file_path, mode="r") as reader:
        for record in reader:
            barcode_id = record.name
            sequence = record.sequence

            # Search for barcodes
            bc_start_match = bc_start.search(sequence)
            bc_end_match = bc_end.search(sequence)

            # Search for orf
            orf_start_match = orf_start.search(sequence)
            orf_end_match = orf_end.search(sequence)

            # Extract sequences if matches are found and in-frame
            # When barcodes are translated, they have to be in-frame.
            if bc_start_match and bc_end_match and (not translate_barcode or is_in_frame(bc_start_match, bc_end_match)):
                barcode_sequence = extract_subsequence(sequence, bc_start_match, bc_end_match)

                if orf_start_match and orf_end_match and is_in_frame(orf_start_match, orf_end_match):
                    orf_sequence = extract_subsequence(sequence, orf_start_match, orf_end_match)
                    sequences[barcode_id] = (barcode_sequence, orf_sequence)

    return sequences


def translate(sequence):
    """
    Translate a nucleotide sequence into an amino acid sequence.

    Args:
        sequence (str): The nucleotide sequence.

    Returns:
        str: The corresponding amino acid sequence.
    """
    return str(Seq(sequence).translate())


def compare_aa_sequences(sequence, reference):
    """
     Compare two amino acid sequences and identify positions where they differ.

    Args:
        sequence (str): The amino acid sequence to compare.
        reference (str): The reference amino acid sequence.

    Returns:
        list: A list of amino acid changes, as tuples (position, refAA, newAA).
    """
    changes = [
        (position, ref_base, seq_base)
        for position, (ref_base, seq_base) in enumerate(zip(reference, sequence))
        if ref_base != seq_base
    ]

    return changes


def get_variants(sequences, reference, translate_barcode):
    """
    Determine amino acid changes for each consensus sequence relative to a the reference.

    Args:
        sequences (dict): Dictionary mapping cluster IDs to extracted (barcode, ORF) sequences.
        reference (dict): Dictionary containing exactly one reference entry of the same format.
        translate_barcode (bool): Translate DNA barcodes, used for flycodes.

    Returns:
       pl.DataFrame: Contains barcode_id, barcode, amino acid change position, reference amino acid, variant amino acid and variant type.
    """
    data = []

    # Getting reference sequences on a protein level
    if len(reference) != 1:
        raise ValueError(
            "The reference file does not contain exactly one record."
        )
    _, reference_sequence = next(iter(reference.items()))
    reference_orf_aa = translate(reference_sequence[1])

    # Looping through all sequences
    for barcode_id, (barcode_nt, orf_nt) in sequences.items():
        if (translate_barcode):
            barcode = translate(barcode_nt)
        else:
            barcode = barcode_nt

        orf_aa = translate(orf_nt)

        if (len(orf_aa) != len(reference_orf_aa)):
            data.append({
                "barcode_id": barcode_id,
                "barcode": barcode,
                "variant_type": "indel",
                "position": None,
                "reference_aa": None,
                "variant_aa": None
            })
            continue

        if (orf_aa == reference_orf_aa):
            data.append({
                "barcode_id": barcode_id,
                "barcode": barcode,
                "variant_type": "wt",
                "position": None,
                "reference_aa": None,
                "variant_aa": None
            })
            continue

        if (orf_aa != reference_orf_aa):
            orf_aa_changes = compare_aa_sequences(orf_aa, reference_orf_aa)
            data.extend({
                "barcode_id": barcode_id,
                "barcode": barcode,
                "variant_type": "change",
                "position": pos,
                "reference_aa": ref_aa,
                "variant_aa": var_aa
            } for pos, ref_aa, var_aa in orf_aa_changes)

    return pl.DataFrame(data)


def main(args):
    """
    Coordinate the extraction, translation, variant calling, and database writing of barcodes.

    Args:
        args (Namespace): Parsed command-line arguments.
    """
    assembly_path = args.assembly_path
    reference_path = args.reference_path
    sample_id = args.sample_id
    barcode_pattern = args.barcode_pattern
    orf_pattern = args.orf_pattern
    translate_barcode = args.translate_barcode

    sequences = extract_sequences(assembly_path, barcode_pattern, orf_pattern, translate_barcode)
    reference = extract_sequences(reference_path, barcode_pattern, orf_pattern, translate_barcode)

    variants = get_variants(sequences, reference, translate_barcode)
    variants.write_csv(f"{sample_id}_variants.csv")

    # writing the barcode map for enrich2
    barcode_map = [
        {"barcode": barcode, "orf": orf}
        for _, (barcode, orf) in sequences.items()
    ]
    pl.DataFrame(barcode_map).write_csv(f"{sample_id}_barcodemap.txt", include_header=False, separator="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--assembly_path", type=str, required=True)
    parser.add_argument("--reference_path", type=str, required=True)
    parser.add_argument("--sample_id", type=str, required=True)
    parser.add_argument("--barcode_pattern", nargs=2, type=str, required=True)
    parser.add_argument("--orf_pattern", nargs=2, type=str, required=True)
    parser.add_argument("--translate_barcode", action='store_true')

    args = parser.parse_args()
    main(args)
