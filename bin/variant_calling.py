#!/usr/bin/env python3
import argparse
import re
import dnaio
import polars as pl
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq


def extract_sequences(file_path, barcode, orf):
    """
    Parse a FASTA file to extract barcode and ORF sequences.

    Args:
        file_path (str): Path to the FASTA file.
        barcode (tuple): A (start_pattern, end_pattern) tuple for the barcode.
        orf (tuple): A (start_pattern, end_pattern) tuple for the ORF.

    Returns:
        dict: Keys are cluster IDs; values are tuples of extracted sequences (barcode_seq, orf_seq) in-frame only.
    """
    def is_in_frame(start_match, end_match):
        return (end_match.start() - start_match.end()) % 3 == 0

    def extract_subsequence(sequence, start_match, end_match):
        return sequence[start_match.end():end_match.start()]

    sequences = {}

    fc_start = re.compile(barcode[0], re.IGNORECASE)
    fc_end = re.compile(barcode[1], re.IGNORECASE)

    orf_start = re.compile(orf[0], re.IGNORECASE)
    orf_end = re.compile(orf[1], re.IGNORECASE)

    with dnaio.open(file_path, mode="r") as reader:
        for record in reader:
            cluster_id = record.name
            sequence = record.sequence

            # Search for barcodes
            fc_start_match = fc_start.search(sequence)
            fc_end_match = fc_end.search(sequence)

            # Search for orf
            orf_start_match = orf_start.search(sequence)
            orf_end_match = orf_end.search(sequence)

            # Extract sequences if matches are found and in-frame
            if fc_start_match and fc_end_match and is_in_frame(fc_start_match, fc_end_match):
                barcode_sequence = extract_subsequence(sequence, fc_start_match, fc_end_match)

                if orf_start_match and orf_end_match and is_in_frame(orf_start_match, orf_end_match):
                    orf_sequence = extract_subsequence(sequence, orf_start_match, orf_end_match)
                    sequences[cluster_id] = (barcode_sequence, orf_sequence)

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
        translate_barcode (bool): Translate DNA barcodes, used for barcodes.

    Returns:
       pl.DataFrame: Contains cluster_ID, barcode, amino acid change position, reference amino acid, variant amino acid and variant type.
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
    for cluster_id, (barcode_nt, orf_nt) in sequences.items():
        if (translate_barcode):
            barcode = translate(barcode_nt)
        else:
            barcode = barcode_nt

        orf_aa = translate(orf_nt)

        if (len(orf_aa) != len(reference_orf_aa)):
            data.append({"cluster_id": cluster_id, "barcode": barcode, "variant_type": "indel"})
            continue

        if (orf_aa == reference_orf_aa):
            data.append({"cluster_id": cluster_id, "barcode": barcode, "variant_type": "wt"})
            continue

        if (orf_aa != reference_orf_aa):
            orf_aa_changes = compare_aa_sequences(orf_aa, reference_orf_aa)
            data.extend({
                "cluster_id": cluster_id,
                "barcode": barcode,
                "position": pos,
                "reference_aa": ref_aa,
                "variant_aa": var_aa,
                "variant_type": "change"
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

    sequences = extract_sequences(assembly_path, barcode_pattern, orf_pattern)
    reference = extract_sequences(reference_path, barcode_pattern, orf_pattern)

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
