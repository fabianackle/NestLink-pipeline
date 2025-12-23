#!/usr/bin/env python3

import argparse
import csv
import re

import dnaio


def read_references(file_path, orf_5p, orf_3p):
    """
    Read reference sequences and extract ORF start/end positions (1-based).
    """
    records = {}

    orf_start = re.compile(orf_5p, re.IGNORECASE)
    orf_end = re.compile(orf_3p, re.IGNORECASE)

    with dnaio.open(file_path) as reader:
        for record in reader:
            start_match = orf_start.search(record.sequence)
            end_match = orf_end.search(record.sequence)

            if not start_match or not end_match:
                continue

            records[record.name] = (
                start_match.end() + 1,
                end_match.start(),
            )

    return records


def write_gff(references, output_path):
    """
    Write GFF3 records from ORF coordinates.
    """
    with open(output_path, "w", newline="") as file:
        file.write("##gff-version 3\n")

        writer = csv.writer(file, delimiter="\t")
        for record, (start, end) in references.items():
            writer.writerow([
                record, ".", "gene", start, end, ".", "+", ".",
                f"ID=gene:{record}_g;biotype=protein_coding",
            ])
            writer.writerow([
                record, ".", "transcript", start, end, ".", "+", ".",
                f"ID=transcript:{record}_t;Parent=gene:{record}_g;biotype=protein_coding;",
            ])
            writer.writerow([
                record, ".", "CDS", start, end, ".", "+", "0",
                f"Parent=transcript:{record}_t;",
            ])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--references", required=True)
    parser.add_argument("--orf_5p", required=True)
    parser.add_argument("--orf_3p", required=True)
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    references = read_references(
        args.references,
        args.orf_5p,
        args.orf_3p,
    )

    write_gff(references, args.output)


if __name__ == "__main__":
    main()
