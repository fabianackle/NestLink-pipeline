#!/usr/bin/env python3
import argparse
import subprocess
import dnaio
import duckdb
import polars as pl


def barcodes_to_dataframe(file_path):
    """
    Read a FASTA of extracted barcodes into a polars DataFrame with read IDs and their corresponding barcodes.

    Args:
        file_path (str): Path to the FASTA file with extracted barcodes.

    Returns:
        pl.DataFrame: A DataFrame with 'read_id' and 'barcode' columns.
    """
    data = []

    with dnaio.open(file_path) as reader:
        for record in reader:
            row = {"read_id": record.name, "barcode": record.sequence}
            data.append(row)

    if not data:
        raise ValueError(
            f"No barcodes found in file '{file_path}'. "
            "Please check the barcode and or sequence extraction processes."
        )

    return pl.DataFrame(data)


def map_barcodes_to_hq_barcodes(file_path, data=[]):
    """
    Parse a SAM file of barcodes aligned to high-quality barcodes into a polars DataFrame mapping read IDs to barcode IDs.

    Args:
        file_path (str): Path to the SAM file with aligned barcodes.

    Returns:
        pl.DataFrame: A DataFrame with 'read_id', 'barcode_id', 'cigar' and 'edit_distance' columns.
    """
    with open(file_path, "r") as f:
        for line in f:
            # Skip header lines
            if line.startswith("@"):
                continue

            # Split the SAM fields
            fields = line.strip().split("\t")

            # Extract necessary fields
            read_id = fields[0]     # read_id (Query template NAME - QNAME)
            barcode_id = fields[2]  # Reference barcode_id (Reference sequence NAME - RNAME)
            flag = int(fields[1])   # SAM flag
            cigar = fields[5]       # CIGAR string

            # Check if the alignment is valid (flag 4 means unaligned)
            if barcode_id != "*" and not (flag & 4):
                edit_distance = fields[12].split(":")[2]
                row = {"read_id": read_id, "barcode_id": barcode_id, "cigar": cigar, "edit_distance": edit_distance}
                data.append(row)

    if not data:
        raise ValueError(
            f"No aligned barcodes found in file '{file_path}'."
        )

    return pl.DataFrame(data, schema_overrides={"edit_distance": pl.UInt8})


def write_references(sample_id, hq_barcodes_df, reference_seq):
    """
    Write reference FASTA file containing all high-quality barcodes.

    Args:
        hq_barcodes_df (pl.DataFrame): DataFrame with barcode IDs and their barcode sequences.
        reference_seq (str): Path to the reference sequence FASTA file.
    """
    # Read the reference sequence.
    with dnaio.open(reference_seq, mode="r") as reader:
        reference_record = None
        for record in reader:
            reference_record = record
            break  # stop after the first record

        if reference_record is None:
            raise ValueError(
                f"No record found in the reference sequence file '{reference_seq}'. "
                "Please check the file."
            )

        reference_sequence = reference_record.sequence

    records = []  # for storing all barcode records.
    for barcode_id in hq_barcodes_df.get_column("barcode_id").unique():
        record = dnaio.SequenceRecord(barcode_id, reference_sequence)
        records.append(record)

    # Writing a reference file for polishing.
    file_name = f"{sample_id}_references.fasta"
    with dnaio.open(file_name, mode="w") as writer:
        for record in records:
            writer.write(record)


def main(sample_id:str, reference_seq:str, barcodes:str, barcode_min_coverage:int, barcode_regex:str, threads:int):
    """
    Main function to determine high-quality barcodes and map reads to these high-quality barcodes.

    Args:
        sample_id (str): The sample id.
        barcodes (str): Path to the FASTA file with extracted barcodes.
        reference_seq (str): Path to the reference sequence FASTA file.
        barcode_regex (srt): Regex that matches the used barcode.
        threads (int): Threads for alignment.
    """
    # Reading in the barcodes.
    barcodes_df = barcodes_to_dataframe(barcodes)

    # Add a new column with True/False indicating a valid barcode.
    barcodes_df = barcodes_df.with_columns(
        pl.col("barcode").str.contains(barcode_regex).alias("is_valid_barcode")
    )
    barcodes_df.write_csv(f"{sample_id}_reads.csv")

    # Getting high-quality barcodes, they have to be valid and be seen more than 10 times.
    hq_barcodes_df = duckdb.sql(
        f"""
        SELECT 
            gen_random_uuid() AS barcode_id, 
            barcode
        FROM barcodes_df
        WHERE is_valid_barcode
        GROUP BY barcode
        HAVING COUNT(barcode) >= {barcode_min_coverage};
        """
    ).pl()
    hq_barcodes_df.write_csv(f"{sample_id}_hq_barcodes.csv")

    high_quality_barcodes = {sequence:barcode_id for barcode_id, sequence in hq_barcodes_df.rows()}

    # Writing the high-quality barcodes into hq_barcodes.fasta to be used as a reference for alignment.
    with dnaio.open("hq_barcodes.fasta", mode="w") as writer:
        for sequence, barcode_id in high_quality_barcodes.items():
            writer.write(dnaio.SequenceRecord(barcode_id, sequence))

    # Writing the low-quality barcodes into lq_barcodes.fasta to be aligned to the high-quality barcodes.
    mapped_hq_barcodes = []
    with dnaio.open(f"lq_barcodes.fasta", mode="w") as writer:
        for read_id, barcode, _ in barcodes_df.rows():
            if barcode in high_quality_barcodes:
                mapped_hq_barcodes.append({"read_id": read_id, "barcode_id": high_quality_barcodes[barcode], "cigar": f'{len(barcode)}M', "edit_distance": 0})
            else:
                writer.write(dnaio.SequenceRecord(read_id, barcode))

    # Indexing the barcodes.
    subprocess.run(["bwa", "index", "hq_barcodes.fasta"], check=True)

    # Aligning low-quality barcodes to the high-quality barcodes.
    with open("lq_to_hq.sai", "w") as sai_file:
        subprocess.run(["bwa", "aln", "-t", str(threads), "-N", "-n 2", "hq_barcodes.fasta", "lq_barcodes.fasta"], stdout=sai_file, check=True)

    # Generating alignments in the SAM format,
    with open("lq_to_hq.sam", "w") as sam_file:
        subprocess.run(["bwa", "samse", "hq_barcodes.fasta", "lq_to_hq.sai", "lq_barcodes.fasta"], stdout=sam_file, check=True)

    # Mapping barcodes to high quality barcodes.
    mapped_barcodes_df = map_barcodes_to_hq_barcodes("lq_to_hq.sam", data=mapped_hq_barcodes)
    mapped_barcodes_df.write_csv(f"{sample_id}_mapped_reads.csv")

    # Filtering mapped_barcodes_df to only contain up to 100 reads with low edit distance barcodes per high quality barcode.
    mapped_barcodes_df = duckdb.sql(
        f"""
        WITH ranked_reads AS (
            SELECT 
                read_id,
                barcode_id,
                cigar,
                edit_distance,
                row_number() OVER (PARTITION BY barcode_id ORDER BY edit_distance ASC) AS row_num,
                COUNT(*) OVER (PARTITION BY barcode_id) AS barcode_coverage
            FROM mapped_barcodes_df
        )
        SELECT 
            read_id,
            barcode_id,
            cigar,
            edit_distance
        FROM ranked_reads
        WHERE
            row_num <= 100 AND
            barcode_coverage >= {barcode_min_coverage};
        """
    ).pl()
    mapped_barcodes_df.write_csv(f"{sample_id}_mapped_reads_filtered.csv")

    # Writing references to disk.
    write_references(sample_id, mapped_barcodes_df, reference_seq)

if __name__ == "__main__":
    # Reading arguments and calling the main function.
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_id", type=str,
        help="The sample id.")
    parser.add_argument("--reference_seq", type=str,
        help="Path to the reference sequence FASTA file.")
    
    parser.add_argument("--barcodes", type=str,
        help="Path to the FASTA file with extracted barcodes.")
    parser.add_argument("--barcode_min_coverage", type=int, default=10,
        help="The minimal amount a barcode has to be seen to be considered a high-quality barcode.")
    parser.add_argument("--barcode_regex", type=str,
        help="Regex that matches the used barcode.")
    
    parser.add_argument("--threads", type=int,
        help="threads (int): Threads for alignment.")

    args = parser.parse_args()
    main(args.sample_id, args.reference_seq, args.barcodes, args.barcode_min_coverage, args.barcode_regex, args.threads)
