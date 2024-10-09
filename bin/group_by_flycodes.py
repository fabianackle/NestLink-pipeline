#!/usr/bin/env python3
import argparse
import os
import re
import uuid
import dnaio
import pandas as pd
from collections import Counter
from collections import defaultdict


def parse_centroids_file(file_path, min_cluster_sizes):
    """
    Reads the flycode_centroids.fasta file and returns a list of clusterids with clustersizes equal or larger than min_cluster_sizes.
    """
    centroids = []
    with dnaio.open(file_path) as reader:
        for record in reader:
            readid, clusterid, size = parse_header(record.name)
            if size >= min_cluster_sizes:
                centroids.append(clusterid)
    return centroids


def parse_header(fasta_header):
    """
    Parses the header of a flycode_centroids file and returns a tuple containg the readid, direction of where the flycode was found, the clusterid and the size of the cluster.
    """
    parts = fasta_header.split(';')
    readid = parts[0]
    clusterid = int(parts[1].split('=')[1])
    size = int(parts[2].split('=')[1])
    return (readid, clusterid, size)


def is_valid_flycode(seq):
    pattern = r"^GGTAGT(GCA|GTT|GAT|CCA|GAA|ACT|GGT|TCT|TAC|CTG|TGG|CAG|TTC|AAC){6,8}(TGGCGG|TGGCTGCGG|TGGCAGTCTCGG|TGGCAGGAAGGAGGTCGG)$"
    return re.match(pattern, seq) is not None


def read_clusters(folder_path, cluster_ids):
    """
    Read clusters from files and return two dictionaries. One dictionary maps cluster IDs to the reference flycode and the most common flycode as a tuple. The other dictionary maps read IDs to their cluster IDs.
    Only clusters with a most abundant flycode that appears at least 10 times are considered.
    Args:
    folder_path (str): The path to the folder containing the cluster files.
    cluster_ids (list): A list of cluster IDs to be processed.
    Returns:
    tuple: A tuple containing two dictionaries. The first dictionary maps maps read IDs to their corresponding cluster IDs. The second dictionary maps cluster IDs to a tuple of the read ID of a the reference flycode and the flycode sequence.
    """
    reads = {}
    references = {}
    # Iterate over each cluster_id
    for cluster_id in cluster_ids:
        file_name = 'cluster.fasta' + str(cluster_id)
        file_path = os.path.join(folder_path, file_name)

        # Open the respective cluster file.
        with dnaio.open(file_path) as reader:
            flycodes_list = []  # all flycodes of the cluster
            reads_temp = {}  # Temporary dictionary to map read IDs to their cluster_id for the current cluster file.

            for record in reader:
                readid, _ = record.name.split(';')
                flycode = record.sequence
                reads_temp[readid] = (flycode, cluster_id)
                # Collect all flycodes in the current cluster to count occurrences later.
                flycodes_list.append(flycode)

            # Count occurrences of each flycode in the cluster.
            flycodes_counter = Counter(flycodes_list)
            # Identify the most frequently occurring flycode in the cluster.
            most_common_flycode, count = flycodes_counter.most_common(1)[0]
            # Check if the most common flycode is valid and occurs at least 10 times.
            if is_valid_flycode(most_common_flycode) and count >= 10:
                cluster_uuid = str(uuid.uuid4())
                # Update the reads dictionary with read IDs from the current cluster.
                reads = reads | reads_temp
                # Record the reference and most sequence of the most common flycode for the current cluster.
                references[cluster_id] = (cluster_uuid, most_common_flycode)
    return references, reads


def bin_reads_by_flycodes(reads_path, reads_dict):
    """
    Returns a dictionary with clusterids as key and a list of corresponding reads (records) as value.
    """
    reads_binned = defaultdict(list)
    with dnaio.open(reads_path) as reader:
        for record in reader:
            readid = record.name.split('	')[0]
            if readid in reads_dict:
                _, cluster_id = reads_dict[readid]
                reads_binned[cluster_id].append(record)
    return reads_binned


def write_binned_reads(folder_path, binned_reads):
    '''
    Writes new fastq.gz files containing reads binned by flycodes for alignment. The direction of the reads is all fwd.
    '''
    for cluster_id, reads in binned_reads.items():
        file_name = 'reads_fc' + str(cluster_id) + '.fastq.gz'
        file_path = os.path.join(folder_path, file_name)
        with dnaio.open(file_path, mode="w") as writer:
            for record in reads:
                writer.write(record)


def write_reference_files(folder_path, reference, references):
    '''
    Create and write reference files for clusters. This function writes two types of files: one for each cluster containing for alignment and a common reference file containing all reference reads for medaka stitching.

    Args:
    folder_path (str): The directory path where the reference files will be saved.
    reference (str): The path to the initial reference file.
    references (dict): A dictionary mapping cluster IDs to tuples of read IDs and corresponding flycodes.
    '''

    # Open and read the reference sequence.
    with dnaio.open(reference, mode="r") as reader:
        for record in reader:
            reference = record.sequence
            break
    # Splitting the reference intwo two parts using the flycode as delimiter.
    reference = reference.split("GGTAGTNNNNNNNNNNNNNNNNNNNNNTGGcgg")

    records = []
    for cluster_id, (read_id, flycode) in references.items():

        # Create and write a fasta file for each cluster with the correct flycode inserted into the reference sequence.
        file_name = f'{folder_path}/reference_fc{cluster_id}.fasta'
        with dnaio.open(file_name, mode="w") as writer:
            reference_sequence = f'{reference[0]}{flycode}{reference[1]}'
            record = dnaio.SequenceRecord(read_id, reference_sequence)
            writer.write(record)
            # Add the newly created reference sequence for this cluster to a list for the common reference file.
            records.append(record)

    # Write all individual cluster references into a single common reference fasta file.
    file_name = f'{folder_path}/reference.fasta'
    with dnaio.open(file_name, mode="w") as writer:
        for record in records:
            writer.write(record)


def write_csv(references_dict, reads_dict):
    df_clusters = pd.DataFrame.from_dict(references_dict, orient='index', columns=['cluster_uuid', 'most_common_flycode'])
    df_clusters.reset_index(inplace=True)
    df_clusters.rename(columns={'index': 'cluster_id'}, inplace=True)
    df_clusters.to_csv('clusters.csv', index=False)

    df_reads = pd.DataFrame.from_dict(reads_dict, orient='index', columns=['flycode', 'cluster_id'])
    df_reads.reset_index(inplace=True)
    df_reads.rename(columns={'index': 'read_id'}, inplace=True)
    df_reads.to_csv('reads.csv', index=False)


def main(centroids, clusters, outdir):
    cluster_ids = parse_centroids_file(centroids, 10)
    references_dict, reads_dict = read_clusters(clusters, cluster_ids)
    write_csv(references_dict, reads_dict)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--centroids", type=str)
    parser.add_argument("--clusters", type=str)
    parser.add_argument("--outdir", type=str)

    args = parser.parse_args()
    main(args.centroids, args.clusters, args.outdir)
