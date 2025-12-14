process GROUP_BY_BARCODES {
    conda "bioconda::bwa=0.7.19 bioconda::samtools=1.22.1 bioconda::dnaio=1.2.3 conda-forge::polars=1.35.1 conda-forge::pyarrow=22.0.0 conda-forge::python-duckdb=1.4.1"
    container "${ workflow.containerEngine == 'apptainer' ?
        'oras://community.wave.seqera.io/library/bwa_dnaio_samtools_polars_pruned:a0bc9e59a0317e84' :
        'community.wave.seqera.io/library/bwa_dnaio_samtools_polars_pruned:20f93e2b62ebb5b6' }"

    tag "${sample_id}"

    input:
    tuple val(sample_id), path(barcodes), path(reference)

    output:
    tuple val(sample_id), path("references.fasta"), emit: references
    tuple val(sample_id), path("${sample_id}_mapped_reads_filtered.csv"), emit: barcode_map
    tuple val(sample_id), path("${sample_id}_reads.csv"), path("${sample_id}_clusters.csv"), path("${sample_id}_mapped_reads.csv"), path("${sample_id}_mapped_reads_filtered.csv"), emit: csv

    script:
    """
    group_by_barcodes.py \
        --sample_id ${sample_id} \
        --reference_seq ${reference} \
        --barcodes ${barcodes} \
        --barcode_regex "${params.barcode_regex}" \
        --barcode_min_coverage ${params.barcode_min_coverage} \
        --threads ${task.cpus}
    """
}
