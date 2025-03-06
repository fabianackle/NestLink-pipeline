process GROUP_BY_FLYCODES {
    cpus 8
    memory '4 GB'
    time '60m'
    conda "bioconda::bwa=0.7.18 bioconda::samtools=1.21 bioconda::dnaio=1.2.3 conda-forge::polars=1.22.0 conda-forge::pyarrow=19.0.1 conda-forge::python-duckdb=1.2.0"
    tag "${sample_id}"

    publishDir params.outdir, mode: 'copy', pattern: '*.csv'

    input:
    tuple val(sample_id), path(flycodes), path(sequences)
    path(reference)

    output:
    tuple val(sample_id), path("clusters/*.fastq.gz"), path("references/*.fasta"), path("references.fasta"), emit:grouped_reads
    tuple path("${sample_id}_reads.csv"), path("${sample_id}_clusters.csv"), path("${sample_id}_mapped_reads.csv"), path("${sample_id}_mapped_reads_filtered.csv"), emit:csv

    script:
    """
    group_by_flycodes.py \
        --sample_id ${sample_id} \
        --flycodes ${flycodes} \
        --sequence ${sequences} \
        --reference_seq ${reference} \
        --reference_flycode ${params.reference_flycode}
    """
}