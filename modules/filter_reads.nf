process FILTER_READS {
    cpus 1
    memory '4 GB'
    time '60m'
    conda "bioconda::filtlong=0.2.1 conda-forge::pigz=2.8"
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_filtered.fastq.gz"), emit: reads

    script:
    """
    filtlong \
        --min_mean_q 96.84 \
        --min_length $params.filter_min_length \
        --max_length $params.filter_max_length \
        ${reads} \
        | pigz -p $task.cpus > ${sample_id}_filtered.fastq.gz
    """
}