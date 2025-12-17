process FILTER_READS {
    conda "bioconda::filtlong=0.3.1 conda-forge::pigz=2.8"
    container "${ workflow.containerEngine == 'apptainer' ?
        'oras://community.wave.seqera.io/library/filtlong_pigz:fab2e27db205f4da' :
        'community.wave.seqera.io/library/filtlong_pigz:a11566667039b4fe' }"

    tag "${sample_id}"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_filtered.fastq.gz"), emit: reads

    script:
    """
    filtlong \
        --min_mean_q ${params.filter_quality} \
        --min_length ${params.filter_min_length} \
        --max_length ${params.filter_max_length} \
        ${reads} \
        | pigz -p ${task.cpus} > ${sample_id}_filtered.fastq.gz
    """
}
