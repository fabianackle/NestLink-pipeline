process DORADO_CONSENSUS {
    tag "${sample_id}"

    publishDir params.outdir, mode: 'copy', pattern: '*.fastq.gz'

    input:
    tuple val(sample_id), path(bam), path(bai), path(reference)

    output:
    tuple val(sample_id), path("${sample_id}_polished.fastq.gz"), emit: consensus

    script:
    """
    dorado polish ${bam} ${reference} \
        --qualities \
        --ignore-read-groups \
        --batchsize 250 \
        --infer-threads ${ task.cpus.intdiv(2) } \
        --threads ${task.cpus} \
        --bacteria \
        > ${sample_id}_polished.fastq

    gzip ${sample_id}_polished.fastq   
    """
}
