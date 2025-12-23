process DORADO_CONSENSUS {
    container "nanoporetech/dorado:shaf2aed69855de85e60b363c9be39558ef469ec365"

    tag "${sample_id}"

    publishDir params.outdir, mode: 'copy', pattern: '*.fastq.gz'

    input:
    tuple val(sample_id), path(bam), path(bai), path(reference)

    output:
    tuple val(sample_id), path("${sample_id}_consensus.fastq.gz"), emit: consensus
    tuple val(sample_id), path("${sample_id}_variants.vcf"), emit: variants

    script:
    """
    dorado polish ${bam} ${reference} \
        --qualities \
        --vcf \
        --ignore-read-groups \
        --batchsize 250 \
        --infer-threads ${ task.cpus.intdiv(2) } \
        --threads ${task.cpus} \
        ${params.polish_bacteria ? '--bacteria' : ''} \
        -o .

    gzip consensus.fastq 
    mv consensus.fastq.gz ${sample_id}_consensus.fastq.gz
    mv variants.vcf ${sample_id}_variants.vcf 
    """
}
