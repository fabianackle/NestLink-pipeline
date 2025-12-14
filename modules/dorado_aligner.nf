process DORADO_ALIGNER {
    container "nanoporetech/dorado:shaf2aed69855de85e60b363c9be39558ef469ec365"
    
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(basecalled), path(reference)

    output:
    tuple val(sample_id), path("${sample_id}.aligned.bam"), emit: alignment

    script:
    """
    dorado aligner ${reference} ${basecalled} > ${sample_id}.aligned.bam
    """
}
