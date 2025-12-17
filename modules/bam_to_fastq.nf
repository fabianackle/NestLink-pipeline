process BAM_TO_FASTQ {
    conda "bioconda::samtools=1.22.1 conda-forge::pigz=2.8"
    container "${ workflow.containerEngine == 'apptainer' ?
        'oras://community.wave.seqera.io/library/samtools_pigz:3e26f333629e20af' :
        'community.wave.seqera.io/library/samtools_pigz:243083689f8d2ea0' }"

    tag "${basecalled.baseName}"

    input:
    tuple val(sample_id), path(basecalled)

    output:
    tuple val(sample_id), path("${sample_id}.fastq.gz"), emit: fastq_gz

    script:
    """
    samtools fastq \
        --threads ${task.cpus} \
        ${basecalled} \
        | pigz -p ${task.cpus} > ${sample_id}.fastq.gz
    """
}
