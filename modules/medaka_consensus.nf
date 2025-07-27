process MEDAKA_CONSENSUS {
    cpus 8
    memory '16 GB'
    time {'60m * task.attempt'}
    errorStrategy 'retry'
    maxRetries 2
    clusterOptions '--gpus=1'
    tag "${sample_id}"

    publishDir params.outdir, mode: 'copy', pattern: '*.fasta'

    input:
    tuple val(sample_id), path(references), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}_assembly.fasta"), emit: consensus
    tuple path("medaka_version.log"), path("medaka_interference.log"), path("medaka_sequence.log"), emit: log

    script:
    """
    medaka --version > medaka_version.log

    medaka inference \
        --batch 200 --threads 2 \
        --model ${params.medaka_dorado_model}  \
        ${bam} results.contigs.hdf \
        2> medaka_interference.log

    medaka sequence \
        results.contigs.hdf ${references} ${sample_id}_assembly.fasta \
        2> medaka_sequence.log
    """
}