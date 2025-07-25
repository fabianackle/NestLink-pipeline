process EXTRACT_FLYCODES {
    cpus 8
    memory '4 GB'
    time '60m'
    conda "bioconda::cutadapt=5.1 bioconda::seqkit=2.10.0"
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(fastq_gz)

    output:
    tuple val(sample_id), path("${sample_id}_flycodes.fasta"), emit: flycodes

    script:
    """
    extract_flycode_adapter_rc=\$(echo '${params.extract_flycode_adapter}' | tr 'ACGTacgt.' 'TGCAtgca.' | rev)

    cutadapt \
        -j $task.cpus \
        -g $params.extract_flycode_adapter \
        --error-rate 0.1 \
        --minimum-length 30 --maximum-length 36 \
        --discard-untrimmed \
        --fasta ${fastq_gz} > ${sample_id}_flycodes_fwd.fasta

    cutadapt \
        -j $task.cpus \
        -g \$extract_flycode_adapter_rc \
        --error-rate 0.1 \
        --minimum-length 30 --maximum-length 36 \
        --discard-untrimmed \
        --fasta ${fastq_gz} > ${sample_id}_flycodes_rc.fasta

    seqkit seq \
        --threads $task.cpus \
        --seq-type dna \
        --reverse --complement ${sample_id}_flycodes_rc.fasta \
        --out-file ${sample_id}_flycodes_rev.fasta

    cat ${sample_id}_flycodes_fwd.fasta ${sample_id}_flycodes_rev.fasta > ${sample_id}_flycodes_concat.fasta

    seqkit rmdup \
        --by-name ${sample_id}_flycodes_concat.fasta \
        --out-file ${sample_id}_flycodes.fasta
    """
}