process EXTRACT_BARCODES {
    conda "bioconda::cutadapt=5.2 bioconda::seqkit=2.10.1"
    container "${ workflow.containerEngine == 'apptainer' ?
        'oras://community.wave.seqera.io/library/cutadapt_seqkit:747e3b8919e43ec5' :
        'community.wave.seqera.io/library/cutadapt_seqkit:e9bcc3ce5550445c' }"

    tag "${sample_id}"

    input:
    tuple val(sample_id), path(fastq_gz)

    output:
    tuple val(sample_id), path("${sample_id}_barcodes.fasta"), emit: barcodes

    script:
    """
    extract_barcode_adapter_rc=\$(echo '${params.extract_barcode_adapter}' | tr 'ACGTacgt.' 'TGCAtgca.' | rev)

    cutadapt \
        -j ${task.cpus} \
        -g ${params.extract_barcode_adapter} \
        --error-rate 0.1 \
        --maximum-length 50 \
        --discard-untrimmed \
        --fasta ${fastq_gz} > ${sample_id}_barcodes_fwd.fasta

    cutadapt \
        -j ${task.cpus} \
        -g \$extract_barcode_adapter_rc \
        --error-rate 0.1 \
        --maximum-length 50 \
        --discard-untrimmed \
        --fasta ${fastq_gz} > ${sample_id}_barcodes_rc.fasta

    seqkit seq \
        --threads ${task.cpus} \
        --seq-type dna \
        --reverse --complement ${sample_id}_barcodes_rc.fasta \
        --out-file ${sample_id}_barcodes_rev.fasta

    cat ${sample_id}_barcodes_fwd.fasta ${sample_id}_barcodes_rev.fasta > ${sample_id}_barcodes_concat.fasta

    seqkit rmdup \
        --by-name ${sample_id}_barcodes_concat.fasta \
        --out-file ${sample_id}_barcodes.fasta
    """
}
