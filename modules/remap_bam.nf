process REMAP_BAM {
    conda "bioconda::dnaio=1.2.3 bioconda::pysam=0.23.3 bioconda::samtools=1.22.1"
    container "${ workflow.containerEngine == 'apptainer' ?
        'oras://community.wave.seqera.io/library/dnaio_pysam_samtools:19e11ea6ed78be12' :
        'community.wave.seqera.io/library/dnaio_pysam_samtools:b619b880eeab0208' }"

    tag "${sample_id}"

    input:
    tuple val(sample_id), path(alignment), path(barcode_map)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: bam

    script:
    """
    remap_bam.py \
        --barcode_map ${barcode_map} \
        --input_bam ${alignment} \
        --output_bam ${sample_id}_unsorted.bam

    samtools sort -@ ${task.cpus} ${sample_id}_unsorted.bam -o ${sample_id}.bam
    samtools index -@ ${task.cpus} ${sample_id}.bam
    """
}
