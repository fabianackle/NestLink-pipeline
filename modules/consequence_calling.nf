process CONSEQUENCE_CALLING {
    conda "bioconda::dnaio=1.2.3 bioconda::bcftools=1.23"

    tag "${sample_id}"

    input:
    tuple val(sample_id), path(variants), path(references)

    script:
    """
    write_gff.py \
        --references ${references} \
        --orf_5p ${params.orf_5p} \
        --orf_3p ${params.orf_3p} \
        --output ${sample_id}.gff

    bcftools csq \
        -f ${references} \
        -g ${sample_id}.gff \
        ${variants} \
        > ${sample_id}_consequences.vcf

    bcftools query \
        -f '[%CHROM,%FILTER,%TBCSQ{0}\n]' \
        ${sample_id}_consequences.vcf \
        > ${sample_id}_consequences.csv
    """
}
