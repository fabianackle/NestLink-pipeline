#!/usr/bin/env nextflow

/* Define parameters */
params.data = "$projectDir/data/barcode15.fastq.gz"
params.name = "barcode15"
params.outdir = "$projectDir/results"

/* Print pipeline info */
log.info """
    =================================
    N E S T L I N K   P I P E L I N E
    =================================
    """
    .stripIndent()

/* Processes */

process ExtractFlycodes {
    cpus 8
    conda "bioconda::cutadapt=4.9 bioconda::seqkit=2.8.2"
    tag "Cutadapt on $fastq_gz"

    input:
    tuple val(sampleName), path(fastq_gz)

    output:
    tuple val(sampleName), path("${sampleName}_flycodes.fasta")

    script:
    """
    cutadapt -j $task.cpus \
        -g CATCATCATTCAAGa...ctggaggtgctcttt \
        --error-rate 0.2 \
        --minimum-length 30 --maximum-length 50 \
        --discard-untrimmed \
        --fasta ${fastq_gz} > fwd.fasta
    
    cutadapt -j $task.cpus \
        -g aaagagcacctccag...tCTTGAATGATGATG \
        --error-rate 0.2 \
        --minimum-length 30 --maximum-length 50 \
        --discard-untrimmed \
        --fasta ${fastq_gz} > rev.fasta

    seqkit seq --threads $task.cpus \
        --reverse --complement rev.fasta \
        --out-file rev_rc.fasta

    cat fwd.fasta rev_rc.fasta > ${sampleName}_flycodes.fasta
    """
}

process ClusterFlycodes {
    cpus 8
    conda "bioconda::vsearch=2.28"
    tag "Vsearch on $flycodes"

    input:
    tuple val(sampleName), path(flycodes)

    output:
    tuple val(sampleName), path("${sampleName}_flycode_centroids.fasta"), path("${sampleName}_flycode_clusters")

    script:
    """
    mkdir ${sampleName}_flycode_clusters
    vsearch -cluster_size \
        $flycodes \
        -id 0.90 \
        -sizeout \
        -clusterout_id \
        -centroids ${sampleName}_flycode_centroids.fasta \
        -clusters ${sampleName}_flycode_clusters/cluster.fasta
    """
}

process GroupSequences {
    cpus 1
    conda "bioconda::dnaio=1.2.1 conda-forge::pandas=2.2.2"
    tag "group_by_flycodes.py on $clusters with $centroids"

    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sampleName), path(centroids), path(clusters)

    output:
    tuple path("${sampleName}_cluster.csv"), path("${sampleName}_reads.csv")

    script:
    """
    mkdir ${sampleName}_binned
    group_by_flycodes.py \
        --centroids $centroids \
        --clusters $clusters \
        --outdir ${sampleName}_binned

    mv clusters.csv ${sampleName}_cluster.csv
    mv reads.csv ${sampleName}_reads.csv
    """
}

/* Workflow */
workflow {
    Channel
        .fromPath(params.data)
        .set { basecalled_ch }
    Channel
        .of(params.name)
        .set { name_ch}

    flycodes_ch = ExtractFlycodes(name_ch.combine(basecalled_ch))
    clusters_ch = ClusterFlycodes(flycodes_ch)
    GroupSequences(clusters_ch)
}