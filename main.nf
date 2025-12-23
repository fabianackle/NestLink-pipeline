#!/usr/bin/env nextflow

/* Validation */
include { validateParameters  } from 'plugin/nf-schema'
include { paramsSummaryLog    } from 'plugin/nf-schema'

/* Modules */
include { BAM_TO_FASTQ        } from './modules/bam_to_fastq.nf'
include { DORADO_ALIGNER      } from './modules/dorado_aligner.nf'
include { FILTER_READS        } from './modules/filter_reads.nf'
include { EXTRACT_BARCODES    } from './modules/extract_barcodes.nf'
include { GROUP_BY_BARCODES   } from './modules/group_by_barcodes.nf'
include { REMAP_BAM           } from './modules/remap_bam.nf'
include { DORADO_CONSENSUS    } from './modules/dorado_consensus.nf'
include { VARIANT_CALLING     } from './modules/variant_calling.nf'
include { CONSEQUENCE_CALLING } from './modules/consequence_calling.nf'
include { DUCKDB              } from './modules/duckdb.nf'

/* Workflows */
workflow {
    log.info(
        """
        ┌─────────────────────────────────────────────┐
        │ O N T   C O N S E N S U S   P I P E L I N E │
        │ by Fabian Ackle                             │  
        └─────────────────────────────────────────────┘
        """.stripIndent()
    )

    validateParameters()
    log.info(paramsSummaryLog(workflow))

    basecalled_ch = channel.fromPath(params.data, checkIfExists: true).map { file -> tuple(file.baseName, file) }
    reference_ch = channel.fromPath(params.reference, checkIfExists: true)

    consensus(basecalled_ch, reference_ch)
}

workflow consensus {
    take:
    basecalled_ch
    reference_ch

    main:
    BAM_TO_FASTQ(basecalled_ch)
    FILTER_READS(BAM_TO_FASTQ.out.fastq_gz)
    EXTRACT_BARCODES(FILTER_READS.out.reads)

    barcodes_ch = EXTRACT_BARCODES.out.barcodes.combine(reference_ch)
    GROUP_BY_BARCODES(barcodes_ch)

    align_ch = basecalled_ch.combine(reference_ch)
    DORADO_ALIGNER(align_ch)

    remap_ch = DORADO_ALIGNER.out.alignment.join(GROUP_BY_BARCODES.out.barcode_map)
    REMAP_BAM(remap_ch)

    consensus_ch = REMAP_BAM.out.bam.join(GROUP_BY_BARCODES.out.references)
    DORADO_CONSENSUS(consensus_ch)

    variant_ch = DORADO_CONSENSUS.out.consensus.combine(reference_ch)
    VARIANT_CALLING(variant_ch)

    consquence_ch = DORADO_CONSENSUS.out.variants.join(GROUP_BY_BARCODES.out.references)
    CONSEQUENCE_CALLING(consquence_ch)

    duck_ch = GROUP_BY_BARCODES.out.csv.join(VARIANT_CALLING.out.variants_db)
    DUCKDB(duck_ch)
}
