process DUCKDB {
    conda "conda-forge::duckdb-cli=1.4.1"
    container "${ workflow.containerEngine == 'apptainer' ?
        'oras://community.wave.seqera.io/library/duckdb-cli:1.4.1--666aeb32eabd82a9' :
        'community.wave.seqera.io/library/duckdb-cli:1.4.1--d924e68d63392ee0' }"

    tag "${sample_id}"

    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample_id), path(reads_csv), path(hq_barcodes_csv), path(mapped_reads_csv), path(mapped_reads_filtered_csv), path(variants_csv)

    output:
    tuple val(sample_id), path("${sample_id}.db"), emit: database

    script:
    """
    cat << EOF > query.sql

    CREATE TABLE reads AS
    FROM read_csv('${reads_csv}',
        columns = {
            'read_id': 'UUID',
            'barcode': 'VARCHAR',
            'is_valid_barcode ': 'BOOLEAN'
        });

    CREATE TABLE hq_barcodes AS
    FROM read_csv('${hq_barcodes_csv}',
        columns = {
            'barcode_id': 'UUID',
            'barcode': 'VARCHAR'
        });

    CREATE TABLE reads_mapped AS
    FROM read_csv('${mapped_reads_csv}',
        columns = {
            'read_id': 'UUID',
            'barcode_id': 'UUID',
            'cigar': 'VARCHAR',
            'edit_distance': UINT8
        });

    CREATE TABLE reads_mapped_filtered AS
    FROM read_csv('${mapped_reads_filtered_csv}',
        columns = {
            'read_id': 'UUID',
            'barcode_id': 'UUID',
            'cigar': 'VARCHAR',
            'edit_distance': UINT8
        });

    CREATE TABLE variants AS
    FROM read_csv('${variants_csv}',
        columns = {
            'barcode_id': 'UUID',
            'barcode': 'VARCHAR',
            'variant_type': VARCHAR,
            'position': 'UINT16',
            'reference_aa': VARCHAR,
            'variant_aa': 'VARCHAR'
        });

    CREATE VIEW reads_summary AS
    FROM reads
    JOIN reads_mapped USING(read_id);

    -- Variant summary: barcodes per variant, reads per variant, reads per barcode
    CREATE VIEW variants_summary AS
    SELECT
        "position",
        reference_aa || "position" || variant_aa AS variant,
        COUNT(DISTINCT barcode) AS barcodes,
        COUNT(read_id) AS reads,
        COUNT(read_id) / COUNT(DISTINCT barcode) AS reads_per_barcode
    FROM variants
    JOIN reads_mapped USING (barcode_id)
    WHERE
        variant_type = 'change'
    GROUP BY
        "position", variant
    ORDER BY
        "position" ASC, variant ASC;

    EOF

    duckdb < query.sql ${sample_id}.db
    """
}
