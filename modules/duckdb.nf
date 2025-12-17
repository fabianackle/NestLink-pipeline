process DUCKDB {
    conda "conda-forge::duckdb-cli=1.4.1"
    container "${ workflow.containerEngine == 'apptainer' ?
        'oras://community.wave.seqera.io/library/duckdb-cli:1.4.1--666aeb32eabd82a9' :
        'community.wave.seqera.io/library/duckdb-cli:1.4.1--d924e68d63392ee0' }"

    tag "${sample_id}"

    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample_id), path(reads_csv), path(clusters_csv), path(mapped_reads_csv), path(mapped_reads_filtered_csv), path(variants_csv)

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

    CREATE TABLE clusters AS
    FROM read_csv('${clusters_csv}',
        columns = {
            'cluster_id': 'UUID',
            'barcode': 'VARCHAR'
        });

    CREATE TABLE mapped_reads AS
    FROM read_csv('${mapped_reads_csv}',
        columns = {
            'read_id': 'UUID',
            'cluster_id': 'UUID',
            'cigar': 'VARCHAR',
            'edit_distance': UINT8
        });

    CREATE TABLE mapped_reads_filtered AS
    FROM read_csv('${mapped_reads_filtered_csv}',
        columns = {
            'read_id': 'UUID',
            'cluster_id': 'UUID',
            'cigar': 'VARCHAR',
            'edit_distance': UINT8
        });

    CREATE TABLE variants AS
    FROM read_csv('${variants_csv}',
        columns = {
            'cluster_id': 'UUID',
            'barcode': 'VARCHAR',
            'variant_type': VARCHAR,
            'position': 'UINT16',
            'reference_aa': VARCHAR,
            'variant_aa': 'VARCHAR'
        });

    CREATE VIEW barcodes AS
    FROM reads
    JOIN mapped_reads USING(read_id);

    -- Variant summary: barcodes per variant, reads per variant, reads per barcode
    CREATE VIEW variant_summary AS
    SELECT
        "position",
        reference_aa || "position" || variant_aa AS variant,
        COUNT(DISTINCT barcode) AS barcodes,
        COUNT(read_id) AS reads,
        COUNT(read_id) / COUNT(DISTINCT barcode) AS reads_per_barcode
    FROM variants
    JOIN mapped_reads USING (cluster_id)
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
