# NestLink-pipeline

> [!NOTE]
> With this adapted version of the NestLink-pipeline flycodes can be extracted and clustered.

## Requirements
- Mamba/Conda ([https://conda-forge.org/](https://conda-forge.org/))
- Nextflow `mamba create -n nextflow bioconda::nextflow`

## Running the pipeline
1. Clone the repository.
2. Merge the basecalled fastq.gz by concatenating the individual files. `cat *.fastq.gz > barcodeNN.fastq.gz `
2. Place the basecalled sequencing data and the reference sequence into `projectDir/data/`. Adjust the filenames in `main.nf` (`params.data` and `params.name`).
3. Run the main workflow of the pipeline:
`nextflow run main.nf`