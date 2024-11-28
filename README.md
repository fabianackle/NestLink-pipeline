# NestLink-pipeline
NestLink-pipeline is a pipeline for processing [NestLink libraries](https://www.nature.com/articles/s41592-019-0389-8) sequenced by nanopore sequencing. Reads are binned according to their flycodes (UMIs). Accurate consensus sequences are calculated using Medaka. Variants are called with the pipeline, resulting in a flycode assignment table that links protein variants to their respective set of flycodes.

> [!WARNING]
> NestLink-pipeline is still in development. Certain library-specific strings are still hard-coded in `main.nf` and have to be edited before running the pipeline.

## Requirements
### Local and cluster execution
- Mamba/ Conda ([https://conda-forge.org/](https://conda-forge.org/))
- Nextflow ([Installation guide](https://www.nextflow.io/docs/latest/install.html))
- mini_align ([mini_align.sh](https://raw.githubusercontent.com/nanoporetech/pomoxis/master/scripts/mini_align) placed in `projectDir/bin/`)
### Cluster execution only
- Slurm workflow manager
- Singularity

## Running the pipeline on the s3it cluster
1. Clone the repository.
2. Edit the params.json file, specify the nanopore reads (bam) and reference sequence.
3. Run the pipeline:
`sbatch run_NL-pipeline.slurm`

## Running the pipeline locally
> [!NOTE]
> Consensus sequence generation with medaka has to be run manually.
1. Prepare the pipeline as described above.
2. Run the pipeline:
`bash run_NL-pipeline.sh`