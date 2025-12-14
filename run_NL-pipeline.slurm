#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=3850
#SBATCH --time=1-00:00:00
#SBATCH --output=log/main_%j
##SBATCH --mail-user=user@example.com
##SBATCH --mail-type=BEGIN,END 

module load miniforge3
module load apptainer
eval "$(conda shell.bash hook)"
conda activate nextflow

nextflow run main.nf \
    -profile s3it,apptainer \
    -params-file params.json \
    -resume \
    -with-report \
    -with-timeline \
    -ansi-log false
