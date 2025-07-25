#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=3850
#SBATCH --time=1-00:00:00
#SBATCH --output=log/main_%j
##SBATCH --mail-user=user@example.com
##SBATCH --mail-type=BEGIN,END 

module load singularityce
module load mamba
eval "$(conda shell.bash hook)"
conda activate nextflow

nextflow run main.nf \
    -profile cluster \
    -params-file MsbA.json \
    -resume \
    -with-report \
    -with-timeline \
    -ansi-log false
