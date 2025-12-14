#!/usr/bin/env bash

nextflow run main.nf \
    -params-file kdelr.json \
    -profile conda
    -resume \
    -with-report \