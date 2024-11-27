#!/usr/bin/env bash

nextflow run main.nf \
    -profile standard \
    -params-file params.json \
    -resume \
    -with-report \
    -with-timeline