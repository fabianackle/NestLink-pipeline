#!/usr/bin/env bash

nextflow run main.nf \
    -params-file kdelr.json \
    -resume \
    -with-report \