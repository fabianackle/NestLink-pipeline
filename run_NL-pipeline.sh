#!/usr/bin/env bash

nextflow run main.nf \
    -profile standard \
    -params-file MsbA.json \
    -resume \
    -with-report \