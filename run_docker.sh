#!/bin/bash

docker run \
    -v "$PWD/.nextflow:/workspace/.nextflow" \
    -v "$PWD/output:/output" \
    -v "$PWD/input:/input" \
    -w /workspace \
    -v "$PWD:/workspace" \
    afdb-kit nextflow run workflow.nf -resume
