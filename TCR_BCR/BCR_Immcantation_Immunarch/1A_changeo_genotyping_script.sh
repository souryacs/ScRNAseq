#!/bin/bash

##================
## script 1A: BCR analysis using immcantation pipeline
## infer V segment genotypes using TIgGER.
## check https://immcantation.readthedocs.io/en/stable/docker/pipelines.html#xpipeline
##================

## Check the tutorial and download the "immcantation_suite-4.1.0.sif" (or the latest version)
## It is immcantation executable
## user needs to edit
imm_exec_path = "/path/to/immcantation_suite-4.1.0.sif"

## output directory to store the results
## user needs to edit
OutDir="/path/to/outdir"

## execute script
cd $OutDir'/out_1_change_o'
singularity exec $imm_exec_path tigger-genotype \
    -d filtered_contig_db-pass.tsv \
    -n newRun_sampledata \
    -o . \
    -p 4 

