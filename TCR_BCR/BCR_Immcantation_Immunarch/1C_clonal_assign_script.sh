#!/bin/bash

##================
## script 1C: BCR analysis using immcantation pipeline
## Assigns Ig sequences into clonally related lineages and builds full germline sequences.
## check https://immcantation.readthedocs.io/en/stable/docker/pipelines.html#xpipeline
##================

## Check the tutorial and download the "immcantation_suite-4.1.0.sif" (or the latest version)
## It is immcantation executable
## user needs to edit
imm_exec_path = "/path/to/immcantation_suite-4.1.0.sif"

## output directory to store the results
## user needs to edit
OutDir="/path/to/outdir"

## check the output file *threshold-values.tab to get the threshold
## user can edit if needed, according to the output file 
DIST=0.174

## we need to move to the directory containing the example data
## as said in the tutorial
cd $OutDir'/out_1_change_o'
singularity exec $imm_exec_path changeo-clone \
    -d newRun_sampledata_genotyped.tsv \
    -x ${DIST} \
    -n newRun_sampledata \
    -o . \
    -p 4 

