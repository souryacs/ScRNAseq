#!/bin/bash 

##================
## script 3: BCR analysis using immcantation pipeline
## Building linage trees
## check https://immcantation.readthedocs.io/en/stable/tutorials/10x_tutorial.html
##================

## Check the tutorial and download the "immcantation_suite-4.1.0.sif" (or the latest version)
## It is immcantation executable
## user needs to edit
imm_exec_path = "/path/to/immcantation_suite-4.1.0.sif"

## output directory to store the results
## user needs to edit
OutDir="/path/to/outdir"

cd $OutDir'/out_2_clonal_groups'

## To run IgPhyML
singularity exec $imm_exec_path BuildTrees.py \
    -d filtered_contig_heavy_germ-pass.tsv \
    --minseq 3 \
    --clean all \
    --igphyml \
    --collapse \
    --nproc 2 \
    --asr 0.9


