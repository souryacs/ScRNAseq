#!/bin/bash

##================
## script 1: BCR analysis using immcantation pipeline
## check https://immcantation.readthedocs.io/en/stable/tutorials/10x_tutorial.html
## Assigns new annotations and infers clonal relationships to 10x Genomics single-cell V(D)J data output by Cell Ranger.
## 10X genomic VD(J) data
## Assign V, D, and J genes and define clonal groups
##================

## Check the tutorial and download the "immcantation_suite-4.1.0.sif" (or the latest version)
## It is immcantation executable
## user needs to edit
imm_exec_path = "/path/to/immcantation_suite-4.1.0.sif"

## output directory to store the results
## user needs to edit
OutDir="/path/to/outdir"

## execute script
## Assigns new annotations and infers clonal relationships
cd $OutDir
singularity exec $imm_exec_path changeo-10x \
        -s filtered_contig.fasta \
        -a filtered_contig_annotations.csv \
        -o . \
        -g human \
        -t ig \
        -x 0.1

