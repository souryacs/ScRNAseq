#!/bin/bash

##================
## script to run Callisto Bus on BD Rhapsody generated fastq files
## For details of Kallisto, check the manual
## https://pachterlab.github.io/kallisto/manual.html
##================

## Kallisto compatible transcriptome index file for the reference genome Human
## user can edit the file path
Indexfile='path/to/Kallisto/Indices/Ensemble_Transcriptome_V96/homo_sapiens/transcriptome.idx'

## output directory - user can edit
OutDir='/path/to/Kallisto_Bus_Out'

## fastq files
## we are assuming paired-end read fastq files
## user needs to edit the paths
fastq_1_R1="File1_R1_001.fastq.gz"
fastq_1_R2="File1_R2_001.fastq.gz"

## BDWTA: BD-Rhapsody technology
## https://github.com/pachterlab/kallisto/releases
kallisto bus \
    -i $Indexfile \
    -o $OutDir \
    -t 8 \
    -x BDWTA \
    $fastq_1_R1 \
    $fastq_1_R2


