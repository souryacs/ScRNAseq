#!/bin/bash

##================
## script to run Callisto quantification on BD Rhapsody generated fastq files
## For details of Kallisto, check the manual
## https://pachterlab.github.io/kallisto/manual.html
##================

## Kallisto compatible transcriptome index file for the reference genome Human
## user can edit the file path
Indexfile='path/to/Kallisto/Indices/Ensemble_Transcriptome_V96/homo_sapiens/transcriptome.idx'

## output directory - user can edit
OutDir='/path/to/Kallisto_Quant_Out'

## GTF file for reference genome - user can edit the file path
GTFFile='/path/to/GTF/hg38/gencode.v38.annotation.gtf'

## chromosome size file corresponding to the reference genome
chrsizefile='/path/to/GRCh38-PhiX-gencodev29/chrNameLength.txt'

## fastq files
## we are assuming paired-end read fastq files
## user needs to edit the paths
fastq_1_R1="File1_R1_001.fastq.gz"
fastq_1_R2="File1_R2_001.fastq.gz"

## script to run the quantification method
kallisto quant \
    -i $Indexfile \
    -o $OutDir \
    --genomebam \
    -t 8 \
    -g $GTFFile \
    -c $chrsizefile \
    $fastq_1_R1 \
    $fastq_1_R2

