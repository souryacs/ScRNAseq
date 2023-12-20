#!/bin/bash

## Running Velocyto to produce the .loom files compatible for scVelo

##===========
## parameters
##===========

## output directory
## user can edit
OutDir='velocyto_out'
mkdir -p $OutDir

## GTF file for reference genome - from CellRanger
## assume human reference genome - hg38
## user needs to edit
GTFFile='/path/to/CellRanger/GRCh38/genes/genes.gtf'

## bam file - scRNA-seq - from Cellranger output
## user needs to edit
bamfile='/path/to/possorted_genome_bam1.bam'

## current sample ID
## user needs to edit
sampleID='currSample'

## directory containing Cellranger output genes, barcodes and matrix
## filtered feature barcode matrix
## user needs to edit
CellRangerDir='/path/to/filtered_feature_bc_matrix'

##===============
## execute
##===============

## Usually the "barcodes.tsv" file in "CellRangerDir" is zipped - unzip
gunzip $CellRangerDir'/barcodes.tsv.gz'

## first copy the original bam file into the specified output directory
cp $bamfile $OutDir'/'
cp $bamfile'.bai' $OutDir'/'

# then use this copied bam file
bamfile=$OutDir'/possorted_genome_bam.bam'

# sorted bam file
sortedfile=$OutDir'/cellsorted_possorted_genome_bam.bam'

# use samtools sort to manually sort the file
# https://github.com/velocyto-team/velocyto.py/issues/212
samtools sort -@ 8 -t CB -O BAM -o $sortedfile $bamfile

## convert the cellranger output to .loom file using velocyto
## this loom file will be required in scvelo
velocyto run \
    -o $OutDir \
    -e $sampleID \
    -b $CellRangerDir'/barcodes.tsv' \
    -@ 8 $bamfile $GTFFile

## again gzip the "barcodes.tsv" file in "CellRangerDir" 
gzip $CellRangerDir'/barcodes.tsv'

 

