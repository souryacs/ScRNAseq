#!/bin/bash

##===================
## BCR data analysis using MixCR
## BCR data is generated from BD-Rhapsody protocol
## check https://github.com/milaboratory/mixcr
##===================

## Output directory
## user can edit this path
OutDir='MixCR_Out_BCR'
mkdir -p $OutDir

## we assume paired end fastq file
## user can edit the file path as required
FastqFile_R1='FastqFile_R1_001.fastq.gz'
FastqFile_R2='FastqFile_R2_001.fastq.gz'

## assumimng human reference genome
## The string "outprefix_BCR" will be the prefix string for every file
## user can edit accordingly
mixcr analyze bd-rhapsody-human-bcr-full-length \
	-t 8 \
	-f \
	$FastqFile_R1 \
	$FastqFile_R2 \
	$OutDir'/outprefix_BCR'

mixcr exportClones \
	$OutDir'/BC_S2_L003.clna' \
	$OutDir'/BC_S2_L003_clna_table_exportClones.tsv' 

## export the alignments
mixcr exportAlignments \
	$OutDir'/outprefix_BCR.clna' \
	$OutDir'/outprefix_BCR_clna_table_exportAlignments.tsv' 

mixcr exportAlignments \
	$OutDir'/outprefix_BCR.vdjca' \
	$OutDir'/outprefix_BCR_vdjca_table_exportAlignments.tsv' 

mixcr exportAlignments \
	$OutDir'/outprefix_BCR.refined.vdjca' \
	$OutDir'/outprefix_BCR_refined_vdjca_table_exportAlignments.tsv' 


