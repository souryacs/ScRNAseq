#!/bin/bash

##===================
## TCR data analysis using MixCR
## TCR data is generated from BD-Rhapsody protocol
## check https://github.com/milaboratory/mixcr
##===================

## Output directory
## user can edit this path
OutDir='MixCR_Out_TCR'
mkdir -p $OutDir

## we assume paired end fastq file
## user can edit the file path as required
FastqFile_R1='FastqFile_R1_001.fastq.gz'
FastqFile_R2='FastqFile_R2_001.fastq.gz'

## assumimng human reference genome
## The string "outprefix_TCR" will be the prefix string for every file
## user can edit accordingly
mixcr analyze bd-rhapsody-human-tcr-full-length \
	-t 8 \
	-f \
	$FastqFile_R1 \
	$FastqFile_R2 \
	$OutDir'/outprefix_TCR'

mixcr exportClones \
	$OutDir'/outprefix_TCR.clna' \
	$OutDir'/outprefix_TCR_table_exportClones.tsv' 

## export the alignments
mixcr exportAlignments \
	$OutDir'/outprefix_TCR.clna' \
	$OutDir'/outprefix_TCR_clna_table_exportAlignments.tsv' 

mixcr exportAlignments \
	$OutDir'/outprefix_TCR.vdjca' \
	$OutDir'/outprefix_TCR_vdjca_table_exportAlignments.tsv' 

mixcr exportAlignments \
	$OutDir'/outprefix_TCR.refined.vdjca' \
	$OutDir'/outprefix_TCR_refined_vdjca_table_exportAlignments.tsv' 


