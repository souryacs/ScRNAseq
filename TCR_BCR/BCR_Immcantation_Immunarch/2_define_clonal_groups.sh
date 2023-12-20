#!/bin/bash

##================
## script 2: BCR analysis using immcantation pipeline
## Define heavy chain clones, Split them in light chains and construct germline V and J sequences
## check https://immcantation.readthedocs.io/en/stable/tutorials/10x_tutorial.html
##================

## Check the tutorial and download the "immcantation_suite-4.1.0.sif" (or the latest version)
## It is immcantation executable
## user needs to edit
imm_exec_path = "/path/to/immcantation_suite-4.1.0.sif"

## output directory to store the results
## user needs to edit
OutDir="/path/to/outdir"

cd $OutDir

## first extract contigs with heavy chain and productive status
## somehow this file is not generated from the first script
awk -F'[\t]' '((NR==1) || (($4=="T") && ($17=="IGH")))' filtered_contig_db-pass.tsv > filtered_contig_heavy_productive-T.tsv

##===== now proceed with the pipeline mentioned in the document
## check https://immcantation.readthedocs.io/en/stable/tutorials/10x_tutorial.html

## define heavy chain clones
## user can check and edit the parameters if needed
singularity exec $imm_exec_path DefineClones.py \
    -d filtered_contig_heavy_productive-T.tsv \
    --act set \
    --model ham \
    --norm len \
    --dist 0.09 \
    --outname filtered_contig_heavy

## split heavy chain clones with different light chains
singularity exec $imm_exec_path light_cluster.py \
    -d filtered_contig_heavy_clone-pass.tsv \
    -e filtered_contig_light_productive-T.tsv \
    -o filtered_contig_heavy_clone-light.tsv

## reconstruct heavy chain germline V and J sequences
## uses the fasta sequences corresponding to the reference genome from /usr/local directory
## saved during installation of this package
singularity exec $imm_exec_path CreateGermlines.py \
    -d filtered_contig_heavy_clone-light.tsv \
    -g dmask \
    --cloned \
    -r /usr/local/share/germlines/imgt/human/vdj/imgt_human_IGHV.fasta \
        /usr/local/share/germlines/imgt/human/vdj/imgt_human_IGHD.fasta \
        /usr/local/share/germlines/imgt/human/vdj/imgt_human_IGHJ.fasta \
    --outname filtered_contig_heavy


