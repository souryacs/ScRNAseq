#!/usr/bin/env Rscript

##================
## script 4: BCR analysis using immcantation pipeline
## plotting linage tree
##================
library(alakazam)
library(ape)
library(dplyr)

## output directory to store the results
## user needs to edit
OutDir="/path/to/outdir"

inpdir <- paste0(OutDir, '/out_2_clonal_groups')

##read in the data
db <- readIgphyml(paste0(inpdir, "/filtered_contig_heavy_germ-pass_igphyml-pass.tab"), format="phylo", branches="mutations")

for (i in 1:length(db$trees)) {
	png(paste0(inpdir, "/out_lineage_tree_", i, ".png"), width=8,height=6,unit="in",res=300)
	plot(db$trees[[i]],show.node.label=TRUE)
	add.scale.bar(length=5)
	dev.off()
}

