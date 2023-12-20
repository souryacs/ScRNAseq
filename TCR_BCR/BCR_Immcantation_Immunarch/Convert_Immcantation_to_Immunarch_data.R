#!/usr/bin/env Rscript

##=======================
# script to convert Immcantation data to Immnarch data
##=======================
library(data.table)
options(scipen = 10)
options(datatable.fread.datatable=FALSE)

## output directory to store the results
## user needs to edit
OutDir="/path/to/outdir"

basedir <- paste0(OutDir, '/out_1_change_o')

## file from the immcantation pipeline
## with genotyoe infor
inpfile <- paste0(basedir, '/newRun_sampledata_germ-pass.tsv')
inpdata <- data.table::fread(inpfile, header=T)
cat(sprintf("\n number of entries : %s ", nrow(inpdata)))

outdata <- cbind.data.frame(inpdata$umi_count, inpdata$junction, 
							inpdata$junction_aa, inpdata$v_call_10x, inpdata$d_call_10x, inpdata$j_call_10x, 
							inpdata$v_sequence_end, inpdata$d_sequence_start, inpdata$d_sequence_end, 
							inpdata$j_sequence_start, inpdata$sequence)
colnames(outdata) <- c('Clones', 'CDR3.nt', 'CDR3.aa', 'V.name', 'D.name', 'J.name', 
						'V.end', 'D.start', 'D.end', 'J.start', 'Sequence')

## actually the "CDR3.nt" and "CDR3.aa" fields contain flanking codons
## for CDR3.nt, remove first three and last three nucleotides
## for CDR3.aa, remove first and last amino acids
outdata$CDR3.nt <- substr(as.character(outdata$CDR3.nt), 4, nchar(as.character(outdata$CDR3.nt)) - 3)
outdata$CDR3.aa <- substr(as.character(outdata$CDR3.aa), 2, nchar(as.character(outdata$CDR3.aa)) - 1)

## get the Proportion of clones
sumclonecnt <- sum(outdata$Clones)
outdata$Proportion <- (outdata$Clones * (1.0 / sumclonecnt))

## number of inserted nucleotides (N-nucleotides) at V-J junction (-1 for receptors with VDJ recombination)
VJ.ins <- rep(-1, nrow(outdata))
## number of inserted nucleotides (N-nucleotides) at V-D junction (-1 for receptors with VJ recombination)
VD.ins <- rep(-1, nrow(outdata))
## number of inserted nucleotides (N-nucleotides) at D-J junction (-1 for receptors with VJ recombination)
DJ.ins <- rep(-1, nrow(outdata))

## get the indices with VDJ recombination
VDJ_idx <- which(inpdata$d_sequence_start > inpdata$v_sequence_end)
if (length(VDJ_idx) > 0) {
	VD.ins[VDJ_idx] <- (inpdata[VDJ_idx, c('d_sequence_start')] - inpdata[VDJ_idx, c('v_sequence_end')] - 1)
	DJ.ins[VDJ_idx] <- (inpdata[VDJ_idx, c('j_sequence_start')] - inpdata[VDJ_idx, c('d_sequence_end')] - 1)
}

VJ_idx <- setdiff(seq(1,nrow(inpdata)), VDJ_idx)
if (length(VJ_idx) > 0) {
	VJ.ins[VJ_idx] <- (inpdata[VJ_idx, c('j_sequence_start')] - inpdata[VJ_idx, c('v_sequence_end')] - 1)
}

cat(sprintf("\n length VDJ_idx : %s ", length(VDJ_idx)))
cat(sprintf("\n length VJ_idx : %s ", length(VJ_idx)))

## insert in the final data structure
outdata$VJ.ins <- VJ.ins
outdata$VD.ins <- VD.ins
outdata$DJ.ins <- DJ.ins

## reorganize the columns
outdata <- cbind.data.frame(outdata$Clones, outdata$Proportion, 
							outdata$CDR3.nt, outdata$CDR3.aa, outdata$V.name, outdata$D.name, outdata$J.name, 
							outdata$V.end, outdata$D.start, outdata$D.end, 
							outdata$J.start, outdata$VJ.ins, outdata$VD.ins, outdata$DJ.ins, outdata$Sequence)
colnames(outdata) <- c('Clones', 'Proportion', 
						'CDR3.nt', 'CDR3.aa', 'V.name', 'D.name', 'J.name', 
						'V.end', 'D.start', 'D.end', 'J.start', 
						'VJ.ins', 'VD.ins', 'DJ.ins', 'Sequence')

## write the output file
outfile <- paste0(basedir, '/newRun_sampledata_germ-pass_compatible_Immunarch.tsv')
write.table(outdata, outfile, row.names=F, col.names=T, sep="\t", quote=F, append=F)

