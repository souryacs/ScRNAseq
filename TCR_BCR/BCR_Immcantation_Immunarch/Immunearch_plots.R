#!/usr/bin/env Rscript

##===============
# script to process BCR data compatible to immunarch pipeline
# for different performance plots
##===============

library(immunarch)
library(ggplot2)

## output directory to store the results
## user needs to edit
OutDir="/path/to/outdir"

## input file containing the immunarch compatible BCR data
inpfile <- paste0(OutDir, '/out_1_change_o/newRun_sampledata_germ-pass_compatible_Immunarch.tsv')

plotdir <- paste0(dirname(inpfile), '/newRun_sampledata_Immunarch_plots')
system(paste("mkdir -p", plotdir))

## load the input data
immdata <- repLoad(inpfile)
names(immdata$data) <- rep('testdata', length(immdata$data))

## number of unique clonotypes
plotfile <- paste0(plotdir, '/num_unique_clonotype.pdf')
exp_vol <- repExplore(immdata$data, .method = "volume")
p1 <- vis(exp_vol)
ggsave(plotfile, plot=p1, width=6, height=6)

## distribution of CDR3 lengths
plotfile <- paste0(plotdir, '/distribution_CDR3_lengths.pdf')
exp_len <- repExplore(immdata$data, .method = "len", .col = "aa")
p1 <- vis(exp_len)
ggsave(plotfile, plot=p1, width=6, height=6)

## distribution of clonotype abundances
plotfile <- paste0(plotdir, '/distribution_clonotype_abundances.pdf')
exp_cnt <- repExplore(immdata$data, .method = "count")
p2 <- vis(exp_cnt)
ggsave(plotfile, plot=p2, width=6, height=6)

##===========
## clonality - estimate the diversity of samples
## kept for future implementation - requires multi sample data
##===========
if (0) {

## The clonal.prop method computes the proportion of repertoire occupied by the pools of cell clones
imm_pr <- repClonality(immdata$data, .method = "clonal.prop")

## The top method considers the most abundant cell clonotypes
imm_top <- repClonality(immdata$data, .method = "top", .head = c(10, 100, 1000, 3000, 10000))
plotfile <- paste0(plotdir, '/top_clonal_proportion.pdf')
p1 <- vis(imm_top)
ggsave(plotfile, plot=p1, width=6, height=6)

## While the rare method deals with the least prolific clonotypes
imm_rare <- repClonality(immdata$data, .method = "rare")
plotfile <- paste0(plotdir, '/rare_clonal_proportion.pdf')
p1 <- vis(imm_rare)
ggsave(plotfile, plot=p1, width=6, height=6)

## the homeo method assesses the clonal space homeostasis, i.e., the proportion of the repertoire occupied by the clones of a given size
imm_hom <- repClonality(immdata$data, .method = "homeo", .clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1))
plotfile <- paste0(plotdir, '/relative_abundance.pdf')
p1 <- vis(imm_hom)
ggsave(plotfile, plot=p1, width=6, height=6)

}

##===========
# diversity measures
##===========
# Chao1 diversity measure
# Chao1 estimator is a nonparameteric asymptotic estimator of species richness (number of species in a population).
div_chao <- repDiversity(immdata$data, "chao1")
plotfile <- paste0(plotdir, '/diversity_measure_chao1.pdf')
p1 <- vis(div_chao)
ggsave(plotfile, plot=p1, width=6, height=6)

# Hill numbers
# Hill numbers are a mathematically unified family of diversity indices (differing only by an exponent q).
div_hill <- repDiversity(immdata$data, "hill")
plotfile <- paste0(plotdir, '/diversity_measure_hill.pdf')
p1 <- vis(div_hill)
ggsave(plotfile, plot=p1, width=6, height=6)

# D50 - number of clonotypes occupying 50% of repetoires
div_d50 <- repDiversity(immdata$data, "d50")
plotfile <- paste0(plotdir, '/diversity_measure_d50.pdf')
p1 <- vis(div_d50)
ggsave(plotfile, plot=p1, width=6, height=6)

# Ecological diversity measure
# True diversity, or the effective number of types, refers to the number of equally-abundant types needed for the average proportional abundance of the types to equal that observed in the dataset of interest where all types may not be equally abundant.
div_div <- repDiversity(immdata$data, "div")
plotfile <- paste0(plotdir, '/diversity_measure_div.pdf')
p1 <- vis(div_div)
ggsave(plotfile, plot=p1, width=6, height=6)

# Gini-Simpson index is the probability of interspecific encounter, i.e., probability that two entities represent different types.
div_gini_simp <- repDiversity(immdata$data, "gini.simp")
plotfile <- paste0(plotdir, '/diversity_measure_gini_symposon_index.pdf')
p1 <- vis(div_gini_simp)
ggsave(plotfile, plot=p1, width=6, height=6)

# Inverse Simpson index is the effective number of types that is obtained when the weighted arithmetic mean is used to quantify average proportional abundance of types in the dataset of interest.
div_inv_simp <- repDiversity(immdata$data, "inv.simp")
plotfile <- paste0(plotdir, '/diversity_measure_inverse_symposon_index.pdf')
p1 <- vis(div_inv_simp)
ggsave(plotfile, plot=p1, width=6, height=6)

if (0) {
	# The Gini coefficient measures the inequality among values of a frequency distribution (for example levels of income). A Gini coefficient of zero expresses perfect equality, where all values are the same (for example, where everyone has the same income). A Gini coefficient of one (or 100 percents ) expresses maximal inequality among values (for example where only one person has all the income).
	div_gini <- repDiversity(immdata$data, "gini")
	plotfile <- paste0(plotdir, '/diversity_measure_gini.pdf')
	p1 <- vis(div_gini)
	ggsave(plotfile, plot=p1, width=6, height=6)
}

# Rarefaction is a technique to assess species richness from the results of sampling through extrapolation
imm_raref <- repDiversity(immdata$data, "raref", .verbose = F)
plotfile <- paste0(plotdir, '/diversity_measure_rarefaction.pdf')
p1 <- vis(imm_raref)
ggsave(plotfile, plot=p1, width=6, height=6)

plotfile <- paste0(plotdir, '/diversity_measure_rarefaction_log.pdf')
p1 <- vis(imm_raref, .log = TRUE)
ggsave(plotfile, plot=p1, width=6, height=6)


##================
## Gene usage computation
## kept for future implementation - requires multi sample data
##================


##================
## Spectratyping
## useful way to represent distributions of genes per sequence length. 
## Parameter .quant controls the quantity that used to compute proportions of genes - either by clonotype (id) or by number of clones per clonotype (count). 
## Parameter .col controls which column to choose, e.g., “nt” for lengths of CDR3 nucleotide sequences only (without grouping by gene segments), “aa+v” for lengths of CDR3 amino acid sequences (grouped by V gene segments).
##================
plotfile <- paste0(plotdir, '/Spectratyping_CDR3_nt.pdf')
p1 <- vis(spectratype(immdata$data[[1]], .quant = "id", .col = "nt"))
ggsave(plotfile, plot=p1, width=6, height=6)

plotfile <- paste0(plotdir, '/Spectratyping_CDR3_aa.pdf')
p2 <- vis(spectratype(immdata$data[[1]], .quant = "count", .col = "aa+v"))
ggsave(plotfile, plot=p2, width=6, height=6)


##=================
## Tracking of clonotypes
## kept for future implementation - requires multi sample data
##=================








##=================
## Kmer statistics computation
##=================
kmers_len3_coding <- getKmers(immdata$data[[1]], 3)
plotfile <- paste0(plotdir, '/kmers_len3_coding.pdf')
p1 <- vis(kmers_len3_coding)
ggsave(plotfile, plot=p1, width=6, height=6)
kp <- kmer_profile(kmers_len3_coding, "self")
p1 <- vis(kp)
p2 <- vis(kp, .plot = "seq")
plotfile <- paste0(plotdir, '/kmers_len3_coding_motif_text_logo.pdf')
ggsave(plotfile, plot=p1, width=6, height=6)
plotfile <- paste0(plotdir, '/kmers_len3_coding_motif_sequence_logo.pdf')
ggsave(plotfile, plot=p2, width=6, height=6)

kmers_len3_all <- getKmers(immdata$data[[1]], 3, .coding = F)
plotfile <- paste0(plotdir, '/kmers_len3_all.pdf')
p1 <- vis(kmers_len3_all)
ggsave(plotfile, plot=p1, width=6, height=6)
kp <- kmer_profile(kmers_len3_all, "self")
p1 <- vis(kp)
p2 <- vis(kp, .plot = "seq")
plotfile <- paste0(plotdir, '/kmers_len3_all_motif_text_logo.pdf')
ggsave(plotfile, plot=p1, width=6, height=6)
plotfile <- paste0(plotdir, '/kmers_len3_all_motif_sequence_logo.pdf')
ggsave(plotfile, plot=p2, width=6, height=6)

kmers_len5_coding <- getKmers(immdata$data[[1]], 5)
plotfile <- paste0(plotdir, '/kmers_len5_coding.pdf')
p1 <- vis(kmers_len5_coding)
ggsave(plotfile, plot=p1, width=6, height=6)
kp <- kmer_profile(kmers_len5_coding, "self")
p1 <- vis(kp)
p2 <- vis(kp, .plot = "seq")
plotfile <- paste0(plotdir, '/kmers_len5_coding_motif_text_logo.pdf')
ggsave(plotfile, plot=p1, width=6, height=6)
plotfile <- paste0(plotdir, '/kmers_len5_coding_motif_sequence_logo.pdf')
ggsave(plotfile, plot=p2, width=6, height=6)

kmers_len5_all <- getKmers(immdata$data[[1]], 5, .coding = F)
plotfile <- paste0(plotdir, '/kmers_len5_all.pdf')
p1 <- vis(kmers_len5_all)
ggsave(plotfile, plot=p1, width=6, height=6)
kp <- kmer_profile(kmers_len5_all, "self")
p1 <- vis(kp)
p2 <- vis(kp, .plot = "seq")
plotfile <- paste0(plotdir, '/kmers_len5_all_motif_text_logo.pdf')
ggsave(plotfile, plot=p1, width=6, height=6)
plotfile <- paste0(plotdir, '/kmers_len5_all_motif_sequence_logo.pdf')
ggsave(plotfile, plot=p2, width=6, height=6)


