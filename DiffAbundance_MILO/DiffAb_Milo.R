#!/usr/bin/env Rscript

##==================
## differential abundance using MILO
## https://rawcdn.githack.com/MarioniLab/miloR/3646391023f600bae00efd9d940b888503d7a536/docs/articles/milo_demo.html
## input: Seurat object with cluster annotation information
##==================

## check the following packages installed or not
library(Seurat)
library(ggplot2)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(patchwork)
library(statmod)

## output directory to contain the results
OutDir <- paste0(getwd(), '/out_Milo')
system(paste("mkdir -p", OutDir))

##==============
## parameters
##==============

## input Seurat object file - user can put custom Seurat object file path here
InpSeuratFile <- 'Inp_Seurat_Object.rds'

## d: number of dimensions used in Milo
NumDim_Milo <- 30

## k: number of neighbors in KNN model
Num_Neigh_Milo <- 25

## prop: proportion of cells to randomly sample to start with
Prop_Cell <- 0.1

## neighborhood computed from MILO
neighborhood_plotfile <- paste0(OutDir, "/neiborhood_size.jpg")

## Out Seurat object with MILO specific differential abundance information
Out_Milo_objfile <- paste0(OutDir, "/Seurat_MILO_obj.RDS")

## cell-specific neighborhood cell count information
cntcell_neigh_textfile <- paste0(OutDir, "/count_cells_neighborhood.txt")

## differential abundance results 
da_sc.obj_milo_outfile <- paste0(OutDir, "/da_results_milo.csv")

## differential abundance plot
da_sc.obj_milo_pvalue_file <- paste0(OutDir, "/DA_pvalue_milo.pdf")

Annotation_Fraction_THR <- 0.7

##==============
## first read the Seurat object file
## if needed, user can also annotate the Seurat object with custom cluster labels, or provide a pre-annotated Seurat object
##==============
sc.obj <- readRDS(InpSeuratFile)

##==============
## apply Milo
##==============
if (file.exists(Out_Milo_objfile) == FALSE) {
	sc.obj_sce <- as.SingleCellExperiment(sc.obj)
	sc.obj_milo <- Milo(sc.obj_sce)

	### d=NumDim_Milo, k=Num_Neigh_Milo
	sc.obj_milo <- buildGraph(sc.obj_milo, k = Num_Neigh_Milo, d = NumDim_Milo)
	sc.obj_milo <- makeNhoods(sc.obj_milo, prop=Prop_Cell, k=Num_Neigh_Milo, d=NumDim_Milo, refined=TRUE)
	jpeg(neighborhood_plotfile, res=600, width=6000, height=4800, pointsize=10)
	plotNhoodSizeHist(sc.obj_milo)
	dev.off()

	## dump the Milo object
	saveRDS(sc.obj_milo, file = Out_Milo_objfile)

} else {
	sc.obj_milo <- readRDS(Out_Milo_objfile)
}

## count the number of cells in each neighborhood and each sample
## sample information is provided in the "orig.ident" field
## of the colData(sc.obj_milo) structure
sc.obj_milo <- countCells(sc.obj_milo, meta.data = data.frame(colData(sc.obj_milo)), sample="orig.ident")
## nXm data frame
## n: number of neighborhoods
## m: number of experimental samples
cntdf <- nhoodCounts(sc.obj_milo)
write.table(cntdf, cntcell_neigh_textfile, row.names=T, col.names=T, sep="\t", quote=F, append=F)

## differential abundance testing
sc.obj_milo_design <- data.frame(colData(sc.obj_milo))[,c("orig.ident", "Disease")]
sc.obj_milo_design <- distinct(sc.obj_milo_design)
rownames(sc.obj_milo_design) <- sc.obj_milo_design$orig.ident
sc.obj_milo <- calcNhoodDistance(sc.obj_milo, d=NumDim_Milo)
saveRDS(sc.obj_milo, Out_Milo_objfile)

## differential abundance results for individual neighborhoods
da_results_sc.obj_milo <- testNhoods(sc.obj_milo, design = ~ Disease, design.df = sc.obj_milo_design)
da_results_sc.obj_milo %>% arrange(SpatialFDR) %>% head()
write.csv(da_results_sc.obj_milo, da_sc.obj_milo_outfile, row.names = TRUE)

pdf(da_sc.obj_milo_pvalue_file)
ggplot(da_results_sc.obj_milo, aes(PValue)) + geom_histogram(bins=50)
ggplot(da_results_sc.obj_milo, aes(logFC, -log10(SpatialFDR))) + geom_point() + geom_hline(yintercept = 1) 
sc.obj_milo <- buildNhoodGraph(sc.obj_milo)

## UMAP + neighborhood
# plotUMAP(sc.obj_milo) + plotNhoodGraphDA(sc.obj_milo, da_results_sc.obj_milo, alpha=0.05) + plot_layout(guides="collect")
## UMAP with cluster annotation (stored in the "ident" field)
plotUMAP(sc.obj_milo, colour_by="ident", text_by = "ident", text_size = 3) + guides(fill="none")
plotNhoodGraphDA(sc.obj_milo, da_results_sc.obj_milo, alpha=0.05)

## differential results 
da_results_sc.obj_milo <- annotateNhoods(sc.obj_milo, da_results_sc.obj_milo, coldata_col="ident")
da_results_sc.obj_milo$ident <- ifelse(da_results_sc.obj_milo$ident_fraction < Annotation_Fraction_THR, "Mixed", da_results_sc.obj_milo$ident)
plotDAbeeswarm(da_results_sc.obj_milo, group.by = "ident")
ggplot(da_results_sc.obj_milo, aes(ident_fraction)) + geom_histogram(bins=50)
dev.off()

## write the results one more time
write.csv(da_results_sc.obj_milo, da_sc.obj_milo_outfile, row.names = TRUE)

saveRDS(sc.obj_milo, Out_Milo_objfile)

