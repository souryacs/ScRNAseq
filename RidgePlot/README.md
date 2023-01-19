Seurat Ridge plot
======================

	For a given list of genes or proteins, generates the ridge plot (https://satijalab.org/seurat/reference/ridgeplot).

	Rscript RidgesPlot.r *SeuratObjFile* *GeneListFile* *OutPlotFile*

	where,

		*SeuratObjFile*: Seurat object file (in .rds format)

		*GeneListFile*: File containing the target genes / antobodies for plotting (one gene per row - assuming no header)

		*OutPlotFile*: Output plot file.


Installation Requirements
============================

	The following R packages need to be installed:

		ggplot2, Seurat, ggridges, patchwork


Example
==========

	For example, let the genelist file contain the following entries (antibodies from CITE-seq data):

				CD11c-TotalSeqC
				CD21-TotalSeqC
				CD27-TotalSeqC 
				CD73-TotalSeqC
				CD80-TotalSeqC
				CXCR3-TotalSeqC
				IgM-TotalSeqC
				IgG-TotalSeqC
				
	The output plot file will contain Ridge plots (https://satijalab.org/seurat/reference/ridgeplot) for these antibodies. 


Queries
=======

For any query, please e-mail:
Sourya Bhattacharyya

sourya@lji.org


