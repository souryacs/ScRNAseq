Seurat Gene Signature
======================

	For a given gene list, generates the gene signature plot.

	Rscript GeneSignaturePlot.r *SeuratObjFile* *GeneListFile* *OutDir*

	where,

		*SeuratObjFile*: Seurat object file (in .rds format)

		*GeneListFile*: File containing the target genes for plotting (one gene per row - assuming no header)

		*OutDir*: Output directory to contain the output files and plots.


Theory
========

	Gene signature score computation is described in Holla et al. Science Advances 7(22), 2021.

	First, gene expression data (Seurat object) is scaled (z-score normalization for individual genes) to remove bias towards highly expressed genes. A cell-specific signature score is computed by first sorting the normalized scaled gene expression values for each cell followed by summing up the indices (ranks) of the signature genes. These signature scores are projected in transcriptomic-based UMAP.



Installation Requirements
============================

	The following R packages need to be installed:

		ggplot2, Seurat, ggplotify


Output Files
============

	Two important output files:

		a) Signature_Score.txt: Cell-specific signature scores.

		b) Gene_Signature.pdf: UMAP plot of signature scores.



Queries
=======

For any query, please e-mail:
Sourya Bhattacharyya

sourya@lji.org

