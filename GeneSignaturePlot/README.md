Seurat Gene Signature
======================

For a given gene list, generates the gene signature plot.

We provide two different scripts for plotting gene signatures.

1. Gene signature score computation, according to Holla et al. Science Advances 7(22), 2021.

	- Theory: 

		- First, gene expression data (Seurat object) is scaled (z-score normalization for individual genes) to remove bias towards highly expressed genes.
		- A cell-specific signature score is computed by first sorting the normalized scaled gene expression values for each cell followed by summing up the indices (ranks) of the signature genes.
		- These signature scores are projected in transcriptomic-based UMAP.

	- Execution:

		Rscript GeneSignaturePlot.r *SeuratObjFile* *GeneListFile* *OutDir*

		where,

			*SeuratObjFile*: Seurat object file (in .rds format)

			*GeneListFile*: File containing the target genes for plotting (one gene per row - assuming no header)

			*OutDir*: Output directory to contain the output files and plots.

2. Gene signature score using Seurat's routine and the recently proposed ANS package (https://www.biorxiv.org/content/10.1101/2023.09.20.558114v1)

	- Theory:

		- Here signature scores are computed by first distributing individual genes into *gene-bins* subject to their expression and then selecting a set of control genes from the specific *gene-bin* containing that gene.
		- Seurat and ANS differ by the selection of control genes.

	- Execution:

		First check the *Gene_Signature_ANS.r* script and edit the parameters section. Then execute the R script.

Installation Requirements
============================

	The following R packages need to be installed:

		ggplot2, Seurat, ggplotify

	To execute the *Gene_Signature_ANS.r* script, user also needs to download the R script *adjusted_neighborhood_scoring.R* from the GitHub repository <https://github.com/lciernik/ANS_signature_scoring> and place it in the current folder.


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

