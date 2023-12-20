
Various scripts for processing single cell RNA-seq data
-----------------------------------------------------

Devloped by : Sourya Bhattacharyya

Supervisors: Dr. Ferhat Ay and Dr. Pandurangan Vijayanand

La Jolla Institute for Immunology

La Jolla, San Diego, CA 92037, USA


Scripts (folders)
==================

	Individual folders contain separate README files for detailed descriptions.

	Pipeline:

		Seurat pipeline to process scRNA-seq data. The input is a folder contain the genes and the count matrices, obtained from CellRanger.

	ViolinPlots:

		Script to plot the violin plots, showing target gene expressions for individual clusters.

	DotPlot:

		Script to plot the dot plots, showing target gene expressions for individual clusters.

	FeaturePlot:

		Script to plot the feature plots, showing target gene expressions for individual clusters (in the complete UMAP).

	MarkerGeneExp:

		Script to dump the marker gene expression per cluster and also per sample, where the sample name information is provided in a metadata.

	GeneSignaturePlot:

		Script to plot the signature scores per cell (check the README for the details) in UMAP.

	RidgePlot:

		Ridge plot using Seurat object, for specific genes and antibodies.
		Here we showed examples of antibodies, using CITE-seq data as an example.

	DiffAbundance_MILO:

		Differential abundance (DA) analysis using Milo. Package: https://marionilab.github.io/miloR/articles/milo_demo.html

	TCR_BCR:

		Various scripts / pipelines for T cell receptor (TCR) and B cell receptor (BCR) data analysis.

		- BCR_Immcantation_Immunarch:

			- Pipeline for processing BCR data and plotting relevant metrics using a combination of Immcantation and Immunarch pipelines. 
			- Details are provided in the corresponding README file.

	BD_Rhapsody_Kallisto:

		Scripts to analyze scRNA-seq data generated from BD Rhapsody protocol, using Kallisto.



Queries
=======

For any query, please e-mail:
Sourya Bhattacharyya

sourya@lji.org
