Single cell RNA-seq analysis using Seurat
=============================================

This script "scRNASeq_Seurat.r" needs to be executed using the command:

	Rscript scRNASeq_Seurat.r


Installation
==============

The following R packages need to be installed:

	plyr; dplyr; Seurat; ggplot2; Matrix; fossil; harmony

	** Note: The package "harmony" is required if the scRNA-seq data requires batch correction.


Execution
================

First, edit the Seurat related parameters in lines 16 - 43.

	** Various Seurat parameters + input data folder (CellRanger data output)

		Specifically, user needs to check the following parameters:

			A) Cluster_Resolution: resolution parameter of Seurat for clustering.

			B) The input data folder associated with the function "Read10X" 

			C) The project name - used in the function "CreateSeuratObject"

			D) Mitochondrial gene name pattern - according to the reference genome.	

			E) Batch correction - User needs to incorporate the field containing the batch information (we assumed the metadata "BatchInfo" contains the batch information - check line 143).

	** Note: If batch correction is performed, the dimension reduction is performed using "harmony". Otherwise, PCA is used for dimension reduction.


Output files
===============

	Important Output files
	-------------------------

		1) Out_UMAP.pdf: Clustering UMAP.

		2) "SeuratObjFile" (input parameter): stores the Seurat object.

		3) ClusterWise_DE_Markers.txt, ClusterWise_DE_Power.txt: Cluster-specific DE genes.

		4) Pairwise_cluster_marker/cluster_*/*_vs_*.csv: Pairwise cluster-specific DE genes.


	Other output files
	--------------------

		1) QC_metrics_Violin_plots.pdf: QC metrics plots

		2) QC_Feature_Relationship.pdf: QC features.

		3) QC_metrics_Violin_plots_after_filtering.pdf: QC metric plots after filtering the cells according to various Seurat parameters.

		4) top_variable_genes.pdf: Plot of HVG genes.

		5) PCA_plot.pdf

		6) PCA_Heatmaps.pdf

		7) GenevsPC_score_plot.pdf: most variable genes and their PC scores

		8) PC_Pvals.pdf: JackStraw and Elbow Plots.


Queries
=======

For any query, please e-mail:
Sourya Bhattacharyya

sourya@lji.org

