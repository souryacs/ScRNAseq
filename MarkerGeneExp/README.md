Seurat Average Expression per cluster
=======================================

	Dumps the average expression of individual genes per cluster and also per sample type (provided the Seurat object is created by merging multiple input samples / replicates).

		** Here we assume that the field "SampleName" of the input Seurat object contains the sample information. 
		** (check line 24 of the .r code). 
		** User can edit this field name for the custom input Seurat object, or even delete this argument if the input Seurat object is created from one replicate / sample.

	Rscript MarkerGeneExpr.R *SeuratObjFile* *GeneListFile* *OutFile*

	where,

		*SeuratObjFile*: Seurat object file (in .rds format)

		*GeneListFile*: File containing the target genes for plotting (one gene per row - assuming no header)

		*OutFile*: File to contain the average expression per cluster and per sample type.


Queries
=======

For any query, please e-mail:
Sourya Bhattacharyya

sourya@lji.org

