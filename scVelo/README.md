RNA velocity analysis of scRNA-seq data (scVelo)
==================================================

scVelo (https://scvelo.readthedocs.io/en/stable/) is the most widely used package to perform RNA velocity analysis from scRNA-seq data.

	- *Convert_SeuratObj_h5ad.R*: R Script to convert the Seurat object (.rds) into .h5ad format, which is required for scVelo.

		- User needs to provide the input Seurat object file and the output file name to store the .h5ad formatted file.

	- *Run_velocyto_script.sh*: Script to run Velocyto to generate the .loom files from the Cellranger output .bam file. This .loom file is an input to scVelo.

		- User needs to check the parameters section of this script and edit the configuration parameters.

	- *scVelo.py*: Script to run scVelo

		- User needs to check the parameter section and provide three inputs: 
			- 1) .h5ad converted Seurat object file, 
			- 2) .loom file from Velocyto, 
			- 3) target gene list (one file with first column containing this list) whose trajectory need to be examined.

		- Output plots will be stored under the directory *Plots*.

		- Check the scVelo manual for the details of the output files.

Queries
=======

For any query, please e-mail:
Sourya Bhattacharyya

sourya@lji.org

