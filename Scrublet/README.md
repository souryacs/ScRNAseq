Doublet detection of scRNA-seq data using Scrublet
==================================================

The Scrublet package (https://github.com/swolock/scrublet) is an effective and widely-used tool to detect doublets from scRNA-seq data.


	- The script *Run_Scrublet.py* runs Scrublet on input 10X Cellranger output data. It reqiures three inputs:

		1. features.tsv (gene filename)		
		2. matrix.mtx (matrix filename)
		3. Outfile (Output file name to store the results)

		*Note*: These files must be unzipped.

	- Output files:

		- The *Outfile* contains two columns: *doublet_score* and *predicted_doublet*
		- if *predicted_doublet* is True, the corresponding cell is a doublet and needs to be discarded before postprocessing

Queries
=======

For any query, please e-mail:
Sourya Bhattacharyya

sourya@lji.org

