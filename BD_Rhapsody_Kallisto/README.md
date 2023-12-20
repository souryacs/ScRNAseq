Using Kallisto to quantify BD-Rhapsody generated scRNA-seq data
==================================================================

Check the Kallisto package and corresponding manual: https://pachterlab.github.io/kallisto/manual.html

	- Script: Kallisto_Quant.sh
		- Runs the quantification algorithm. 
		- Kallisto can process either single-end or paired-end reads. The default running mode is paired-end and requires an even number of FASTQ files represented as pairs.
		- Important note: only supply one sample at a time to kallisto. The multiple FASTQ (pair) option is for users who have samples that span multiple FASTQ files.
		- Note: User needs to check and edit the parameters.

	 - Script: Kallisto_Bus.sh
	 	- Works with raw FASTQ files for single-cell RNA-Seq datasets. 
		- For each read the cell barcode and UMI information and the equivalence class resulting from pseudoalignment are stored in a BUS file output.bus stored in the output directory directory, along with matrix.ec and transcripts.txt which store information about the equivalence classes and transcript names for downstream processing.
		- Note: User needs to check and edit the parameters.

Queries
=======

For any query, please e-mail:
Sourya Bhattacharyya

sourya@lji.org


