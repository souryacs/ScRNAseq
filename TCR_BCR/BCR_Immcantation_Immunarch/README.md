BCR analysis using Immcantation and Immunarch
==============================================

We present a step by step pipeline to analyze BCR data from 10X.

For Immcantation pipeline, user needs to check the following tutorials:

	- https://immcantation.readthedocs.io/en/stable/docker/pipelines.html
	- https://immcantation.readthedocs.io/en/stable/tutorials/10x_tutorial.html

For Immunarch, user needs to check the following link for package details and tutorials:

	- https://immunarch.com/

Basically, Immcantation aligns and processes the BCR data and constructs lineage trees.
Immunarch uses thoese processed data and plots various metrics and quality scores.
	
User needs to run the following scripts as per the order mentioned below. 

	*Note*: They need to check individual scripts, edit the relevant parameters.

		- 1_change_o_script.sh: Assigns new annotations and infers clonal relationships to 10x Genomics single-cell V(D)J data output by Cell Ranger.
		- 1A_changeo_genotyping_script.sh: infer V segment genotypes using TIgGER.
		- 1B_clonal_threshold_script.sh: Performs automated detection of the clonal assignment threshold
		- 1C_clonal_assign_script.sh: Assigns Ig sequences into clonally related lineages and builds full germline sequences.
		- 2_define_clonal_groups.sh: Define heavy chain clones, Split them in light chains and construct germline V and J sequences
		- 3_build_lineage_tree.sh: Building linage trees
		- 4_plot_lineage_tree.R: R script to plot lineage trees using the output file from the earlier script.
		- Convert_Immcantation_to_Immunarch_data.R: R script to convert the Immcantation pipeline output into Immunarch Compatible format.
		- Immunearch_plots.R: R script to plot various metrics using the input Immunarch Compatible file.

Queries
=======

For any query, please e-mail:
Sourya Bhattacharyya

sourya@lji.org

