##=================
## script implementing RNA velocity on its test data

## before running this script
## check this link: https://github.com/basilkhuder/Seurat-to-RNA-Velocity/issues/8

## python3 -m pip uninstall numba
## python3 -m pip install numba==0.52.0
##=================

import os
import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import matplotlib
import matplotlib.pyplot as plt

## to download data and not worry about
## ssl.SSLCertVerificationError: [SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed: unable to get local issuer certificate
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

scv.logging.print_version()

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization

## create output directory, and make it the current directory
outdir = os.getcwd()
plotdir = outdir + '/Plots'
os.makedirs(plotdir, exist_ok=True)

##================
## parameters
##================
## .loom file containing gene expression
## user needs to edit according to the output of Velocyto
GeneExprFile_1 = '/path/to/velocyto_out/currSample.loom'

## .h5ad file containing clustering information (Seurat object)
## generated from the script Convert_SeuratObj_h5ad.R
ClusterFile = "/path/to/SeuratObj.h5ad"

## target genes whose characteristics is to be plotted
## user needs to edit according to the the target gene list
## the first column of this file contains the list of genes
TargetGenelistFile = "/path/to/Genes.txt"

##================
## execution
##================

## read this gene list
genedf = pd.read_table(TargetGenelistFile, header=None)
gene_list = genedf[0]	# get the first column - gene names

## this fiele stores the processed (velocity estimated) gene expression data
processed_file = outdir + "/scvelo_processed.h5ad"
if not os.path.exists(processed_file):
# if 1:

	## read cluster information from Seurat object
	clusterdata = scv.read(ClusterFile)

	## read input gene expression
	geneexprdata = scv.read(GeneExprFile_1)
	geneexprdata.var_names_make_unique()

	## merge with clustering information
	adata = scv.utils.merge(geneexprdata, clusterdata)

	## this statement is to be checked
	## according to the input data
	## copy the cluster information
	adata.obs['clusters'] = adata.obs['seurat_clusters']

	## required due to this error
	## https://forums.fast.ai/t/attributeerror-can-only-use-cat-accessor-with-a-category-dtype/53198/2
	adata.obs['clusters'] = adata.obs['clusters'].astype('category')

	## scVelo is based on adata, 
	## an object that stores a data matrix adata.X, 
	## annotation of observations adata.obs, 
	## variables adata.var, 
	## and unstructured annotations adata.uns. 
	## Names of observations and variables can be accessed via adata.obs_names and adata.var_names, respectively. 
	## AnnData objects can be sliced like dataframes, for example, adata_subset = adata[:, list_of_gene_names]
	## for details: https://anndata.readthedocs.io/en/latest/

	## plot the proportion as a file
	scv.pl.proportions(adata)
	outfile = plotdir + '/1_Plot_Proportions.png'
	plt.savefig(outfile)

	##===============
	## Preprocess the Data
	##===============
	## gene selection by detection (with a minimum number of counts) and high variability (dispersion), 
	## normalizing every cell by its total size and logarithmizing X. 
	## Filtering and normalization is applied in the same vein to spliced/unspliced counts and X. 
	## Logarithmizing is only applied to X. 
	## If X is already preprocessed from former analysis, it will not be touched.

	## All of this is summarized in a single function scv.pp.filter_and_normalize, which essentially runs the following:
	## 1. scv.pp.filter_genes(adata, min_shared_counts=20)
	## 2. scv.pp.normalize_per_cell(adata)
	## 3. scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
	## 4. scv.pp.log1p(adata)

	if 1:
		## original code
		# scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
		## modified code
		scv.pp.filter_and_normalize(adata, min_shared_counts=0) #, n_top_genes=10000)
		print("\n\n **** after filtering - No of genes in adata : " + str(len(adata.var.index)))

	# ##===============================
	# ## the following code is to be used 
	# ## only when we do not have reference clustering information (from Seurat / Signac)
	# ## currently tthe code is not fully functional
	# ##===============================
	# if not os.path.exists(ClusterFile):

	# 	scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

	# 	## Further preprocessing (such as batch effect correction) may be used to remove unwanted sources of variability
	# 	## Note: we did not perform batch correction

	# 	##================
	# 	## clustering and embedding
	# 	##================
	# 	## Computing the neighborhood graph
	# 	sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

	# 	## Embedding the neighborhood graph
	# 	sc.tl.umap(adata)

	# 	## Clustering the neighborhood graph (Leiden algorithm)
	# 	## similar to Seurat 
	# 	## implemented in Scanpy
	# 	sc.tl.leiden(adata)

	# 	## plot using scanpy umap
	# 	sc.pl.umap(adata, color=['leiden', 'CST3', 'NKG7'])
	# 	outfile = plotdir + '/1A_cluster_umap.png'
	# 	plt.savefig(outfile)	

	# 	## plot using scanpy umap - normalized count
	# 	sc.pl.umap(adata, color=['leiden', 'CST3', 'NKG7'], use_raw=False)
	# 	outfile = plotdir + '/1B_cluster_umap_normalized.png'
	# 	plt.savefig(outfile)

	##================
	## Estimate RNA velocity
	##================

	## Velocities are vectors in gene expression space and represent the direction and speed of movement of the individual cells.
	## velocities are obtained by modeling transcriptional dynamics of splicing kinetics, either stochastically (default) or deterministically (by setting mode='deterministic'). 
	## For each gene, a steady-state-ratio of pre-mature (unspliced) and mature (spliced) mRNA counts is fitted, which constitutes a constant transcriptional state. 
	## Velocities are then obtained as residuals from this ratio. 
	## Positive velocity indicates that a gene is up-regulated, which occurs for cells that show higher abundance of unspliced mRNA for that gene than expected in steady state. 
	## Conversely, negative velocity indicates that a gene is down-regulated.

	## here we'll employ velocity analysis by stochastically (default) 
	scv.tl.velocity(adata)

	## The computed velocities are stored in adata.layers just like the count matrices.

	##================
	## estimate cell to cell transitions using velocity
	##================

	## The combination of velocities across genes can then be used to estimate the future state of an individual cell. 
	## In order to project the velocities into a lower-dimensional embedding, transition probabilities of cell-to-cell transitions are estimated. 
	## That is, for each velocity vector we find the likely cell transitions that are accordance with that direction. 
	## The transition probabilities are computed using cosine correlation between the potential cell-to-cell transitions and the velocity vector, and are stored in a matrix denoted as velocity graph. 
	## The resulting velocity graph has dimension n_obs X n_obs 
	## and summarizes the possible cell state changes that are well explained through the velocity vectors 
	## (for runtime speedup it can also be computed on reduced PCA space by setting approx=True).

	scv.tl.velocity_graph(adata)

	## For a variety of applications, the velocity graph can be converted to a transition matrix by applying a Gaussian kernel to transform the cosine correlations into actual transition probabilities. 
	## You can access the Markov transition matrix via scv.utils.get_transition_matrix.

	##================	
	# Save the result for use later on
	adata.write(processed_file)

else:
	# read the existing results
	adata = sc.read_h5ad(processed_file)

##================
## Project the velocities
##================

## Note: this test data has an already pre-computed UMAP embedding, and annotated clusters. 
## When applying to your own data, these can be obtained with scv.tl.umap and scv.tl.louvain
## For more details, see the scanpy tutorial (https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html)
## plotting functions are defaulted to using basis='umap' and color='clusters', which you can set accordingly.
## The velocity vector field displayed as streamlines yields fine-grained insights into the developmental processes. It accurately delineates the cycling population of ductal cells and endocrine progenitors. Further, it illuminates cell states of lineage commitment, cell-cycle exit, and endocrine cell differentiation.

scv.pl.velocity_embedding_stream(adata, basis='umap', color='clusters')
outfile = plotdir + '/2_velocity_embedding_umap.png'
plt.savefig(outfile)

## The most fine-grained resolution of the velocity vector field we get at single-cell level, with each arrow showing the direction and speed of movement of an individual cell.
scv.pl.velocity_embedding(adata, color='clusters', arrow_length=3, arrow_size=2, dpi=120)
outfile = plotdir + '/2A_velocity_embedding_umap_1.png'
plt.savefig(outfile)


##================
## Interprete the velocities
##================

## check the gif
## https://user-images.githubusercontent.com/31883718/80227452-eb822480-864d-11ea-9399-56886c5e2785.gif

## Gene activity is orchestrated by transcriptional regulation. 
## Transcriptional induction for a particular gene results in an increase of (newly transcribed) precursor unspliced mRNAs while, conversely, repression or absence of transcription results in a decrease of unspliced mRNAs. 
## Spliced mRNAs is produced from unspliced mRNA and follows the same trend with a time lag. 
## Time is a hidden/latent variable. 
## Thus, the dynamics needs to be inferred from what is actually measured: spliced and unspliced mRNAs as displayed in the phase portrait.

## examine the phase portraits of some marker genes, 
## visualized with scv.pl.velocity(adata, gene_names) 
## or scv.pl.scatter(adata, gene_names).

## The black line corresponds to the estimated ‘steady-state’ ratio, i.e. the ratio of unspliced to spliced mRNA abundance which is in a constant transcriptional state. 
## RNA velocity for a particular gene is determined as the residual, i.e. how much an observation deviates from that steady-state line. 
## Positive velocity indicates that a gene is up-regulated, which occurs for cells that show higher abundance of unspliced mRNA for that gene than expected in steady state. 
## Conversely, negative velocity indicates that a gene is down-regulated.

scv.pl.velocity(adata, gene_list, ncols=2)
outfile = plotdir + '/3_phase_portrait_velocity_marker_genes.png'
plt.savefig(outfile, bbox_inches="tight")

for i in range(len(gene_list)):
	if gene_list[i] in adata.var.index:
		# scv.pl.velocity(adata, gene_list, basis='umap', ncols=1)      #2)
		scv.pl.velocity(adata, gene_list[i], ncols=1)
		outfile = plotdir + '/3_phase_portrait_velocity_marker_genes_' + gene_list[i] + '.png'
		plt.savefig(outfile, bbox_inches="tight")


## scatterplot depicting the cell transition
## w.r.t cluster and velocity
scv.pl.scatter(adata, gene_list, color=['clusters', 'velocity'])        #, add_outline='Ngn3 high EP, Pre-endocrine, Beta')
outfile = plotdir + '/4_scatter_velocity_marker_genes.png'
plt.savefig(outfile, bbox_inches="tight")

for i in range(len(gene_list)):
	if gene_list[i] in adata.var.index:
		scv.pl.scatter(adata, gene_list[i], color=['clusters', 'velocity'])
		outfile = plotdir + '/4_scatter_velocity_marker_genes_' + gene_list[i] + '.png'
		plt.savefig(outfile, bbox_inches="tight")

##=============== 
## Identify important genes
##===============

## we can test which genes have cluster-specific differential velocity expression, being siginificantly higher/lower compared to the remaining population. 
## The module scv.tl.rank_velocity_genes runs a differential velocity t-test and outpus a gene ranking for each cluster. 
## Thresholds can be set (e.g. min_corr) to restrict the test on a selection of gene candidates.

scv.tl.rank_velocity_genes(adata, groupby='clusters', min_corr=.3)

df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head()

outfile = plotdir + '/5_top_ranked_genes_differential_velocity_expression.csv' 
df.to_csv(outfile)


# # ## now check individual clusters and plot the corresponding top genes specific to that cluster
# # ## Note: we do have cluster information for this data
# # ## otherwise, we need to first define and annotate the clusters
# # ## one cluster has name: Ngn3 high EP
# # ## other cluster has name: Pre-endocrine

# # kwargs = dict(frameon=False, size=10, linewidth=1.5, add_outline='Ngn3 high EP, Pre-endocrine, Beta')

# # scv.pl.scatter(adata, df['Ngn3 high EP'][:5], ylabel='Ngn3 high EP', **kwargs)
# # scv.pl.scatter(adata, df['Pre-endocrine'][:5], ylabel='Pre-endocrine', **kwargs)

# # outfile = plotdir + '/6_top_ranked_genes_various_clusters_differential_velocity_expression.png' 
# # plt.savefig(outfile)


# # ##================ 
# # ## Velocities in cycling progenitors
# # ##================

# # ## The cell cycle detected by RNA velocity, is biologically affirmed by cell cycle scores 
# # ## (standardized scores of mean expression levels of phase marker genes).

# # scv.tl.score_genes_cell_cycle(adata)
# # scv.pl.scatter(adata, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[5, 95])

# # outfile = plotdir + '/7_velocity_cell_cycle_scores.png' 
# # plt.savefig(outfile)

# # ## For the cycling Ductal cells, we may screen through S and G2M phase markers. 
# # ## The previous module also computed a spearmans correlation score, 
# # ## which we can use to rank/sort the phase marker genes to then display their phase portraits.

# # s_genes, g2m_genes = scv.utils.get_phase_marker_genes(adata)
# # s_genes = scv.get_df(adata[:, s_genes], 'spearmans_score', sort_values=True).index
# # g2m_genes = scv.get_df(adata[:, g2m_genes], 'spearmans_score', sort_values=True).index

# # kwargs = dict(frameon=False, ylabel='cell cycle genes')
# # scv.pl.scatter(adata, list(s_genes[:2]) + list(g2m_genes[:3]), **kwargs)

# # outfile = plotdir + '/8_velocity_cell_cycle_scores_top_genes.png' 
# # plt.savefig(outfile)

# # ## Particularly Hells and Top2a are well-suited to explain the vector field in the cycling progenitors. 
# # ## Top2a gets assigned a high velocity shortly before it actually peaks in the G2M phase. 
# # ## There, the negative velocity then perfectly matches the immediately following down-regulation.

# # scv.pl.velocity(adata, ['Hells', 'Top2a'], ncols=2, add_outline=True)

# # outfile = plotdir + '/9_velocity_cell_cycle_scores_selected_genes_outline.png' 
# # plt.savefig(outfile)

##==============
## Speed and coherence
##==============

## Two more useful stats: 
## The speed or rate of differentiation is given by the length of the velocity vector. 
## The coherence of the vector field (i.e., how a velocity vector correlates with its neighboring velocities) provides a measure of confidence.
## These provide insights where cells differentiate at a slower/faster pace, and where the direction is un-/determined.

scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95])

outfile = plotdir + '/10_speed_coherence_velocity.png' 
plt.savefig(outfile)


# # ## On cluster-level, we find that differentiation substantially speeds up after cell cycle exit (Ngn3 low EP), 
# # ## keeping the pace during Beta cell production while slowing down during Alpha cell production.
# # df = adata.obs.groupby('clusters')[keys].mean().T
# # df.style.background_gradient(cmap='coolwarm', axis=1)

# # outfile = plotdir + '/11_cluster_level_velocity_length_confidence.csv' 
# # df.to_csv(outfile)

##============== 
## Velocity graph and pseudotime
##==============

## We can visualize the velocity graph to portray all velocity-inferred cell-to-cell connections/transitions.
scv.pl.velocity_graph(adata, threshold=.1)
outfile = plotdir + '/12_velocity_graph.png' 
plt.savefig(outfile)

## Further, the graph can be used to draw descendents/anscestors coming from a specified cell. 
## Here, a pre-endocrine cell is traced to its potential fate.
x, y = scv.utils.get_cell_transitions(adata, basis='umap', starting_cell=70)
ax = scv.pl.velocity_graph(adata, c='lightgrey', edge_width=.05, show=False)
ax = scv.pl.scatter(adata, x=x, y=y, s=120, c='ascending', cmap='gnuplot', ax=ax)
outfile = plotdir + '/13_velocity_cell_transitions.png' 
plt.savefig(outfile)

## Finally, based on the velocity graph, a velocity pseudotime can be computed. 
## After inferring a distribution over root cells from the graph, 
## it measures the average number of steps it takes to reach a cell after 
## walking along the graph starting from the root cells.
scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot')
outfile = plotdir + '/14_velocity_pseudotime.png' 
plt.savefig(outfile)


# # ##=====================
# # ## PAGA velocity graph
# # ##=====================

# # if 0:
# # 	## PAGA graph abstraction has benchmarked as top-performing method for trajectory inference. 
# # 	## It provides a graph-like map of the data topology with weighted edges corresponding to the connectivity between two clusters. 
# # 	## Here, PAGA is extended by velocity-inferred directionality.

# # 	# this is needed due to a current bug - bugfix is coming soon.
# # 	adata.uns['neighbors']['distances'] = adata.obsp['distances']
# # 	adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']

# # 	scv.tl.paga(adata, groups='clusters')
# # 	df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
# # 	df.style.background_gradient(cmap='Blues').format('{:.2g}')

# # 	outfile = plotdir + '/15_PAGA_transitions_confidence.csv' 
# # 	df.to_csv(outfile)

# # 	## This table can be summarized by a directed graph superimposed onto the UMAP embedding.
# # 	scv.pl.paga(adata, basis='umap', size=50, alpha=.1,
# # 	            min_edge_width=2, node_size_scale=1.5)

# # 	outfile = plotdir + '/16_PAGA_transitions_graph_UMAP.png' 
# # 	plt.savefig(outfile)

