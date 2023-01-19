library(ggplot2)
library(ggridges)
library(Seurat)
library(patchwork)

sink('out.log')

args <- commandArgs(trailingOnly = TRUE)
SeuratObjFile <- args[1]
GeneListFile <- args[2]
OutPlotFile <- args[3]

GeneListData <- read.table(GeneListFile, header=F, sep="\t", stringsAsFactors=F)
GeneList <- GeneListData[,1]
cat(sprintf("\n Input number of genes : %s ", length(GeneList)))

## here we assume that input Seurat object corresponds to CITE-seq data
## so it has also an antibody matrix
sc.obj <- readRDS(SeuratObjFile)
## Get antibody matrix
abmat <- as.matrix(sc.obj$`Antibody Capture`)
cat(sprintf("\n\n ===>>> Dimension of abmat : %s X %s ", nrow(abmat), ncol(abmat)))

##===== insert the cluster information
sc.obj$clusterID <- sc.obj@meta.data$seurat_clusters
nclust <- length(unique(sc.obj$clusterID))
cat(sprintf("\n\n ** nclust : %s ", nclust))

##==== check the antibody capture data and only select those cells which are present in the "sc.obj" gene expression data
m <- match(colnames(abmat), colnames(sc.obj))
idx_s <- which(!is.na(m))
cat(sprintf("\n Number of entries in colnames(sc.obj) : %s ", length(colnames(sc.obj))))
cat(sprintf("\n Number of entries in colnames(abmat) : %s ", length(colnames(abmat))))
cat(sprintf("\n Number of matching entries : %s ", length(idx_s)))

##=== extract the subset of antibody capture matrix
abmat_sub <- abmat[, idx_s]

##===== create antibody Capture assay in the existing Seurat object
##===== integrated data from all four samples
sc.obj[["ADT"]] <- CreateAssayObject(counts = as.sparse(abmat_sub))

##==== normalize and scale the antibody Capture object
sc.obj <- NormalizeData(sc.obj, assay = "ADT", normalization.method = "CLR")
sc.obj <- ScaleData(sc.obj, assay = "ADT")

##==== visualize the Ridge plots for the integrated data
pdf(OutPlotFile)
print(RidgePlot(sc.obj, features = GeneList, group.by = 'clusterID', ncol = 2))
dev.off()

