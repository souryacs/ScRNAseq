library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
SeuratObjFile <- args[1]
GeneListFile <- args[2]
OutFile <- args[3]

# load the Seurat object
sc.obj <- readRDS(SeuratObjFile)

## get number of clusters
nclust <- length(unique(as.vector(Idents(object = sc.obj))))
cat(sprintf("\n Number of clusters : %s ", nclust))

# read input gene list file
GeneListData <- read.table(GeneListFile, header=F, sep="\t", stringsAsFactors=F) 

List_of_Genes <- GeneListData[,1]
cat(sprintf("\n Input gene list - number of genes : %s  gene list : %s ", length(List_of_Genes), paste(List_of_Genes, collapse=" ")))

# get average gene expression for each cluster and for each of these sample categories
# output: gene X cluster matrix
# sample categories are stored in the ident "SampleName"
AvgExpr_Category_Default <- AverageExpression(object = sc.obj, assays = "RNA", return.seurat = F, add.ident = "SampleName", slot="data")
write.table(AvgExpr_Category_Default, OutFile, row.names=T, col.names=T, sep="\t", quote=F, append=F)


