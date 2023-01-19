library(Seurat)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
SeuratObjFile <- args[1]
GeneListFile <- args[2]
Outplotfile <- args[3]

system(paste("mkdir -p", dirname(Outplotfile)))

# read Seurat object
sc.obj <- readRDS(SeuratObjFile)
DefaultAssay(sc.obj) <- 'RNA'

# read input gene list file
GeneListData <- read.table(GeneListFile, header=F, sep="\t", stringsAsFactors=F) 

GeneList <- GeneListData[,1]
cat(sprintf("\n Input gene list - number of genes : %s  gene list : %s ", length(GeneList), paste(GeneList, collapse=" ")))

pdf(Outplotfile, width = 8, height = 20)
print(FeaturePlot(object = sc.obj, features = intersect(GeneList, rownames(sc.obj)), reduction = "umap", cols = c("grey", "blue"), ncol = 3, label = TRUE, label.size = 1, pt.size = 0.1, repel = TRUE))
dev.off()

