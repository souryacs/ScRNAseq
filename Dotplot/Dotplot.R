library(Seurat)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
SeuratObjFile <- args[1]
GeneListFile <- args[2]
OutPlotFile <- args[3]

sc.obj <- readRDS(SeuratObjFile)

GeneListData <- read.table(GeneListFile, header=F, sep="\t", stringsAsFactors=F)
GeneList <- GeneListData[,1]
cat(sprintf("\n Input gene list - number of genes : %s  gene list : %s ", length(GeneList), paste(GeneList, collapse=" ")))

p <- DotPlot(object = sc.obj, features = GeneList, assay="RNA") + scale_colour_gradient2(low="royalblue3", mid="khaki1", high="tomato3") + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
ggsave(filename=OutPlotFile, plot=p, width=12, height=6)
