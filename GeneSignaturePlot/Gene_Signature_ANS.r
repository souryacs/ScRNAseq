##============
## Gene Signature plot using the package ANS (adjusted neighborhood scoring)
## Paper: ANS: Adjusted Neighborhood Scoring to improve assessment of gene signatures in single-cell RNA-seq data
## Laure Ciernik et al.
## https://www.biorxiv.org/content/10.1101/2023.09.20.558114v1
##============

library(ggplot2)
library(Seurat)
library(ggplotify)

##================
## parameters
##================

## important - download this source code of the ANS package from GitHub
## https://github.com/lciernik/ANS_signature_scoring
## and place in the current directory
source("adjusted_neighborhood_scoring.R")

## input files - user can edit according to their custom files
SeuratObjFile <- "Custom_Seurat_object.rds"
GeneListFile <- "Sample_Genelist_file.txt"	## sample gene list file - first column should contain the list of genes
OutDir <- getwd()	## assuming current working directory is the output directory - user can edit

system(paste("mkdir -p", OutDir))
sink(paste0(OutDir, '/out.log'))


##================
## Here we execute both ANS as well as the gene signature from Seurat's routine
##================

sc.obj <- readRDS(SeuratObjFile)

GeneListData <- read.table(GeneListFile, header=F, sep="\t", stringsAsFactors=F)
GeneList <- GeneListData[,1]
cat(sprintf("\n Input number of genes : %s ", length(GeneList)))

## create target gene list if some genes are missing in the Seurat object
Target_Genelist <- c()
for (i in 1:length(GeneList)) {
	if (GeneList[i] %in% rownames(sc.obj)) {
		Target_Genelist <- c(Target_Genelist, GeneList[i])
	}
}
cat(sprintf("\n Filtered number of genes : %s ", length(Target_Genelist)))

# #===============================
# # using Seurat's function "AddModuleScore" with default parameters
# #===============================
# sc.obj <- AddModuleScore(object = sc.obj, features = list(Target_Genelist), ctrl = 30, name = "ctrl_30_Seurat")

##===============================
## using Seurat's function "AddModuleScore" with number of control cells = 5
##===============================
sc.obj <- AddModuleScore(object = sc.obj, features = list(Target_Genelist), ctrl = 5, name = "ctrl_5_Seurat")

# ##===============================
# ## using ANS function with number of control cells = 100
# ##===============================
# sc.obj <- AdjustedNeighborhoodScoring(object = sc.obj, features = list(Target_Genelist), ctrl = 30, name = "ctrl_30_ANS")

##===============================
## using ANS function with number of control cells = 5
##===============================
sc.obj <- AdjustedNeighborhoodScoring(object = sc.obj, features = list(Target_Genelist), ctrl = 5, name = "ctrl_5_ANS")

## write the metadata
write.table(sc.obj@meta.data, paste0(OutDir, '/Complete_metadata.txt'), row.names=T, col.names=T, sep="\t", quote=F, append=F)

## plot the cells according to gene signature score
# color palette 
# https://www.r-bloggers.com/2016/07/creating-color-palettes-in-r/
colvec1 = rainbow(50, start=rgb2hsv(col2rgb('brown'))[1], end=rgb2hsv(col2rgb('orange'))[1])
colvec2 = rainbow(50, start=rgb2hsv(col2rgb('yellow'))[1], end=rgb2hsv(col2rgb('green'))[1])
colvec3 = rainbow(50, start=rgb2hsv(col2rgb('cyan'))[1], end=rgb2hsv(col2rgb('darkblue'))[1])
cols = rev(c(colvec1, colvec2, colvec3))
mypalette <- colorRampPalette(cols)(255)

# plotfile <- paste0(OutDir, '/Gene_Signature_ctrl_30_Seurat.pdf')
# pdf(plotfile, width=4, height=4)
# print(FeaturePlot(object = sc.obj, pt.size=0.2, features = "ctrl_30_Seurat1") + scale_colour_gradientn(colours = mypalette))
# dev.off()

plotfile <- paste0(OutDir, '/Gene_Signature_ctrl_5_Seurat.pdf')
pdf(plotfile, width=4, height=4)
print(FeaturePlot(object = sc.obj, pt.size=0.2, features = "ctrl_5_Seurat1") + scale_colour_gradientn(colours = mypalette))
dev.off()

# plotfile <- paste0(OutDir, '/Gene_Signature_ctrl_30_ANS.pdf')
# pdf(plotfile, width=4, height=4)
# print(FeaturePlot(object = sc.obj, pt.size=0.2, features = "ctrl_30_ANS1") + scale_colour_gradientn(colours = mypalette))
# dev.off()

plotfile <- paste0(OutDir, '/Gene_Signature_ctrl_5_ANS.pdf')
pdf(plotfile, width=4, height=4)
print(FeaturePlot(object = sc.obj, pt.size=0.2, features = "ctrl_5_ANS1") + scale_colour_gradientn(colours = mypalette))
dev.off()




