#!/share/apps/R/3.4.3/bin/Rscript

library(ggplot2)
library(Seurat)
library(ggplotify)

args <- commandArgs(trailingOnly = TRUE)
SeuratObjFile <- args[1]
GeneListFile <- args[2]
OutDir <- args[3]

system(paste("mkdir -p", OutDir))

sc.obj <- readRDS(SeuratObjFile)

GeneListData <- read.table(GeneListFile, header=F, sep="\t", stringsAsFactors=F)
GeneList <- GeneListData[,1]
cat(sprintf("\n Input number of genes : %s ", length(GeneList)))

##===== get the normalized + scaled gene expression data
sc.obj <- ScaleData(sc.obj, assay = "RNA")
GeneExpr <- GetAssayData(object=sc.obj, assay="RNA", slot="scale.data")

##==== extract gene expression of only those genes which are provided in the given gene list
idx <- which(rownames(GeneExpr) %in% GeneList)
cat(sprintf("\n Matching number of rows (genes) in sc.obj : %s ", length(idx)))

GeneExpr_sub <- GeneExpr[idx, ]
GeneExpr_RankMat <- matrix(0, nrow=nrow(GeneExpr_sub), ncol=ncol(GeneExpr_sub))
rownames(GeneExpr_RankMat) <- rownames(GeneExpr_sub)
colnames(GeneExpr_RankMat) <- colnames(GeneExpr_sub)
for (i in 1:nrow(GeneExpr_sub)) {
	currgene_exprvec <- GeneExpr_sub[i, ]
	rank_currgene_exprvec <- rank(currgene_exprvec)
	GeneExpr_RankMat[i, ] <- rank_currgene_exprvec
}

## dump the gene expression and the rank information
write.table(as.data.frame(t(GeneExpr_sub)), paste0(OutDir, '/GeneExpr_subset.txt'), row.names=T, col.names=T, sep="\t", quote=F, append=F)
write.table(as.data.frame(t(GeneExpr_RankMat)), paste0(OutDir, '/GeneExpr_Rank.txt'), row.names=T, col.names=T, sep="\t", quote=F, append=F)

## get per cell total rank - signature score
Sig_Score <- as.data.frame(apply(GeneExpr_RankMat, 2, sum))
rownames(Sig_Score) <- colnames(GeneExpr_RankMat)
colnames(Sig_Score) <- 'Signature_score'
write.table(Sig_Score, paste0(OutDir, '/Signature_Score.txt'), row.names=T, col.names=T, sep="\t", quote=F, append=F)

##===== add gene signature score as a metadata to the Seurat object
sc.obj$GeneSignature <- Sig_Score[,1]

## plot the cells according to gene signature score
# color palette - https://www.r-bloggers.com/2016/07/creating-color-palettes-in-r/
colvec1 = rainbow(50, start=rgb2hsv(col2rgb('brown'))[1], end=rgb2hsv(col2rgb('orange'))[1])
colvec2 = rainbow(50, start=rgb2hsv(col2rgb('yellow'))[1], end=rgb2hsv(col2rgb('green'))[1])
colvec3 = rainbow(50, start=rgb2hsv(col2rgb('cyan'))[1], end=rgb2hsv(col2rgb('darkblue'))[1])
cols = rev(c(colvec1, colvec2, colvec3))
mypalette <- colorRampPalette(cols)(255)
plotfile <- paste0(OutDir, '/Gene_Signature.pdf')
pdf(plotfile, width=8, height=6)
FeaturePlot(sc.obj, features = "GeneSignature") + scale_colour_gradientn(colours = mypalette)
dev.off()


