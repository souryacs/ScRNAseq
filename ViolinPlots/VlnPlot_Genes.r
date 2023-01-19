library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
SeuratObjFile <- args[1]
OutDir <- args[2]

sc.obj <- readRDS(SeuratObjFile)

pdf(paste0(OutDir, "/Target_Genes_VlnPlot.pdf"))
print(VlnPlot(sc.obj, features = c("BACH2"), pt.size = 0, idents=c('0','1','2','3','4','5','6','7','8')))
dev.off()

