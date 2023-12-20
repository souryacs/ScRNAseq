#!/share/apps/R/3.4.3/bin/Rscript

# https://github.com/mojaveazure/seurat-disk

library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)

## user needs to call the script with two arguments
## 1: Seurat object (.rds) file
## 2. Output .h5ad file

args <- commandArgs(trailingOnly = TRUE)
InpObjFile <- args[1]
outfile <- args[2]

sc.obj <- readRDS(InpObjFile)

outfile_temp <- gsub("h5ad", "h5Seurat", outfile)
SaveH5Seurat(sc.obj, filename = outfile_temp, overwrite=T)
Convert(outfile_temp, dest = "h5ad", overwrite = TRUE)

