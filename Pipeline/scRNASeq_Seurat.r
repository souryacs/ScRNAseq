library(plyr)
library(dplyr)
library(Seurat)
library(ggplot2)
library(Matrix)
library(fossil)
library(harmony)

## output log file
sink("out_log.txt")

##=============
## Parameters - user can edit
##=============

## Output file storing the Seurat object
SeuratObjFile <- "out_Seurat_obj.rds"

min_cells <- 10	#3
min_features <- 500 #200
nFeature_RNA_lower <- 250	#200
nFeature_RNA_upper <- 2500
percent_mt <- 15	#5

# number of dimensions from PCA which will be used
NUM_PCA_DIM <- 30 # 16

## Seurat clustering resolution   
Cluster_Resolution <- 0.4

## Input data folder - from cellranger
CellRangerDataFolder <- 'Input/raw_feature_bc_matrix'

## mitochondrial gene pattern
## if human genome, pattern = "^MT-"
## if mouse genome, pattern = "^mt-"
MitoPattern <- "^MT-"

## project name
ProjName <- "CurrProject"

## whether batch correction is to be performed
Batch_Correction <- FALSE

##=============
## Main code
##=============

inp.data <- Read10X(data.dir = CellRangerDataFolder))

sc.obj <- CreateSeuratObject(counts = inp.data, project = ProjName, min.cells = min_cells, min.features = min_features)
cat(sprintf("\n\n ==>> total number of cells in complete Seurat object : %s ", length(colnames(sc.obj))))

## filter mitochondrial genes
sc.obj[["percent.mt"]] <- PercentageFeatureSet(sc.obj, pattern = MitoPattern)

# dump the QC metrics
cat(sprintf("\n sc.obj@meta.data dimension - nrow : %s ncol : %s colnames : %s ", nrow(sc.obj@meta.data), ncol(sc.obj@meta.data), paste(colnames(sc.obj@meta.data), sep="\t")))
write.table(sc.obj@meta.data, "QC_metrics.txt", row.names=T, col.names=T, sep="\t", quote=F, append=F)

## Visualize QC metrics as a violin plot
pdf("QC_metrics_Violin_plots.pdf", width=8, height=6)
VlnPlot(sc.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

## Visualize the other combinations
pdf("QC_Feature_Relationship.pdf", width=8, height=8)
plot1 <- FeatureScatter(sc.obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sc.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(sc.obj, feature1 = "nFeature_RNA", feature2 = "percent.mt")
CombinePlots(plots = list(plot1, plot2, plot3), ncol=1)
dev.off()

#========= filtering the samples
cat(sprintf("\n\n ************* Seurat object (before filtering) - number of cells : %s  *********** \n\n", length(colnames(sc.obj))))
sc.obj <- subset(sc.obj, subset = nFeature_RNA > nFeature_RNA_lower & nFeature_RNA < nFeature_RNA_upper & percent.mt < percent_mt)
cat(sprintf("\n\n ************* Seurat object (after filtering) - number of cells : %s  *********** \n\n", length(colnames(sc.obj))))

# Visualize QC metrics as a violin plot (after filtering the data)
pdf("QC_metrics_Violin_plots_after_filtering.pdf", width=8, height=6)
VlnPlot(sc.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

##=========== Normalizing the data
sc.obj <- NormalizeData(sc.obj, normalization.method = "LogNormalize", scale.factor = 10000)

##=========== Identification of highly variable features (feature selection)
sc.obj <- FindVariableFeatures(sc.obj, selection.method = "vst", nfeatures = 2000)

##=========== Identify the 30 most highly variable genes
top_genes <- head(VariableFeatures(sc.obj), 30)

##=========== plot variable features with and without labels
pdf("top_variable_genes.pdf", width=10, height=8)
plot1 <- VariableFeaturePlot(sc.obj)
plot2 <- LabelPoints(plot = plot1, points = top_genes, repel = TRUE)
CombinePlots(plots = list(plot1, plot2), ncol=2)
dev.off()

## Write the gene expression result
x <- HVFInfo(object = sc.obj)
write.table(x, file="gene_expression_info.txt", sep="\t", quote=F)

## Scale the data and removing unwanted sources of variation
## Check this paper https://www.nature.com/articles/nbt.3102.pdf
sc.obj <- ScaleData(object = sc.obj, vars.to.regress = c("percent.mt"))

## Perform dimension reduction
sc.obj <- RunPCA(sc.obj, features = VariableFeatures(object = sc.obj))

## Examine and visualize PCA results a few different ways
pdf("PCA_plot.pdf")
DimPlot(object = sc.obj, dims = c(1, 2), reduction = "pca")
dev.off()

# Heatmaps
pdf("PCA_Heatmaps.pdf")
DimHeatmap(object = sc.obj, dims = 1, reduction = "pca", cells = 500, balanced = TRUE)
DimHeatmap(object = sc.obj, dims = 2, reduction = "pca", cells = 500, balanced = TRUE)
DimHeatmap(object = sc.obj, dims = 3, reduction = "pca", cells = 500, balanced = TRUE)
dev.off()

## Plot the most variable genes and their PC scores
pdf("GenevsPC_score_plot.pdf")
VizDimLoadings(sc.obj, dims = 1, reduction = "pca")
VizDimLoadings(sc.obj, dims = 2, reduction = "pca")
VizDimLoadings(sc.obj, dims = 3, reduction = "pca")
dev.off()

## Determine statistically significant principal components
sc.obj <- JackStraw(object = sc.obj, reduction = "pca", dims = NUM_PCA_DIM, num.replicate = 100,  prop.freq = 0.1, verbose = TRUE)
sc.obj <- ScoreJackStraw(object = sc.obj, dims = 1:NUM_PCA_DIM, reduction = "pca")

pdf("PC_Pvals.pdf")
JackStrawPlot(object = sc.obj, dims = 1:NUM_PCA_DIM, reduction = "pca")
ElbowPlot(object = sc.obj)
dev.off()

##======= batch correct using Harmony package
if (Batch_Correction == TRUE) {
  ## suppose the "BatchInfo" field in metadata contains the batch information
  ## perform batch correction using Harmony
  sc.obj <- RunHarmony(sc.obj, "BatchInfo")
}

##=======
## UMAP
##=======
## currently we perform UMAP using PCA reduced space
## if batch correction is performed, use HARMONY
if (Batch_Correction == FALSE) {
  sc.obj <- RunUMAP(sc.obj, reduction = "pca", dims = 1:NUM_PCA_DIM)
} else {
  sc.obj <- RunUMAP(sc.obj, reduction = "harmony", dims = 1:NUM_PCA_DIM)
}

##=======
## Cluster the cells
##=======
# algorithm = 1 means original Louvain algorithm is used
if (1) {
  sc.obj <- FindNeighbors(sc.obj, reduction = "pca", dims = 1:NUM_PCA_DIM)
} else {
  sc.obj <- FindNeighbors(sc.obj, reduction = "harmony", dims = 1:NUM_PCA_DIM)
}

## find the clusters
sc.obj <- FindClusters(sc.obj, resolution = Cluster_Resolution, algorithm = 1)

## save current seurat object into a file
saveRDS(sc.obj, file = SeuratObjFile)

# plot UMAP
pdf("Out_UMAP.pdf")
print(DimPlot(sc.obj, reduction = "umap"))
print(DimPlot(sc.obj, reduction = "umap", label = TRUE, repel = TRUE, label.size = 4))
if (Batch_Correction == TRUE) {
  ## if batch correction by Harmony is done
  ## using the field "BatchInfo"
  print(DimPlot(sc.obj, reduction = "umap", group.by = "BatchInfo"))
}
dev.off()

#====================
## Finding differentially expressed genes (cluster biomarkers)
## find markers for every cluster compared to all remaining cells, report
#====================
DE_Marker_File <- "ClusterWise_DE_Markers.txt"
DE_Power_File <- "ClusterWise_DE_Power.txt"

if ((file.exists(DE_Marker_File) == FALSE) | (file.exists(DE_Power_File) == FALSE)) {
  sc.markers <- FindAllMarkers(object = sc.obj, test.use = "wilcox", only.pos = FALSE, min.pct = 0.2)
  sc.markers.power <- FindAllMarkers(object = sc.obj, test.use = "wilcox", only.pos = FALSE, min.pct = 0.2)
  sc.markers[,"power"] <- NA
  for (i in 1:nrow(sc.markers)) {
    if (length(sc.markers.power[sc.markers.power$gene==sc.markers$gene[i] & sc.markers.power$cluster==sc.markers$cluster[i],]$power) > 0) {
      sc.markers$power[i] <- sc.markers.power[sc.markers.power$gene==sc.markers$gene[i] & sc.markers.power$cluster==sc.markers$cluster[i],]$power
    }
  }
  write.table(sc.markers, file=DE_Marker_File, sep="\t", quote=F)
  write.table(sc.markers.power, file=DE_Power_File, sep="\t", quote=F)

} else {
  sc.markers <- read.table(DE_Marker_File, header=T, sep="\t", stringsAsFactors=F, row.names=1)
  sc.markers.power <- read.table(DE_Power_File, header=T, sep="\t", stringsAsFactors=F, row.names=1)
}

##===============
## Plot top genes per cluster
##===============
n <- length(table(sc.obj@active.ident))
for (i in 1:n) {
  top_DE_gene_filename <- paste0("Clusters_",(i-1),"_TOP_DE_Genes.pdf")
  top_DE_expr_filename <- paste0("Clusters_",(i-1), "_TOP_DE_Expression.pdf")
  if ((file.exists(top_DE_gene_filename) == FALSE) | (file.exists(top_DE_expr_filename) == FALSE)) {  
    pdf(top_DE_gene_filename)
    sc.markers.clus <- sc.markers[sc.markers$cluster == (i - 1),]
    sc.markers.clus <- sc.markers.clus[sc.markers.clus$avg_logFC > 0,]
    if (nrow(sc.markers.clus) > 0) {
      sc.markers.clus <- sc.markers.clus[order(-sc.markers.clus$avg_logFC),]
      print (head(sc.markers.clus))
      genes <- sc.markers.clus$gene
      print (genes)
      for (j in 1:length(genes)) {
        print(FeaturePlot(object = sc.obj, features = genes[j], cols = c("grey", "blue"), reduction = "umap"))
      }
      dev.off()

      pdf(top_DE_expr_filename)
      for (j in 1:length(genes)) {
        print(VlnPlot(object = sc.obj, features = genes[j]))
      }
      dev.off()
    } else {
      cat ("Cluster ",i - 1," dont have any upregulated genes\n")
    }
  }
}

top10 <- sc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC > 0)
pdf("Clusters_TOP_DE_Heatmap.pdf")
DoHeatmap(object = sc.obj, features = top10$gene, label = TRUE, size=10) + theme(text = element_text(size = 1))
dev.off()


#============
# find markers for individual pairs of clusters
# process each cluster, and then process sequentially one from the rest of the clusters
#============
n <- length(table(sc.obj@active.ident)) # number of clusters
all.clusters <- (seq(1,n) - 1)  # clusters are from 0 to n-1
cat(sprintf("\n\n **** Started to process pairwise clusters to find markers \n number of clusters : %s ", length(all.clusters)))

for (cluster in all.clusters) {
  ClusterOutDir <- paste0('Pairwise_cluster_marker/cluster_', cluster)
  system(paste("mkdir -p", ClusterOutDir))  
  rest.clusters <- all.clusters[all.clusters != cluster] 
  cat(sprintf("\n ---- analyzing one cluster --- %s ", cluster))
  for (next.cluster in rest.clusters) {
    cat(sprintf("\n -- processing pair -- cluster : %s next.cluster : %s ", cluster, next.cluster))
    cluster.markers <- FindMarkers(object=sc.obj, ident.1=cluster, ident.2=next.cluster, only.pos=FALSE, min.pct=0.20)
    tmp.file.name <- paste0(ClusterOutDir, '/', cluster, '_vs_', next.cluster, '.csv')
    write.csv(x=cluster.markers, file=tmp.file.name)  
  } # end next.cluster loop
} # end cluster loop


