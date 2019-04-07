## Load libraries
library(dplyr)
library(Seurat)

##########################################################################
## Load the raw count dataset
count.data <- Read10X(data.dir = "path/to/countfile")

## Initialize the Seurat object with the raw (non-normalized data).
countDF <- CreateSeuratObject(counts = count.data, project = "Insert_project_name", min.cells = 3, min.features = 200)
countDF

##########################################################################

## QC and selecting cells for further analysis

## The [[ operator can add columns to object metadata. This is a great place to stash QC stats
## The following code detects percentage of reads that map to the mitochondrial genome
## Higher mitochondrial percentage indicates low quality/dying cells
## The pattern to identify mitochondrial genes is set to "MT."
countDF[["percent.mt"]] <- PercentageFeatureSet(object = countDF, pattern = "^MT.")

## Visualize QC metrics as a violin plot
## Use these plots to filter cells
pdf("VlnPlot.pdf")
VlnPlot(object = countDF, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

## FeatureScatter is typically used to visualize feature-feature relationships, but can be used
## for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
pdf("FeaturePlots.pdf")
plot1 <- FeatureScatter(object = countDF, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = countDF, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off

## The following code filters cells that have unique feature counts over 2,500 or less than 200.
## Cells with greater than 5% mitochondrial counts are also filtered
countDF <- subset(x = countDF, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

##########################################################################

## Normalizing the data

## The parameters for normalization are set at default
countDF <- NormalizeData(object = countDF, normalization.method = "LogNormalize", scale.factor = 10000)

##########################################################################

## Identification of highly variable features

## Features that exhibit high cell-to-cell variation in the dataset (i.e. highly expressed in some cells and lowly expressed in others)
## By default, 2,000 features per dataset are returned.
countDF <- FindVariableFeatures(object = countDF, selection.method = "vst", nfeatures = 2000)

## Identify the 10 most highly variable genes
top10 <- head(x = VariableFeatures(object = countDF), 10)

## plot variable features with and without labels
pdf("Variablefeatures.pdf")
plot1 <- VariableFeaturePlot(object = countDF)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()

##########################################################################
## Scaling the data

## This applies a linear transformation, a standard pre-processing step prior to dimensional reduction
## The results of scaling are stored in countDF[["RNA"]]@scale.data
all.genes <- rownames(x = countDF)
countDF <- ScaleData(object = countDF, features = all.genes)

##########################################################################
## Perform linear dimensional reduction
countDF <- RunPCA(object = countDF, features = VariableFeatures(object = countDF))

## Examine and visualize PCA results a few different ways
print(x = countDF[["pca"]], dims = 1:5, nfeatures = 5)

pdf("VizDimLoadings.pdf")
VizDimLoadings(object = countDF, dims = 1:2, reduction = "pca")
dev.off()

pdf("DimPlot.pdf")
DimPlot(object = countDF, reduction = "pca")
dev.off()

## DimHeatmap allows for easy exploration of primary sources of heterogeneity in dataset
## Useful for when trying to decide which PCs to include for further downstream analyses.
## Cells and features are ordered according to their PCA scores
pdf("DimHeatmap.pdf")
DimHeatmap(object = countDF, dims = 1, cells = 500, balanced = TRUE)
dev.off()

pdf("DimHeatmap_multdim.pdf")
DimHeatmap(object = countDF, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

##########################################################################
## Determine the 'dimensionality' of the dataset
## Seurat clusters cells based on PCA scores
## The top principal components represent a robust compression of the dataset.
## First plot the PCs and identify p-values for each PCs
## Choose PCs that are significant in downstream analysis
countDF <- JackStraw(object = countDF, num.replicate = 100)
countDF <- ScoreJackStraw(object = countDF, dims = 1:20)

pdf("JackStrawPlot_PCs.pdf")
JackStrawPlot(object = countDF, dims=1:15)
dev.off()

##########################################################################
## Cluster cells
## Use signficant PCs in "dims" argument. For example here the first 10 PCs are used
countDF <- FindNeighbors(object = countDF, dims = 1:10)
countDF <- FindClusters(object = countDF, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(x = Idents(object = countDF), 5)

##########################################################################
## Run non-linear dimensional reduction (tSNE)
## Use the same significant PCs as input
countDF <- RunUMAP(object = countDF, dims = 1:10)

pdf("tSNE.pdf")
DimPlot(object = countDF, reduction = "umap")
dev.off()

##########################################################################
## Save object at this point
saveRDS(countDF, file = "../output/pbmc_tutorial.rds")

##########################################################################
## Finding differentially expressed features (genes)
##Find markers that define clusters via differential expression
## "min.pct" argument requires a feature to be detected at a minimum percentage in either of the two groups of cells
## "thresh.test" argument requires a feature to be differentially expressed by some amount between the two groups
## both these arguments can be set to 0. But it will lead to a dramatic increase in time and will test a large number of 
## features that are unlikely to be highly discriminatory

## find all markers of cluster 1
cluster1.markers <- FindMarkers(object = countDF, ident.1 = 1, min.pct = 0.25)
head(x = cluster1.markers, n = 5)

## find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(object = countDF, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(x = cluster5.markers, n = 5)

## find markers for every cluster compared to all remaining cells, report only the positive ones
countDF.markers <- FindAllMarkers(object = countDF, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
countDF.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

cluster1.markers <- FindMarkers(object = pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

##  Plot DEGs across clusters
pdf("VlnPlot_DEG.pdf")
VlnPlot(object = countDF, features = c("MS4A1", "CD79A"))
dev.off()

## you can plot raw counts as well
pdf("VlnPlot_DEG_rawcount.pdf")
VlnPlot(object = countDF, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
dev.off()

## Plot gene expression across tSNE
pdf("FeaturePlot_tSNE.pdf")
FeaturePlot(object = pdf(VlnPlot_DEG.pdf"), features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))
dev.off()

## Generate heatmap for given cells and features
## Example shows plot of top 10 markers
top10 <- countDF.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

pdf("DoHeatmap.pdf")
DoHeatmap(object = countDF, features = top10$gene) + NoLegend()
dev.off()

##########################################################################
##Assigning cell type identity to clusters
## Order of cell types will mirror order of clusters
## cluster 0 = Memory CD4 T

new.cluster.ids <- c("Memory CD4 T", "Naive CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
    "NK", "DC", "Mk")
names(x = new.cluster.ids) <- levels(x = countDF)
countDF <- RenameIdents(object = countDF, new.cluster.ids)

pdf("Label_clusters.pdf")
DimPlot(object = countDF, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

saveRDS(pbmc, file = "../output/pbmc3k_final.rds")
