##########################################################################
## Load libraries
library(Seurat)
library(cowplot)

##########################################################################
## Set up seurat objects
## Use raw count data
ctrl.data <- read.table(file = "../data/immune_control_expression_matrix.txt.gz", sep = "\t")
stim.data <- read.table(file = "../data/immune_stimulated_expression_matrix.txt.gz", sep = "\t")

## Set up control object
## QC, normalize and find variable features
ctrl <- CreateSeuratObject(counts = ctrl.data, project = "IMMUNE_CTRL", min.cells = 3)
ctrl$stim <- "CTRL"

ctrl[["percent.mt"]] <- PercentageFeatureSet(object = ctrl, pattern = "^MT.")

pdf("VlnPlot.pdf")
VlnPlot(object = ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

pdf("FeaturePlots.pdf")
plot1 <- FeatureScatter(object = ctrl, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = ctrl, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()

ctrl <- subset(x = ctrl, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
ctrl <- NormalizeData(object = ctrl, normalization.method = "LogNormalize", scale.factor = 10000)
ctrl <- FindVariableFeatures(object = ctrl, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(x = ctrl)
ctrl <- ScaleData(object = ctrl, features = all.genes)




## Set up stimulated object
## QC, normalize and find variable features
stim <- CreateSeuratObject(counts = stim.data, project = "IMMUNE_STIM", min.cells = 3)
stim$stim <- "STIM"
stim[["percent.mt"]] <- PercentageFeatureSet(object = stim, pattern = "^MT.")

pdf("VlnPlot.pdf")
VlnPlot(object = stim, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

pdf("FeaturePlots.pdf")
plot1 <- FeatureScatter(object = stim, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = stim, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()

stim <- subset(x = stim, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
stim <- NormalizeData(object = stim, normalization.method = "LogNormalize", scale.factor = 10000)
stim <- FindVariableFeatures(object = stim, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(x = stim)
stim <- ScaleData(object = stim, features = all.genes)


##########################################################################
## Perform integration
immune.anchors <- FindIntegrationAnchors(object.list = list(ctrl, stim), dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)

##########################################################################
## Perform an integrated analysis
DefaultAssay(object = immune.combined) <- "integrated"

## Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(object = immune.combined, verbose = FALSE)
immune.combined <- RunPCA(object = immune.combined, npcs = 30, verbose = FALSE)

## t-SNE and Clustering
immune.combined <- RunUMAP(object = immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(object = immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

## Visualization
pdf("tSNE.pdf")
p1 <- DimPlot(object = immune.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(object = immune.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
dev.off()

## Visualize the two conditions side-by-side
pdf("DimPlot.pdf")
DimPlot(object = immune.combined, reduction = "umap", split.by = "stim")
dev.off()

##########################################################################
## Identify conserved cell type markers
## Performs differential gene expression testing for each group.
## Example shows genes that are conserved irrespective of stimulation condition in cluster 7 (NK cells)
DefaultAssay(object = immune.combined) <- "RNA"
nk.markers <- FindConservedMarkers(object = immune.combined, ident.1 = 7, grouping.var = "stim", verbose = FALSE)
head(x = nk.markers)

## Gene markers for each cluster can be used to annotate each cluster as specific cell types
FeaturePlot(object = immune.combined, features = c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A", "CCL2", "PPBP"), min.cutoff = "q9")
immune.combined <- RenameIdents(object = immune.combined, `0` = "CD14 Mono", `1` = "CD4 Naive T", `2` = "CD4 Memory T", `3` = "CD16 Mono", `4` = "B", `5` = "CD8 T", `6` = "T activated", `7` = "NK", `8` = "DC", `9` = "B Activated", `10` = "Mk", `11` = "pDC", `12` = "Eryth", `13` = "Mono/Mk Doublets")

##Visualize
pdf("tSNE_clusters_labeled.pdf")
DimPlot(object = immune.combined, label = TRUE)
dev.off()

##Plot 2-3 strong marker genes for each cluster
Idents(object = immune.combined) <- factor(x = Idents(object = immune.combined), levels = c("Mono/Mk Doublets", "pDC", "Eryth", "Mk", "DC", "CD14 Mono", "CD16 Mono", "B Activated", "B", "CD8 T", "NK", "T activated", "CD4 Naive T", "CD4 Memory T"))
markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5", "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1", "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ")

pdf("DotPlot_genemarkers.pdf")
DotPlot(object = immune.combined, features = rev(x = markers.to.plot), cols = c("blue", "red"), dot.scale = 8, split.by = "stim") + RotatedAxis()
dev.off()

##########################################################################
## Identify differential expressed genes across conditions

## Plot average expression of both stimulated and control cells and look for genes that are visual outliers
## Takes average expression of both groups and generate scatter plots, highlighting genes that exhibit dramatic responses to stimulation
## Example is shown for one of the clusters that was labeled as T cells
t.cells <- subset(x = immune.combined, idents = "CD4 Naive T")
Idents(object = t.cells) <- "stim"
avg.t.cells <- log1p(x = AverageExpression(object = t.cells, verbose = FALSE)$RNA)
avg.t.cells$gene <- rownames(x = avg.t.cells)
genes.to.label = c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8")

pdf("scatter_STIM_CTRL.pdf")
p1 <- ggplot(avg.t.cells, aes(CTRL, STIM)) + geom_point() + ggtitle("CD4 Naive T Cells")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
plot_grid(p1)
dev.off()


## Find genes that are different between stimulated and control groups
## Example looks for DEGs between stimulated and control CD4 Naive T cells
immune.combined$celltype.stim <- paste(Idents(object = immune.combined), immune.combined$stim, sep = "_")
immune.combined$celltype <- Idents(object = immune.combined)
Idents(object = immune.combined) <- "celltype.stim"
response <- FindMarkers(object = immune.combined, ident.1 = "CD4 Naive T_STIM", ident.2 = "CD4 Naive T_CTRL", verbose = FALSE)
head(x = response, n = 15)
write.table(response, file="DEGlist.xls", header=TRUE, quote=FALSE, col.names=TRUE, row.names=FALSE)

## Visualize difference in gene expression for select genes between control and stimulated in tSNE
pdf(FeaturePlot_DEGs.pdf")
FeaturePlot(object = immune.combined, features = c("CD3D", "GNLY", "IFI6"), split.by = "stim", max.cutoff = 3, cols = c("grey", "red"))
dev.off()

## Visualize difference in gene expression in violin
pdf("Vln_DEGs.pdf")
plots <- VlnPlot(object = immune.combined, features = c("LYZ", "ISG15", "CXCL10"), split.by = "stim", group.by = "celltype", pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)
dev.off()
