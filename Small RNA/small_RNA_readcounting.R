#The code below uses modified systemPiper for counting of reads mapped using a modified GTF.
#NB: Annotations from MiRbase are off by a 100bp - hence if the alignment has been done using a recent version of fasta, please use the ENSEMBL annotations

library(systemPipeR)
library(GenomicFeatures)
library(BiocParallel)

# Make sure the sqlite file is in the right destination
txdb <- loadDb("./data/Macaca_mulatta_smallRNA.sqlite")

eByg <- exonsBy(txdb, by="gene")
bfl <- BamFileList(outpaths(args), yieldSize=50000, index=character())
multicoreParam <- MulticoreParam(workers=4); register(multicoreParam); registered()

#Counting reads - assuming strand non-specific. inter.feature=FALSE as pre- and mature miRNA positions overlap
countemiR <- bplapply(bfl, function(x) summarizeOverlaps(eByg, x, mode="Union", ignore.strand=TRUE, inter.feature=FALSE, singleEnd=TRUE))
countDFmiR <- sapply(seq(along=countemiR), function(x) assays(countemiR[[x]])$counts)
rownames(countDFmiR) <- names(rowRanges(countemiR[[1]]))
colnames(countDFmiR) <- names(outpaths(args))
rpkmDFmiR <- apply(countDFmiR, 2, function(x) returnRPKM(counts=x, ranges=eByg))
                    
#Output counts and RPKM
write.table(countDFmiR, "./results/countDFmiR.xls", col.names=NA, quote=FALSE, sep="\t")
write.table(rpkmDFmiR, "./results/rpkmDFmiR.xls", col.names=NA, quote=FALSE, sep="\t")

#Sample outlier analysis
library(ape)
rpkmDFmiR <- read.delim("./results/rpkmDFmiR.xls", row.names=1, check.names=FALSE)[,-19]
rpkmDFmiR <- rpkmDFmiR[rowMeans(rpkmDFmiR) > 50,]
d <- cor(rpkmDFmiR, method="spearman")
hc <- hclust(as.dist(1-d))
pdf("./results/sample_tree.pdf")
plot.phylo(as.phylo(hc), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, no.margin=TRUE)
dev.off()

library(DESeq2)
countDFmiR <- as.matrix(read.table("./results/countDFmiR.xls"))
colData <- data.frame(row.names=targets$SampleName, condition=targets$Factor)
dds <- DESeqDataSetFromMatrix(countData = countDFmiR, colData = colData, design = ~ condition)
d <- cor(assay(rlog(dds)), method="spearman")
hc <- hclust(dist(1-d))
pdf("./results/sample_tree_rlog.pdf")
plot.phylo(as.phylo(hc), type="p", edge.col=4, edge.width=3, show.node.label=TRUE, no.margin=TRUE)
dev.off()
                   
#Group-wise and sample wise PCA following rlog transformation
rld <- rlog(dds)
pdf("./results/PCA_group.pdf")
plotPCA(rld)
dev.off()

colData <- data.frame(row.names=targets$SampleName, condition=targets$SampleName)
dds <- DESeqDataSetFromMatrix(countData = countDFmiR, colData = colData, design = ~ condition)
rld <- rlog(dds)
pdf("./results/PCA_sample.pdf")
plotPCA(rld)
dev.off()
                   
#Optional group-wise and sample-wise PCA following variance stabilization
colData <- data.frame(row.names=targets$SampleName, condition=targets$Factor)
dds <- DESeqDataSetFromMatrix(countData = countDFmiR, colData = colData, design = ~ condition)
vsd <- varianceStabilizingTransformation(dds)
pdf("results/PCA_group_vsd.pdf")
plotPCA(vsd)
dev.off()

colData <- data.frame(row.names=targets$SampleName, condition=targets$SampleName)
dds <- DESeqDataSetFromMatrix(countData = countDFmiR, colData = colData, design = ~ condition)
vsd <- varianceStabilizingTransformation(dds)
pdf("results/PCA_sample_vsd.pdf")
plotPCA(vsd)
dev.off()
                   
#Done
