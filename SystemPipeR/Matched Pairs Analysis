library(edgeR)
counts <- read.delim("countDFeByg_100719.xls", sep="\t", header=T)
head(counts)
ecoli <- subset(counts, select = c(ME23NS, ME23Ecoli, ME63NS, ME63Ecoli,  ME80NS, ME80Ecoli, ML06NS, ML06Ecoli))
head(ecoli)
rownames(ecoli) <- counts$X
head(ecoli)
keep <- rowSums(cpm(ecoli)>1) >= 4
head(keep)
ecoli_filtered <- ecoli[keep,]
nrow(ecoli)
nrow(ecoli_filtered)
y <- DGEList(ecoli)
keep <- rowSums(cpm(y)>1) >= 4
y <- y[keep,]
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y)
pdf("Leans_Ecoli_MDS.pdf")
plotMDS(y)
dev.off()
Baby <- factor(c(23,23,63,63,80,80,6,6))
Infection <- factor(c("U","I","U","I","U","I", "U", "I"))
data.frame(Sample=colnames(y),Baby,Infection)
design <- model.matrix(~Baby+Infection)
rownames(design) <- colnames(y)
design
y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
topTags(lrt)
summary(de <- decideTestsDGE(lrt))
detags <- rownames(y)[as.logical(de)]
pdf("Leans_Ecoli_MA.pdf")
plotSmear(lrt, de.tags=detags)
abline(h=c(-1, 1), col="blue")
dev.off()
topTags(lrt)
nrow(lrt)
write.table(topTags(lrt, n=18134), file="Leans_Ecoli_pairwise_100719.xls", sep="\t", quote=F)
leans_ecoli_deg <- as.data.frame(topTags(lrt, n=18134)) 

# Merging Annotation File with leans_ecoli_deg Object
annotations <- read.delim("hg38_annotations_100719.xls", sep="\t", header=T, row.names=1)
head(annotations)
leans_ecoli_deg <- as.data.frame(leans_ecoli_deg)
leans_ecoli_deg <- cbind(leans_ecoli_deg, annotations[rownames(leans_ecoli_deg),])
write.table(leans_ecoli_deg, file="Leans_Ecoli_pairwise_annotated_100719.xls", sep="\t", quote=FALSE)
