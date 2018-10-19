library(systemPipeR)
library(edgeR)
countDFmiR <- read.delim("./results/countDFmiR.xls", row.names=1, check.names=FALSE)
targets <- read.delim("targets.txt", comment="#")
cmp <- readComp(file="targets.txt", format="matrix", delim="-")

#Running edgeR
edgeDF <- run_edgeR(countDF=countDFmiR, targets=targets, cmp=cmp[[1]], independent=TRUE, mdsplot="")
desc <- read.delim("Rhesus_annotations.xls", row.names=1)
edgeDF <- cbind(edgeDF, desc[rownames(edgeDF),])
write.table(edgeDF, "./results/edgeRglm_allcomp.xls", quote=FALSE, sep="\t", col.names = NA)

edgeDF <- read.delim("results/edgeRglm_allcomp.xls", row.names=1, check.names=FALSE)
pdf("results/DEGcounts.pdf")
DEG_list <- filterDEGs(degDF=edgeDF, filter=c(Fold=2, FDR=5))
dev.off()
