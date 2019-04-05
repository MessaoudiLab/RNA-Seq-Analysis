library(systemPipeR)
targets <- read.delim("targets.txt", comment.char = "#")
#Check if targets file is correct
targets
#Create args
args <- systemArgs(sysma="hisat2.param", mytargets="targets.txt")

library("GenomicFeatures")
library(BiocParallel)
txdb <- loadDb("./data/Homo_sapiens_hg38.95.sqlite")
eByg <- exonsBy(txdb, by=c("gene"))

# Read the bamfile list
bfl <- BamFileList(outpaths(args), yieldSize=50000, index=character())
multicoreParam <- MulticoreParam(workers=8); register(multicoreParam); registered()
counteByg <- bplapply(bfl, function(x) summarizeOverlaps(eByg, x, mode="Union", ignore.strand=FALSE, inter.feature=FALSE, singleEnd=TRUE))
countDFeByg <- sapply(seq(along=counteByg), function(x) assays(counteByg[[x]])$counts)
rownames(countDFeByg) <- names(rowRanges(counteByg[[1]])); colnames(countDFeByg) <- names(bfl)
rpkmDFeByg <- apply(countDFeByg, 2, function(x) returnRPKM(counts=x, ranges=eByg))
write.table(countDFeByg, "results_032619/countDFeByg_all_032719.xls", col.names=NA, quote=FALSE, sep="\t")
write.table(rpkmDFeByg, "results_032619/rpkmDFeByg_all_032719.xls", col.names=NA, quote=FALSE, sep="\t")
