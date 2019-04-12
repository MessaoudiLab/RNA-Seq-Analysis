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

# TPM Conversion 
# Genes length table need to be generated only once per annotation file
# The original version of the script below also calculates GC content (https://github.com/dpryan79/Answers/blob/master/SEQanswers_42420/GTF2LengthGC.R)                  
library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)
GTFfile = "./data/Macaca_fascicularis.Macaca_fascicularis_5.0.94.gtf"                    
#Load the annotation and reduce it                    
GTF <- import.gff(GTFfile, format="gtf", feature.type="exon")             
grl <- reduce(split(GTF, elementMetadata(GTF)$gene_id))
reducedGTF <- unlist(grl, use.names=T)                    
elementMetadata(reducedGTF)$gene_id <- rep(names(grl), elementNROWS(grl))
elementMetadata(reducedGTF)$widths <- width(reducedGTF)                    
#Create a list of the ensembl_id/length                    
calc_length <- function(x) {   
  sum(elementMetadata(x)$widths)
}
output <- t(sapply(split(reducedGTF, elementMetadata(reducedGTF)$gene_id), calc_length))
output <- t(output)                    
colnames(output) <- c("Length")                    
write.table(output, file="./results/lengths.txt", sep="\t")
genes <- read.delim("./results/lengths.txt", sep="\t", header=T)
genes$Length <- genes$Length/1000
# Table below gives you length of genes (exons) in kpb.                    
write.table(genes, file="./results/lengths.txt", sep="\t")
# Read the raw counts and gene lengths files, and TPM formula
counts <- read.delim("results/countDFeByg.xls", row.names=1, check.names=FALSE)
genelengths <- read.delim("./results/lengths.txt", sep="\t", header=T)
tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}
tpms <- apply(counts, 2, function(x) tpm(x, genelengths$Length))
write.table(tpms, "./results/TPM.xls", quote=FALSE, sep="\t", col.names = NA)
              
# You can check if TPM values are calculated correclty using               
# 1) Because TPM is normalize across different samples, column sums are equal among samples (1e+06)
colSums(tpms)
# 2) Similarly, sample means are equal among samples               
colMeans(tpms)
# 3) Generate normalized RPKM and compare to that generated with SystemPipeR              
rpkm <- function(counts, lengths) {
  rate <- counts / lengths 
  rate / sum(counts) * 1e6
}
rpkms <- apply(counts, 2, function(x) rpkm(x, genelengths$Length))              
write.table(rpkms, "./results/RPKM.xls", quote=FALSE, sep="\t", col.names = NA)
