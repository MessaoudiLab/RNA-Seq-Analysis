# RNA-Seq Analysis Protocol
## Basic linux Commands

ls: list files in current directory

cd: change directory

mkdir: make directory

cd ../: go back one directory

pwd -P: absolute path of present working directory

cp: copy

wget {url}: download data from internet



## 1. Setup working directory for RNAseq
- Open the Terminal and Login to UCR cluster
- Make a main for the project
```
mkdir “name of project”
```
- Create 2 sub-directories under the main directory titled “data” and “results”
```
mkdir data
mkdir results
```

In addition to the data and results directories, you'll need the following files in the main directory: slurm.tmpl, .Batchjobs.R, tophat.param, and targets.txt (see examples)

- .Batchjobs and slurm.tmpl are required for submitting jobs to the cluster. Copy these files exactly

- tophat.param is required for defining alignment parameters. Change GTF and FASTA reference genome information in this file. Everything else will be the same

- targets.file outlines the experimental design. See example as a guideline. $Filename lists the absolute pathway to fastq file, $SampleName is a name given to each fastq file and must be unique, $Factor is the experimental condition (i.e. STIM, baseline, day0, etc). <CMP> defines the comparisons you want to make - for example, DPI7-DPI0 will give differentially expressed genes at condition DPI7 relative to DPI0

To write/edit any of these files, use text editor nano (i.e. nano targets.txt) 

## 2. FASTQC
Run FASTQC on all fastq files. Information on how to run FASTQC can be found in the repository "NGS-Pre-Processing"

## 3. Trim files
Based on FASTQC metrics, Trim fastq files using Trim galore. Information on how to run Trim galore can be found in repository "NGS-Pre-Processing"

## 4. Alignment
- In the data directory, create symbolic link for the reference genome file, annotation file and index file. You should have 8 files: fasta, GTF, and 6 index files ending in .bt2 if using bowtie alignment.
```
ln -s absolute/path/to/reference/genome/files .
```

- Start R in the main directory 
- Load the required packages
```
R
library(systemPipeR)
library(GenomicFeatures)
```

- Read in the targets file and save it as object "targets"
```
targets <- read.delim("targets.txt", comment.char = "#")
targets
```

- Create "args" object, which saves information from tophat.param and targets.txt in order to run alignment
```
args <- systemArgs(sysma="tophat.param", mytargets="targets.txt")
moduleload(modules(args))
```
- Check alignment script for the first sequence
```
sysargs(args[1])
```
- Allocate the desired resources for alignment
- Note this assigns 20GB of memory and will run for max 20 hours
```
resources <- list(walltime="20:00:00", ntasks=1, ncpus=cores(args), memory="20G") 
```
- Submit the jobs to the cluster for alignment to take place
- Change Njobs to actual number of sequences to be aligned
```
reg <- clusterRun(args, conffile=".BatchJobs.R", template="slurm.tmpl", Njobs=18, runid="01", resourceList=resources)
```

Alignment typically takes ~6 hours and is usually run overnight. In order to monitor your alignment:
- quit R
```
q()
```
- In the linux environment use qstat
```
qstat | grep "username"
```
- Alignment is complete once jobs are cleared from the queue and tophat directories exist in the results directory

After alignment is complete, obtain alignment summary using readstats
* For all remaining steps, enter a highmem subnode

```
srun -p highmem --mem=100g --time=24:00:00 --pty bash -l
```

Start R in the main directory and Load the required packages again
```
R
library(systemPipeR)
library(GenomicFeatures)
```
- Read in targets file once again
```
targets <- read.delim("targets.txt", comment.char = "#")
targets
```
- Create an object
```
args <- systemArgs(sysma="tophat.param", mytargets="targets.txt")
file.exists(outpaths(args))
```
- Write the table for alignment stats
```
read_statsDF <- alignStats(args=args)
write.table(read_statsDF, file="results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")
```
### 4. Counting and Normalization
The .gtf file, .fa file, sqlite file and organism given down below are subject to change according to project. Include the appropriate files needed for your particular project.
- Load the required packages
```
library(GenomicFeatures)
```
- Create a txdb object
```
txdb <- makeTxDbFromGFF(file="data/Macaca_mulatta.MMUL_1.78.gtf", format="gtf", dataSource="ENSEMBL", organism="Macaca mulatta")
saveDb(txdb, file="./data/Macaca_mulatta.sqlite")
```
- The following performs read counting with summarizeOverlaps in parallel mode with multiple cores
```
library("GenomicFeatures"); library(BiocParallel)
txdb <- loadDb("./data/Macaca_mulatta.sqlite")
eByg <- exonsBy(txdb, by=c("gene"))
bfl <- BamFileList(outpaths(args), yieldSize=50000, index=character())
multicoreParam <- MulticoreParam(workers=8); register(multicoreParam); registered()
counteByg <- bplapply(bfl, function(x) summarizeOverlaps(eByg, x, mode="Union", ignore.strand=FALSE, inter.feature=TRUE, singleEnd=TRUE))
countDFeByg <- sapply(seq(along=counteByg), function(x) assays(counteByg[[x]])$counts)
```
- Wait until read counting is done, then run to make countDFeByg excel file
```
countDFeByg <- sapply(seq(along=counteByg), function(x) assays(counteByg[[x]])$counts)
rownames(countDFeByg) <- names(rowRanges(counteByg[[1]])); colnames(countDFeByg) <- names(outpaths(args))
write.table(countDFeByg, "results/countDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
```
- Generates RPKM normalized expression values from the countDFeByg file
```
rpkmDFeByg <- apply(countDFeByg, 2, function(x) returnRPKM(counts=x, ranges=eByg))
write.table(rpkmDFeByg, "results/rpkmDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
```
A counts file and RPKM normalized expression values are now generated and can be found in the "results" directory. These counts are normalized to remove biases introduced in the preparation steps such as length of the reads and sequencing depth (coverage) of a sample. During RPKM normalization, The total number of reads in a sample is divided by 1,000,000 (this is "per million" factor). Read counts are divided by this per million factor (normalizing for coverage giving you reads per millions). Then, these reads per million values are divided by the length of the gene in kilobases. Because RPKM normalization involves total number of reads in each sample (not just the counts of each individual reads), you might see RPKM values different between different samples/time points. This RPKM normalization step is done separately from EdgeR, which generates FC, p-value and FDR. EdgeR does its own normalization. 

### 5. Correlation Analysis
- The following computes the sample-wise Spearman correlation coefficients from the RPKM normalized expression values. 
```
library(ape)
rpkmDFeByg <- read.delim("./results/rpkmDFeByg.xls", row.names=1, check.names=FALSE)[,-19]
rpkmDFeByg <- rpkmDFeByg[rowMeans(rpkmDFeByg) > 50,]
d <- cor(rpkmDFeByg, method="spearman")
hc <- hclust(as.dist(1-d))
pdf("results/sample_tree.pdf")
plot.phylo(as.phylo(hc), type="p", edge.col="blue", edge.width=2, show.node.label=TRUE, no.margin=TRUE)
dev.off()
```
- The following computes the sample-wise Spearman correlation coefficients from the rlog (regularized-logarithm) transformed expression values generated with the DESeq2 package to make correlation dendrogram of samples for rlog values.
```
library(DESeq2)
countDF <- as.matrix(read.table("./results/countDFeByg.xls"))
colData <- data.frame(row.names=targets$SampleName, condition=targets$Factor)
dds <- DESeqDataSetFromMatrix(countData = countDF, colData = colData, design = ~ condition)
d <- cor(assay(rlog(dds)), method="spearman")
hc <- hclust(dist(1-d))
pdf("results/sample_tree_rlog.pdf")
plot.phylo(as.phylo(hc), type="p", edge.col=4, edge.width=3, show.node.label=TRUE, no.margin=TRUE)
dev.off()
```

### 6. PCA Plots
- Follow the commands down below to create group-wise and sample-wise PCA plots:
```
rld <- rlog(dds)
pdf("results/PCA_group.pdf")
plotPCA(rld)
dev.off()
colData <- data.frame(row.names=targets$SampleName, condition=targets$SampleName)
dds <- DESeqDataSetFromMatrix(countData = countDF, colData = colData, design = ~ condition)
rld <- rlog(dds)
pdf("results/PCA_sample.pdf")
plotPCA(rld)
dev.off()
```
- Follow the commands down below to create group-wise and sample-wise PCA plots based on vsd:
```
colData <- data.frame(row.names=targets$SampleName, condition=targets$Factor)
dds <- DESeqDataSetFromMatrix(countData = countDF, colData = colData, design = ~ condition)
vsd <- varianceStabilizingTransformation(dds)
pdf("results/PCA_group_vsd.pdf")
plotPCA(vsd)
dev.off()
colData <- data.frame(row.names=targets$SampleName, condition=targets$SampleName)
dds <- DESeqDataSetFromMatrix(countData = countDF, colData = colData, design = ~ condition)
vsd <- varianceStabilizingTransformation(dds)
pdf("results/PCA_sample_vsd.pdf")
plotPCA(vsd)
dev.off()
```

### 7. DEG Analysis with edgeR
- Follow the commands down below to run DEG analysis:
```
library(systemPipeR)
library(edgeR)
countDF <- read.delim("results/countDFeByg.xls", row.names=1, check.names=FALSE)
targets <- read.delim("targets.txt", comment="#")
cmp <- readComp(file="targets.txt", format="matrix", delim="-")
edgeDF <- run_edgeR(countDF=countDF, targets=targets, cmp=cmp[[1]], independent=TRUE, mdsplot="")
desc <- read.delim("/bigdata/messaoudilab/abotr002/References/Rhesus_Macaque/Rhesus_annotations.xls", row.names=1)
edgeDF <- cbind(edgeDF, desc[rownames(edgeDF),])
write.table(edgeDF, "./results/edgeRglm_allcomp.xls", quote=FALSE, sep="\t", col.names = NA)
edgeDF <- read.delim("results/edgeRglm_allcomp.xls", row.names=1, check.names=FALSE)
pdf("results/DEGcounts.pdf")
DEG_list <- filterDEGs(degDF=edgeDF, filter=c(Fold=2, FDR=5))
dev.off()
write.table(DEG_list$Summary, "./results/DEGcounts.xls", quote=FALSE, sep="\t", row.names=FALSE)
```
- Merge RPKM file (rpkmDFeByg.xls) and edgeR file (edgeRglm_allcomp.xls) (these files should be in the ‘results’ directory)
- nano rpkmDFeByg.xls and write a header “ENSEMBL_ID” on top of ensembl column. (Header is subject to change depending on project)
- nano edgeRglm_allcomp.xls and a header “RhesusEnsembl” on top of ensembl column. (Header is subject to change depending on project)
- Make sure you are in the ‘results’ directory when running the following in R
- Follow the commands down below to merge the data in the two files to create ‘edgeR_rpkm.xls’
```
edgeR <- read.delim("edgeRglm_allcomp.xls", sep="\t", header=TRUE)
rpkm <- read.delim("rpkmDFeByg.xls", sep="\t", header=TRUE)
edgeR_rpkm <- merge(edgeR, rpkm, by.x="RhesusEnsembl", by.y="ENSEMBL_ID", all=TRUE)
write.table(edgeR_rpkm, file="edgeR_rpkm.xls", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
```
