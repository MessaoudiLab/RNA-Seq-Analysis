# RNA-Seq Analysis Protocol

## 1. Setup working directory for RNAseq
- Open the Terminal and Login to UCR cluster
- Make a main directory for the project
```
mkdir “name of project”
```
- Create 2 sub-directories under the main directory titled “data” and “results”
```
mkdir data
mkdir results
```

In addition to the data and results directories, you'll need the following files in the main directory: slurm.tmpl, .Batchjobs.R, tophat.param or hisat2.param, and targets.txt (see examples)

- ".Batchjobs.R" and "slurm.tmpl" are required for submitting jobs to the cluster. Copy these files exactly as it is written in the examples provided.

- "tophat.param" or "hisat2.param" is required for defining alignment parameters depending on whether you use tophat or hisat2 for alignment. Change the GTF and FASTA reference genome information in this file (line 8 and 14, respectively). Everything else will be the same

- targets.txt outlines the experimental design. See example as a guideline. 
  - $Filename lists the absolute pathway to fastq file, 
  - $SampleName is a name given to each fastq file and must be unique, 
  - $Factor is the experimental condition (i.e. STIM, baseline, day0, etc). 
  - "CMP" defines the comparisons you want to make - for example, DPI7-DPI0 will give differentially expressed genes at condition DPI7   
  relative to DPI0

To write/edit any of these files, use text editor nano 
```
#example
nano targets.txt
```
## 2. Set up Data files
FASTQ files
-Move all fastq files to the data directory

Reference Genome
- Information on how to obtain reference FASTA and GTF file can be found in the repository "Reference Genomes"
- In the data directory, create symbolic link for the reference genome file, annotation file and index file. You should have 8 files: fasta, GTF, and 6 index files ending in .bt2 if using bowtie alignment.

```
ln -s {absolute/path/to/reference/genome/files} .
```

## 3. FASTQC
Run FASTQC on all fastq files. Information on how to run FASTQC can be found in the repository "NGS-Pre-Processing"

## 4. Trim files
Based on FASTQC metrics, Trim fastq files using Trim galore. Information on how to run Trim galore can be found in repository "NGS-Pre-Processing"

## 5. Alignment

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
- Change Njobs to actual number of sequences to be aligned (example shown is 18 jobs)
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
- Read in targets file  again
```
targets <- read.delim("targets.txt", comment.char = "#")
targets
```
- Create the args object again
```
args <- systemArgs(sysma="tophat.param", mytargets="targets.txt")
file.exists(outpaths(args))
```
- Write the table for alignment stats
```
read_statsDF <- alignStats(args=args)
write.table(read_statsDF, file="results/alignStats.xls", row.names=FALSE, quote=FALSE, sep="\t")
```
## 6. Counting and Normalization
The .gtf file, .fa file, sqlite file and organism given down below are subject to change according to project. Include the appropriate files needed for your  project.

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
```
- Wait until read counting is done, then write countDFeByg into an excel file
```
countDFeByg <- sapply(seq(along=counteByg), function(x) assays(counteByg[[x]])$counts)
rownames(countDFeByg) <- names(rowRanges(counteByg[[1]])); colnames(countDFeByg) <- names(outpaths(args))
write.table(countDFeByg, "results/countDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
```
- Generate RPKM normalized expression values from the countDFeByg file
```
rpkmDFeByg <- apply(countDFeByg, 2, function(x) returnRPKM(counts=x, ranges=eByg))
write.table(rpkmDFeByg, "results/rpkmDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
```
A counts file and RPKM normalized expression values are now generated and can be found in the "results" directory. These counts for each gene are normalized to account for differences in the library size of the sample as well as gene length. During RPKM normalization, The total number of reads in a sample is divided by 1,000,000 (this is the "per million" factor). Read counts are divided by this per million factor (normalizing for coverage giving you reads per millions). Then, these reads per million values are divided by the length of the gene in kilobases. This RPKM normalization step is done independently from EdgeR, which uses TMM for normalization. 

## 7. Correlation/clustering Analysis
Before running DEG analysis with edgeR, it's important to visualize the transcriptional profiles of each sample in order to determine if any samples are outliers and need to be removed. Hierarchical clustering or PCA clustering will give you an indication of the magnitude of transcriptional changes following edgeR analysis. An outlier should be removed if the library size or alignment rate is poor (less than 10 million reads and less than 50% alignment)

The next sections go over generating  graphs of either tree clusters or PCAs.

NOTE: The minimum required files you need to generate tree clusters or PCAs is the "targets.txt", "countDFeByg.xls", and "rpkmDFeByg.xls", which was generated in step 6.

### Sample-wise Spearman correlation using RPKM
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
### Sample-wise Spearman correlation using rlog transformed counts data
The following computes the sample-wise Spearman correlation coefficients from the rlog (regularized-logarithm) transformed expression values generated with the DESeq2 package to make correlation dendrogram of samples for rlog values.
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

### PCA Plot to create group-wise PCA plots using rlog normalized counts
```
rld <- rlog(dds)
pdf("results/PCA_group.pdf")
plotPCA(rld)
dev.off()
```
### PCA plot to create sample-wise PCA plots using rlog normalized counts
```
colData <- data.frame(row.names=targets$SampleName, condition=targets$SampleName)
dds <- DESeqDataSetFromMatrix(countData = countDF, colData = colData, design = ~ condition)
rld <- rlog(dds)
pdf("results/PCA_sample.pdf")
plotPCA(rld)
dev.off()
```
### PCA Plot to create group-wise PCA plots using vsd normalized counts
```
colData <- data.frame(row.names=targets$SampleName, condition=targets$Factor)
dds <- DESeqDataSetFromMatrix(countData = countDF, colData = colData, design = ~ condition)
vsd <- varianceStabilizingTransformation(dds)
pdf("results/PCA_group_vsd.pdf")
plotPCA(vsd)
dev.off()
```
### PCA plot to create sample-wise PCA plots using vsd normalized counts
```
colData <- data.frame(row.names=targets$SampleName, condition=targets$SampleName)
dds <- DESeqDataSetFromMatrix(countData = countDF, colData = colData, design = ~ condition)
vsd <- varianceStabilizingTransformation(dds)
pdf("results/PCA_sample_vsd.pdf")
plotPCA(vsd)
dev.off()
```

## 8. DEG Analysis with edgeR
Load required libraries
```
library(systemPipeR)
library(edgeR)
```
Read in raw counts file and targets.txt file
```
countDF <- read.delim("results/countDFeByg.xls", row.names=1, check.names=FALSE)
targets <- read.delim("targets.txt", comment="#")
```
Define comparisons to be made
```
cmp <- readComp(file="targets.txt", format="matrix", delim="-")
```
Run edgeR
```
edgeDF <- run_edgeR(countDF=countDF, targets=targets, cmp=cmp[[1]], independent=TRUE, mdsplot="")
```
Add descriptions to the edgeDF file
  - The output of edgeDF will give you statistical information for each gene which is identified by an ensembl ID. 
  - In order to annotate each ensembl ID with additional information (HGNC symbol, description, gene type, etc), use Biomart:     http://uswest.ensembl.org/biomart/martview/f63abf59cf05faef9ee3b9c3a175acf3
  - In Biomart: 1) choose the species: 2) choose information you want HGNC symbols, description, gene types, etc; 3) download excel file and save the corresponding directory where reference genome is located. See Cynomolgus_genes_5.0.94.txt as an example
  - in R, read in the annotation file and save as object "desc"
```
desc <- read.delim("/bigdata/messaoudilab/abotr002/References/Rhesus_Macaque/Rhesus_annotations.xls", row.names=1)
```
Bind annotation file to edgeDF file
```
edgeDF <- cbind(edgeDF, desc[rownames(edgeDF),])
```
Write annotated edgeDF file
```
write.table(edgeDF, "./results/edgeRglm_allcomp.xls", quote=FALSE, sep="\t", col.names = NA)
```

Obtain summary data
```
edgeDF <- read.delim("results/edgeRglm_allcomp.xls", row.names=1, check.names=FALSE)
pdf("results/DEGcounts.pdf")
DEG_list <- filterDEGs(degDF=edgeDF, filter=c(Fold=2, FDR=5))
dev.off()
write.table(DEG_list$Summary, "./results/DEGcounts.xls", quote=FALSE, sep="\t", row.names=FALSE)
```
## 9 Merge EdgeR file with RPKM file
- Merge RPKM file (rpkmDFeByg.xls) and edgeR file (edgeRglm_allcomp.xls) (these files should be in the ‘results’ directory)
- nano rpkmDFeByg.xls and write a header “RhesusEnsembl” on top of ensembl column. (Header is subject to change depending on project)
```
nano rpkmDFeByg.xls
```
- nano edgeRglm_allcomp.xls and a header “RhesusEnsembl” on top of ensembl column. (Header is subject to change depending on project)
```
nano edgeRglm_allcomp.xls
```

- Make sure you are in the ‘results’ directory when running the following in R
- Follow the example commands down below to merge the data in the two files to create ‘edgeR_rpkm.xls’
```
edgeR <- read.delim("edgeRglm_allcomp.xls", sep="\t", header=TRUE)
rpkm <- read.delim("rpkmDFeByg.xls", sep="\t", header=TRUE)
edgeR_rpkm <- merge(edgeR, rpkm, by.x="RhesusEnsembl", by.y="RhesusEnsembl", all=TRUE)
write.table(edgeR_rpkm, file="edgeR_rpkm.xls", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
```
