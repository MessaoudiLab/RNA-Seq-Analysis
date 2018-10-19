# RNA-Seq Analysis Protocol
## Basic Commands

### ls: list
### cd: change directory
### mkdir: make directory
### cd ../: go back one directory
### pwd: present working directory
### cp: copy
### wget: retrieves information from internet



### 1. Making Directories
- Open the Terminal and Login
- Make a parent directory for the project
```
mkdir “name of project”
```
- Create 2 sub-directories under the parent directory titled “data” and “results”
```
mkdir data
mkdir results
```

In addition to the data and results directories, it is important for most projects to include parameters such as slurm.tmpl, targets.txt, and tophat.param. Include these under the parent directory by copying the parameters from another lab member. Note:To open parameter files, use the command nano. Ex:  nano slurm.tmpl or nano targets.txt

### 2. Downloading Data on Terminal
Libraries are sequenced through the website http://illumina.bioinfo.ucr.edu/ht/

- Once in the parent directory, change directories to the sub-directory “data” 
```
cd data
```
- In the data directory create a download script
```
nano “projectname”_download.sh
```
- Use the wget command to link sequence addresses from the illumina.bioinfo website to the terminal

Example:
wget http://illumina.bioinfo.ucr.edu/illumina_runs/358/flowcell34_lane99_pair1_CGATGT.fastq.gz

- Run the download script
```
sh“projectname”_download.sh
```

### 3. FastQC Files
- After sequences are downloaded, generate fastqc files using fastqc_dir.sh script found in    /bigdata/messaoudilab/arivera/Scripts/fastqc/fastqc_dir.sh
- Run the script 
```
sh “absolute path for fastqc_dir.sh” /”absolute pathway  for the ‘sequences’ directory” /”absolute pathway for the output directory”
```

Example: sh /bigdata/messaoudilab/abotr002/Scripts/fastqc_dir.sh /bigdata/messaoudilab/abotr002/data/sequences/ /bigdata/messaoudilab/abotr002/data/sequences/

- Analyze generated fastQC; check total sequences, sequence length (101 or 75 for NextSeq), %GC (low 40s)
- Check quality control (Phred score: 30 or higher)

### 4. Trimming
Once you have checked for quality control, you must trim sequences if certain regions vary greatly from data in order to not distort sequence assembly. Trimming is done through Trim Galore which will quality filter and trim the reads. You must use a script and signify the phred score and number of base pairs wanted at the end of the command line.

- In data directory, make new directory: 
```
mkdir trim_sequences
```
- move sequence files (.gz)  from “sequences” to “trim_sequences” directory 
```
mv *.gz ../trim_sequences
```
- Check trim_galore_directory.sh script to make sure it’s linked to the correct script (trim_galore.sh) before loading the following. 
```
sh /”absolute pathway for trim_galore_directory.sh” /”absolutepathwayfor‘sequences’directory 30 75”
```
If trimming NextSeq data, Run: 
```
sh /”absolute pathway for trim_galore_directory.sh” /”absolutepathwayfor‘sequences’directory 30 50”
```
- Check if trimming is done
```
qstat | grep “username” (if done, nothing will show)
```
- When trimming is complete, check FastQC reports.

### 5. Alignment
- In the data directory, create symbolic link for the reference genome file (ending in .fasta), index files (ending in .bt2), and annotation file (ending in .gtf) from arivera/viral_genomes/Rhesus_“ virus name ” using:
 ln -s “absolute pathway where the fa and gtf files are”/”file name.fa or .gtf”
- Copy slurm.tmpl, .BatchJobs.R, tophat.param and targets.txt files.nano tophat.param script. Make sure you are using the correct annotation file (ending in .gtf), and correct reference genome file (ending in .fa or .fasta).
- Edit targets.txt as follows: nano targets.txt
- Run R in the main directory where the alignment files are (targets.txt file, tophat.param and etc.)
- Load the required packages
```
R
library(systemPipeR)
library(GenomicFeatures)
```
- Name and read in the targets file
```
targets <- read.delim("targets.txt", comment.char = "#")
targets
```
- Create object in order to run alignment
```
args <- systemArgs(sysma="tophat.param", mytargets="targets.txt")
moduleload(modules(args))
sysargs(args[1])
```
- Allocate the desired resources for alignment
```
resources <- list(walltime="20:00:00", ntasks=1, ncpus=cores(args), memory="20G") # note this assigns 1Gb of Ram per core. If ncpus is   4, then this will amount to 4Gb total
```
- Submit the jobs to the cluster for alignment to take place
```
reg <- clusterRun(args, conffile=".BatchJobs.R", template="slurm.tmpl", Njobs=18, runid="01", resourceList=resources)
waitForJobs(reg)
```

Alignment typically takes 10 hours. In order to ensure that your sequences are being aligned:
- use the command qstat | grep arivera
- Jobs should be running (R)

### 6. Alignment Statistics
- To get alignment stats and counts faster, run a sub node using srun
```
srun  --mem=20gb --cpus-per-task 1 --ntasks 1 --time 10:00:00 --pty bash -l
```
- Load the required packages 
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

### 7. Counting and Normalization
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

### 8. Correlation Analysis
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

### 9. PCA Plots
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

### 10. DEG Analysis with edgeR
- Follow the commands down below to run DEG analysis:
```
library(edgeR)
countDF <- read.delim("results/countDFeByg_OR_colon.txt
", row.names=1, check.names=FALSE)
targets <- read.delim("targetsORTC.txt", comment="#")
cmp <- readComp(file="targetsORTC.txt", format="matrix", delim="-")
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
- nano rpkmDFeByg.xls and write a header “ENSEMBL_ID” on top of ensembl column. R
- nano edgeRglm_allcomp.xls and a header “RhesusEnsembl” on top of ensembl column.
- Make sure you are in the ‘results’ directory when running the following in R
- Follow the commands down below to merge the data in the two files to create ‘edgeR_rpkm.xls’
```
edgeR <- read.delim("edgeRglm_allcomp.xls", sep="\t", header=TRUE)
rpkm <- read.delim("rpkmDFeByg.xls", sep="\t", header=TRUE)
edgeR_rpkm <- merge(edgeR, rpkm, by.x="RhesusEnsembl", by.y="ENSEMBL_ID", all=TRUE)
write.table(edgeR_rpkm, file="edgeR_rpkm_no_TC45768_noloc.xls", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
```
