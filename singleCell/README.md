
# Trim
Trim paired end reads with TrimGalore

- Trim GalorePE.sh is an sbatch script that will run TrimGalore on multiple files

- Make sure to cutomize TrimGalorePE.sh (i.e. output file name ".stdout", email, jobname, directory, file name, trimming parameters)

- This script will submit the job to be run in a subnode. All STDOUT will be redirected to a file called “my.stdout” as well as an email sent to the user when the status of the job changes.

```
sbatch TrimGalorePE.sh
```

# STAR alignment

## Install STAR
```
wget https://github.com/alexdobin/STAR/archive/2.7.0f.tar.gz
tar -xzf 2.7.0f.tar.gz
cd STAR-2.7.0f/source
make STAR
```

## Build indices
Note: For argument "--genomeDir" Create a directory where genome files will go
Note: For argument "--sjdbOverhang" read length refers to the sequences to be aligned

Basic code to build indices
```
STAR --runThreadN {number of cores} --limitGenomeGenerateRAM 92799002666 --runMode genomeGenerate --genomeDir /path/to/resulting/STAR/genome/ --genomeFastaFiles /path/to/genome/fasta/file --sjdbGTFfile /path/to/GTF/or/GFF --sjdbOverhang {read length - 1}
```
Example
```
/bigdata/messaoudilab/arivera/Scripts/STAR-2.7.0f/bin/Linux_x86_64/STAR --runThreadN 4 --limitGenomeGenerateRAM 92799002666 --runMode genomeGenerate --genomeDir /bigdata/messaoudilab/arivera/Reference_genomes/Rhesus_SVV/STAR_genome/ --genomeFastaFiles /bigdata/messaoudilab/arivera/Reference_genomes/Rhesus_SVV/Macaca_SVV.fasta --sjdbGTFfile /bigdata/messaoudilab/arivera/Reference_genomes/Rhesus_SVV/Macaca_mulatta_SVV.gtf --sjdbOverhang 62
```

## Align PE reads

- Basic code
```
STAR --runMode alignReads --genomeLoad NoSharedMemory --limitBAMsortRAM 15000000000 --outSAMtype BAM SortedByCoordinate --outSAMmapqUnique 255 --readFilesCommand zcat --genomeDir /absolute/path/to/STARgenome/ --runThreadN {number of cores} --readFilesIn {read1.fq.gz read2.fq.gz} --outFileNamePrefix {basename} 
```
- STAR_align.sh is an sbatch script that runs the alignment in a subnode.

- The script uses a For Loop to run the alignment on all samples that end in ".read1_val_1.fq.gz". Change the script accordingly, if your samples are named differently.

- In the script, the argument "--readFilesIn" is based on the assumption that samples end in ".read1_val_1.fq.gz". Change accordingly. 

- Make sure to cutomize STAR_align.sh (i.e. output file name ".stdout", email, jobname, directory)

- This script will submit the job into the cluster. All STDOUT will be redirected to a file called “my.stdout” as well as an email sent to the user when the status of the job changes.

## Usage
```
sbatch STAR_align.sh
```
# Counting using Rsubread
- Rsubread is an R package
- Make sure bam files exist in the working directory

Start R

```
library(Rsubread)
```
save bam file names as an object
```
fls <- dir(".",".out.bam")
```
Run featureCounts
```
counts <- featureCounts(files=fls, annot.ext="/bigdata/messaoudilab/arivera/Reference_genomes/Rhesus_SVV/Macaca_mulatta_SVV.gtf", isGTFAnnotationFile=TRUE, GTF.featureType="exon", GTF.attrType="gene_id", isPairedEnd=TRUE, requireBothEndsMapped=TRUE)
```
Add descriptions (i.e. hgnc symbols)
```
desc <- read.delim("/bigdata/messaoudilab/arivera/Reference_Macaque/Rhesus_annotations.xls", row.names=1)
counts2 <- counts$counts
counts2 <- cbind(counts2, desc[rownames(counts2),])
```
Write excel file
```
write.table(counts2, file="Counts_gene.xls", quote=FALSE, sep="\t")
```

Edit/clean up excel file to look like the example counts file (counts.txt)

```
countDF <- read.table("counts.txt", sep="\t", header=TRUE)
countDF2 <- countDF[,-1]
names <- countDF$GENE
rownames(countDF2)=make.names(names,unique=TRUE)
write.table(countDF2, file="Counts_uniquegenes.xls", quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")
```

The final counts file (e.g. Counts_uniquegenes.xls) will be the input for Seurat analysis


# Identify clusters and DEGs between clusters using Seurat
## Run Seurat_clusters_DEGs.R line by line

Will run QC metrics, tSNE clustering, identify differentially expressed genes between clusters, and assign cell type identities to clusters

# Identify DEGs between stimulated/infected versus control
## Run Seurat_stim_vs_control.R line by line

Integrates 2 datasets to look at gene expression differences between conditions (i.e. stimulated versus control, or infected vs. uninfected)

 - Aim to identify cell types that are present in both datasets
 
 - Obtain cell type markers that are conserved in both control and stimulated cells
 
 - Compare datasets to find cell-type specific responses to stimulation
 
