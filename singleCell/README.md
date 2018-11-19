
# Trim
Trim paired end reads with TrimGalore

- Trim GalorePE.sh is an sbatch script that will run TrimGalore on multiple files using a for loop

- Make sure to cutomize TrimGalorePE.sh (i.e. output file name ".stdout", email, jobname, directory, file name, trimming parameters)

- This script will submit the job into the cluster. All STDOUT will be redirected to a file called “my.stdout” as well as an email sent to the user when the status of the job changes.

```
sbatch TrimGalorePE.sh
```

# STAR alignment
## Build indices
```
STAR --runThreadN {number of cores} --runMode genomeGenerate --genomeDir /path/to/resulting/STAR/genome/ --genomeFastaFiles /path/to/genome/fasta/file --sjdbGTFfile /path/to/GTF/or/GFF --sjdbOverhang {read length - 1}
```
Example
```
/bigdata/messaoudilab/arivera/Scripts/STAR-2.5.3a/bin/Linux_x86_64/STAR --runThreadN 4 --limitGenomeGenerateRAM 92799002666 --runMode genomeGenerate --genomeDir /bigdata/messaoudilab/arivera/Reference_genomes/Rhesus_SVV/STAR_genome_62/ --genomeFastaFiles /bigdata/messaoudilab/arivera/Reference_genomes/Rhesus_SVV/Macaca_SVV.fasta --sjdbGTFfile /bigdata/messaoudilab/arivera/Reference_genomes/Rhesus_SVV/Macaca_mulatta_SVV.gtf --sjdbOverhang 62
```

## Align PE reads

- Basic code
```
STAR --runMode alignReads --genomeLoad NoSharedMemory --limitBAMsortRAM 15000000000 --outSAMtype BAM SortedByCoordinate --outSAMmapqUnique 255 --readFilesCommand zcat --genomeDir /absolute/path/to/STARgenome/ --runThreadN {number of cores} --readFilesIn {read1.fq.gz read2.fq.gz} --outFileNamePrefix {basename} 
```
- STAR_align.sh is an sbatch script that runs the code above in a for loop to align multiple paired end reads.

- Make sure to cutomize STAR_align.sh (i.e. output file name ".stdout", email, jobname, directory)

- This script will submit the job into the cluster. All STDOUT will be redirected to a file called “my.stdout” as well as an email sent to the user when the status of the job changes.

## Usage
```
sbatch STAR_align.sh
```
# Counting using Rsubread
- Rsubread is an R package
- Make sure bam files exist in the working directory

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
