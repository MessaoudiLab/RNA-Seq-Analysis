
# small-RNA Expression Analysis Pipeline
# Expression analysis of cellular and plasma micro RNA
The following pipeline allows you to trim small RNA reads generated from Qia-Seq protocol.Reads are then aligned using very sensitive settings in Bowtie2 and counted using GenomicRanges. Differential gene expression analysis using edgeR.

## Step1: Trim files one by one using the following command: 
sbatch -p gpu --gres=gpu:1 --mem=10g --time=2:00:00 small_RNA_trimming.sh Input_Fastq_Prefix
## Step2: Generate the targets file, data, and results directory as described in systemPipeR
## Step3: Run small_RNA_alignment.R line by line
## Step4: Run small_RNA_readcounting.R line by line
## Step5: Run small_RNA_edgeR.R line by line
