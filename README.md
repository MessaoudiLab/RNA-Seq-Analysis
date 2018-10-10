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
- Make a parent directory for the project using mkdir “name of project”
- Create 2 sub-directories under the parent directory titled “data” and “results”
- Do this by mkdir data and mkdir results

In addition to the data and results directories, it is important for most projects to include parameters such as slurm.tmpl, targets.txt, and tophat.param. Include these under the parent directory by copying the parameters from another lab member. Note:To open parameter files, use the command nano. Ex:  nano slurm.tmpl or nano targets.txt

### 2. Transferring Data to Terminal
Libraries are sequenced through the website http://illumina.bioinfo.ucr.edu/ht/

- Once in the parent directory, change directories to the sub-directory “data” by cd data
- In the data directory create a download script by nano “projectname”_download.sh
- Use the wget command to link sequence addresses from the illumina.bioinfo website to the terminal

Example:
wget http://illumina.bioinfo.ucr.edu/illumina_runs/358/flowcell34_lane99_pair1_CGATGT.fastq.gz

- Now run the download script by sh“projectname”_download.sh

### 3. FastQC Files
- After sequences are downloaded, generate fastqc files using fastqc_dir.sh script found in    /bigdata/messaoudilab/arivera/Scripts/fastqc/fastqc_dir.sh
- Run the script as follows: sh “absolute path for fastqc_dir.sh” /”absolute pathway  for the ‘sequences’ directory” /”absolute pathway for the output directory”

Example: sh /bigdata/messaoudilab/abotr002/Scripts/fastqc_dir.sh /bigdata/messaoudilab/abotr002/data/sequences/ /bigdata/messaoudilab/abotr002/data/sequences/

- Analyze generated fastQC; check total sequences, sequence length (101 or 75 for NextSeq), %GC (low 40s)
- Check quality control (Phred score: 30 or higher)

### 4. Trimming
Once you have checked for quality control, you must trim sequences if certain regions vary greatly from data in order to not distort sequence assembly. Trimming is done through Trim Galore which will quality filter and trim the reads. You must use a script and signify the phred score and number of base pairs wanted at the end of the command line.

- In data directory, make new directory: mkdir trim_sequences
- move sequence files (.gz)  from “sequences” to “trim_sequences” directory using mv *.gz ../trim_sequences
- Check trim_galore_directory.sh script to make sure it’s linked to the correct script (trim_galore.sh) before loading the following. Run: sh /”absolute pathway for trim_galore_directory.sh” /”absolutepathwayfor‘sequences’directory 30 75”
If trimming NextSeq data, Run: sh /”absolute pathway for trim_galore_directory.sh” /”absolutepathwayfor‘sequences’directory 30 50”
- Check if trimming is done: qstat | grep “username” (if done, nothing will show)
- When trimming is complete, check FastQC reports.

### 5. Alignment
- In the data directory, create symbolic link for the reference genome file (ending in .fasta), index files (ending in .bt2), and annotation file (ending in .gtf) from arivera/viral_genomes/Rhesus_“ virus name ” using:
 ln -s “absolute pathway where the fa and gtf files are”/”file name.fa or .gtf”
- Copy slurm.tmpl, .BatchJobs.R, tophat.param and targets.txt files.nano tophat.param script. Make sure you are using the correct annotation file (ending in .gtf), and correct reference genome file (ending in .fa or .fasta).
- Edit targets.txt as follows: nano targets.txt
- Run R in the main directory where the alignment files are (targets.txt file, trim.param and etc.)
- Follow the commands down below to complete the alignment
```
R
library(systemPipeR)
library(GenomicFeatures)
targets <- read.delim("targets.txt", comment.char = "#")
targets
args <- systemArgs(sysma="tophat.param", mytargets="targets.txt")
moduleload(modules(args))
sysargs(args[1])
resources <- list(walltime="20:00:00", ntasks=1, ncpus=cores(args), memory="20G") # note this assigns 1Gb of Ram per core. If ncpus is 4, then this will amount to 4Gb total
reg <- clusterRun(args, conffile=".BatchJobs.R", template="slurm.tmpl", Njobs=18, runid="01", resourceList=resources)
waitForJobs(reg)

Alignment typically takes 10 hours. In order to ensure that your sequences are being aligned:
- use the command qstat | grep arivera
- Jobs should be running (R)

### 6. Alignment Statistics





