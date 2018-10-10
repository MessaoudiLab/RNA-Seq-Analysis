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
- move sequence files (.gz)  from “sequences” to “trim_sequences” directory using mv *.gz
