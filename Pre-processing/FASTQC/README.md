Used to generate quality reports of fastq data

# Required files
1. fastqc_dir.sh and 2. fastqc.sh (see examples)

fastqc.sh is a shell script that will load fastqc and run it on zipped fastq files. 
fastqc_dir.sh is a shell wrapper that will execute "fastqc.sh" on a directory containing fastq files

# Modifications
## fastqc.sh
fastqc.sh is enough to run basic FASTQC
but if needed, you can change/add additional arguments. See Fastqc manual

## fastqc_dir.sh
Line 12: In the script make sure to change "scriptsdir" to be equal to the directory (given as absolute path) where the script exists

Line 15: for f in `$directory`*.gz: zipped fastq files should end in .gz. Otherwize choose a basename/suffix that is common for all files

# Usage

```
sh absolutepath/fastqc_dir.sh absolutepath/directory absolutepath/directory
```
sh = will execute shell script
fastqc_dir.sh - input absolute path and name of script
directory - absolute path of directory you will run fastqc_dir.sh on and output directory

example
```
sh /bigdata/messaoudilab/abotr002/Scripts/fastqc_dir.sh /bigdata/messaoudilab/abotr002/data/sequences/ /bigdata/messaoudilab/abotr002/data/sequences/
```

# Interpreting results

- Analyze generated fastQC; check total sequences, sequence length (101 for HiSeq or 75 for NextSeq), %GC (low 40s)
- Check quality control (Phred score: 30 or higher)
- See FASTQC_interpretation.pdf
