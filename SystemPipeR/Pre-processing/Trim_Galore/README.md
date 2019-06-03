Trimming is done through Trim Galore which will quality filter and trim the reads.

# Require files
## trim_galore.sh
shell script that runs trim_galore

## trim_galore_directory.sh 
shell wrapper that runs trim_galore.sh on a directory. This script will also submit the job into a worker node (sbatch)

# Modifications
## trim_galore.sh
additional arguments can be added or modified depending on how you'd like sequences to be trimmed. 
See manual

## trim_galore_directory.sh
Line 15: "scriptsdir" should be equal to the directory (given as absolute path) where the script exists

Line 18: "for f in $directory*.gz": zipped fastq files should end in .gz. Othersize choose a basename/suffix that is common for all files
Check script to make sure it’s linked to the correct script (trim_galore.sh) before executing

# Usage
You only need to run trim_galore_directory.sh since it is already linked to trim_galore.sh
absolutepath/directory = The first argument is the directory containing fastq files you want to trim 
score = The second argument is the phred score cutoff
length = The third argument is the minimum length cutoff

```
sh absolutepath/trim_galore_directory.sh absolutepath/directory score length
```
example 

```
sh /bigdata/messaoudilab/abotr002/Scripts/trim_galore_directory.sh /bigdata/messaoudilab/abotr002/data/sequences/ 30 50
```
