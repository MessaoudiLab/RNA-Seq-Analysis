#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=1-00:15:00     # 1 day and 15 minutes
#SBATCH --output=trim_CD8aria.stdout
#SBATCH --mail-user=riveraa8@uci.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="CD8_aria_trim"
#SBATCH -p highmem # This is the default partition, you can use any of the following; intel, batch, highmem, gpu

#The above job will request 1 node, 10 task (1 cpu core per task), 10GB of memory (1GB per task), for 1 day and 15 minutes. 
#All STDOUT will be redirected to a file called “my.stdout” as well as an email sent to the user when the status of the job changes.

# Print current date
date

# Load fastqc
module load trim_galore/0.4.2

# Change directory to where you submitted the job from, so that relative paths resolve properly
SLURM_SUBMIT_DIR=/bigdata/messaoudilab/arivera/Projects/SVV/CD8_aria/day7/data/Run292-fixed/
cd $SLURM_SUBMIT_DIR

for f1 in *.read1.fastq.gz
do
  	i=`echo $f1 | cut -f 1 -d "."`
        trim_galore --clip_R1 13 --clip_R2 15 --three_prime_clip_R1 2 --three_prime_clip_R2 2 -q 20 -length 52 --paired --trim1 --fastqc -o /bigdata/messaoudilab/arivera/Projects/SVV/CD8_aria/day7/data/Run292-fixed/ $i.read1.fastq.gz $i.read2.fastq.gz
done
