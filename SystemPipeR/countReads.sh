#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=1-00:15:00     # 1 day and 15 minutes
#SBATCH --output=countReads.stdout
#SBATCH --mail-user=useremail@address.com
#SBATCH --mail-type=ALL
#SBATCH --job-name="RNA-Seq readcounting"
#SBATCH -p intel # This is the default partition, you can use any of the following; intel, batch, highmem, gpu


# Print current date
date
Rscript countReads.R
