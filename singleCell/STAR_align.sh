#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=15G
#SBATCH --time=1-00:15:00     # 1 day and 15 minutes
#SBATCH --output=CD8_align.stdout
#SBATCH --mail-user=arive019@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="CD8_aria_align"
#SBATCH -p highmem # This is the default partition, you can use any of the following; intel, batch, highmem, gpu


# Print current date
date

# Change directory to where you submitted the job from, so that relative paths resolve properly
SLURM_SUBMIT_DIR=/bigdata/messaoudilab/arivera/Projects/SVV/CD8_aria/day7/Trimmed_reads/
cd $SLURM_SUBMIT_DIR

# Align
for i in *.read1_val_1.fq.gz; do
        /bigdata/messaoudilab/arivera/Scripts/STAR-2.5.3a/bin/Linux_x86_64/STAR --runMode alignReads --genomeLoad NoSharedMemory --limitBAMsortRAM 15000000000 --outSAMtype BAM SortedByCoordinate --outSAMmapqUnique 255 --readFilesCommand zcat --genomeDir /bigdata/messaoudilab/arivera/Reference_genomes/Rhesus_SVV/STAR_genome_62 --runThreadN 10 --readFilesIn $i ${i%.read1_val_1.fq.gz}.read2_val_2.fq.gz --outFileNamePrefix ${i%.read1_val_1.fq.gz}
done

