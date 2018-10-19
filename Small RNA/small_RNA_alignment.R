#SystemPipeR template for read alignment using bowtie with modified parameters.
library(systemPipeR)
library(GenomicFeatures)

#Read the targets file
targets <- read.delim("targets.txt", comment.char = "#")
#Check if targets file is correct
targets
#Create args
args <- systemArgs(sysma="bowtie-smallRNA.param", mytargets="targets.txt")
moduleload(modules(args))

#Check if alignment files exist
file.exists(outpaths(args))

#Assign resources. Note this assigns 1Gb of Ram per core. If ncpus is 4 then this will amount to 4Gb total
resources <- list(walltime="20:00:00", ntasks=1, ncpus=cores(args), memory="5G") 

#Submit jobs on the cluster
reg <- clusterRun(args, conffile=".BatchJobs.R", template="slurm.tmpl", Njobs=18, runid="01", resourceList=resources)
#Done!
