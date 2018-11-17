
# Trim
Trim paired end reads with TrimGalore

Make sure to cutomize TrimGalorePE.sh (i.e. output file name ".stdout", email, jobname, directory, file name, trimming parameters)

```
sbatch TrimGalorePE.sh
```

# STAR alignment
Build indices
```
STAR --runThreadN {number of cores} --runMode genomeGenerate --genomeDir /path/to/resulting/STAR/genome/ --genomeFastaFiles /path/to/genome/fasta/file --sjdbGTFfile /path/to/GTF/or/GFF --sjdbOverhang {read length - 1}
```
Example
```
/bigdata/messaoudilab/arivera/Scripts/STAR-2.5.3a/bin/Linux_x86_64/STAR --runThreadN 4 --limitGenomeGenerateRAM 92799002666 --runMode genomeGenerate --genomeDir /bigdata/messaoudilab/arivera/Reference_genomes/Rhesus_SVV/STAR_genome/ --genomeFastaFiles /bigdata/messaoudilab/arivera/Reference_genomes/Rhesus_SVV/Macaca_SVV.fasta --sjdbGTFfile /bigdata/messaoudilab/arivera/Reference_genomes/Rhesus_SVV/Macaca_mulatta_SVV.gtf --sjdbOverhang 62
```
