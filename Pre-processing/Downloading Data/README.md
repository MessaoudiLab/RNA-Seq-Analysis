# Downloading Data to Terminal
Libraries are sequenced through the website http://illumina.bioinfo.ucr.edu/ht/

In the desired directory, create a download script (see example)
```
nano example_download.sh
```
In blank text file, write the download code for each sample or directory 

 - Use the wget command to download sequences from URL address 
 - Example:
```
wget http://illumina.bioinfo.ucr.edu/illumina_runs/358/flowcell34_lane99_pair1_CGATGT.fastq.gz
```

Run the  script
```
sh “projectname”_download.sh
