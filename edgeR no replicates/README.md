
# edgeR with no replicates
## Load package
```
library(edgeR)
```

## Read raw counts file
```
countDF <- read.delim("countfile", sep="\t", header=TRUE, row.names=1)
```

## Subset count file to only include the samples to be compared (i.e. control column and experimental column)
```
countDF2 <- subset(countDF, select=c(control, experimental))
counts <- as.matrix(countDF2)
```
```
bcv <- 0.4
y <- DGEList(counts=counts, group=1:2)
et <- exactTest(y, dispersion=bcv^2)
```

## Add descriptions (same as in systemPipeR pipeline)
```
desc <- read.delim("bigdata/messaoudilab/arivera/Reference_Macaque/Rhesus_annotations.xls", row.names=1)
```

## Run edgeR
```
edgeDF <- cbind(et$table, desc[rownames(et$table),])
write.table(edgeDF, "edgeRglm_allcomp.xls", quote=FALSE, sep="\t", col.names=NA)
```
