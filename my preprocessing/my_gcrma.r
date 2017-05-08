# A R script that performs GCRMA preprocessing to selected Cel Files
# Author: Isidro C. Medina Jr., 2008
# Input: Cel files              Individual Cel Files containing expression values
# Output: gcrma_data.txt        Expression Set (tab delimited)

#load the required libraries
library(gcrma); library(affy); library(hu6800cdf); # library(hu6800); library(hu6800probe); 

# Set the Working directory, the working directory should contain all the Affy CEL files
setwd = ("E:\\My Thesis\\Training Set\\01 INPUT Cel Files\\Normal")

## Compute for affinity info
affinity.info.hu6800 <- compute.affinities("hu6800")
#save(affinity.info.hu6800,file = "affinity.hu6800.RData")
load("affinity.hu6800.RData")

# ReadAffy(), Reads all *.CEL (*.cel) files in your current working directory
system.time({
  #targets <- readTargets("targets.txt")         # Import targets information. One entry per line.
  #data <- ReadAffy(filenames=targets$FileName)  # Import expression raw data stores them as AffyBatch object 'data'.
  data <- ReadAffy(celfile.path="E:\\My Thesis\\Training Set\\01 INPUT Cel Files\\Normal")
})
#save(data,file = "raw.RData")

# [3] GC-RMA Normalization
ptm <- proc.time()
   eset.gcrma <- gcrma(data,affinity.info=affinity.info.hu6800)
proc.time() - ptm
save(eset.gcrma,file = "gcrma_normal.RData")

#export the GC-RMA expression values to text file
unlink("gcrma_normal.txt")       #delete existing files first just in case
write.exprs(eset.gcrma, file="gcrma_normal.txt")

pData(eset.gcrma) # Lists the analyzed file names.