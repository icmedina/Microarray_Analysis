# A script that performs MAS5, RMA & GCRMA   preprocessing
# Isidro C. Medina Jr., 2008

#load the required libraries
library(gcrma); library(affy); library(matchprobes); library(hu6800); library(hu6800probe); library(hu6800.db);
library(hu6800cdf);library(hu6800subacdf); library(hu6800subbcdf); library(hu6800subccdf); library(hu6800subdcdf);

# Set the Working directory, the working directory should contain all the Affy CEL files
setwd = ("E:\\My Thesis\\Training Set\\learning_set\\PROCESSED")

# ReadAffy(), Reads all *.CEL (*.cel) files in your current working directory 
system.time({
  #targets <- readTargets("targets.txt")         # Import targets information. One entry per line.
  #data <- ReadAffy(filenames=targets$FileName)  # Import expression raw data stores them as AffyBatch object 'data'.
  data <- ReadAffy(celfile.path="E:\\My Thesis\\Training Set\\learning_set\\01 INPUT Cel Files\\Normal")
})

QCReport(data)

# Generates boxplot for un-normalized log intensity values.
png(file="OUTPUT//raw_OC_boxplot.png") 
boxplot(data, col="red", main="Raw Data"); dev.off()

# [1] MAS 5 Normalization
ptm <- proc.time()
  eset.mas5 <- mas5(data)
proc.time() - ptm	

# Generates boxplot for MAS5 normalized log intensity values.
png(file="OUTPUT//mas5_oc_boxplot.png")
boxplot(data.frame(exprs(eset.mas5)), col="blue", main="MAS 5 Normalized Data"); dev.off() 

# Export the MAS5 expression values to text file mydata.txt 
unlink("OUTPUT//mydata_oc_mas5.txt")   #delete existing files first just in case
write.exprs(eset.mas5, file="OUTPUT//mas5_oc_data.txt")

# [2] RMA Normalization
  ptm <- proc.time() 	
    eset.rma <- rma(data)
  proc.time() - ptm				
# Generates boxplot for RMA normalized log intensity values.
png(file="OUTPUT//rma_oc_boxplot.png")
boxplot(data.frame(exprs(eset.rma)), col="blue", main="RMA Normalized Data"); dev.off()

# Export the RMA expression values to text file mydata.txt (tab delimited format)
unlink("OUTPUT//rma_oc_data.txt")            # delete existing files first just in case 
write.exprs(eset.rma, file="OUTPUT//rma_oc_data.txt")
  
# [3] GC-RMA Normalization
ptm <- proc.time()                                    
  #my.affinity.info <- compute.affinities.local(data, verbose=TRUE)
  eset.gcrma <- gcrma(data,affinity.info=NULL)
proc.time() - ptm

# Generates boxplot for GC-RMA normalized log intensity values.
png(file="OUTPUT//gcrma_oc_boxplot.png")
boxplot(data.frame(exprs(eset.gcrma)), col="blue", main="GC-RMA Normalized Data"); dev.off()
				
#export the GC-RMA expression values to text file mydata.txt 
unlink("OUTPUT//gcrma_oc_data.txt")       #delete existing files first just in case
write.exprs(eset.gcrma, file="OUTPUT//gcrma_oc_data.txt")

pData(eset.gcrma) # Lists the analyzed file names.

# Compare RMA, Mas 5 & GC-RMA
  ptm <- proc.time() 
    comp.mat <- data.frame(
    MAS5=log2(exprs(eset.mas5)[,1]),
    RMA=exprs(eset.rma)[,1],
    GCRMA=exprs(eset.gcrma)[,1])
  proc.time() - ptm	 
plot(comp.mat,pch=16)