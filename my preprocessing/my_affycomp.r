# A script that assess affy preprocessing
# Isidro C. Medina Jr., 2008

# Load the required libraries
library(affy); library(affycomp)

# Set the Working directory, the working directory should contain all the Affy CEL files
setwd=("E:\\My Thesis\\Training Set\\learning_set\\PROCESSED\\OUTPUT")
load("raw_os.RData")
load("gcrma_oc.RData")


# [1] Un-normalized data
png(file="raw_os_boxplot.png")      # Generates boxplot for un-normalized log intensity values.
boxplot(data, col="red", main="Raw OS Data",xlab="Samples", ylab="Intensities")
dev.off()
         
png(file="raw_os_histogram.png")    # Generates histogram for un-normalized log intensity values.
hist(data, col="red", main="Density Estimates of Raw OS Data",xlab="Samples", ylab="Intensities")
dev.off()

# [2] GC-RMA normalized data 
png(file="gcrma_os_boxplot.png")     # Generates boxplot for log intensity values.
boxplot(data.frame(exprs(eset.gcrma)), col="blue", main="GC-RMA Normalized OS Data", xlab="Samples", ylab="Intensities")
dev.off()

png(file="gcrma_os_histogram.png")    # Generates histogram for un-normalized log intensity values.
hist(exprs(eset.gcrma), col="blue", main="Density Estimates of GC-RMA Normalized OS Data", xlab="Samples", ylab="Intensities") 
dev.off()

