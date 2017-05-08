# Microarray analysis in Bioconductor/R
# The file is not executable. Just a collection of commands.
# written by
# Morten Mattingsdal (with good help from Josef Thingnes)
# Bioinformatic Core Facility
# Rikshospitalet-Radiumhospitalet HF
# morten.mattingsdal@medisin.uio.no
# PLEASE: be sure to use parameters/settings appropriate to your data
#
# Worked on:
# R		2.1.0
# Limma		1.9.6
# Maanova	0.98-5
# Resourcerer	1.1.2
# Annotate	1.15.6
# AnnBuilder    1.8.0
# biomaRt       1.4.0
# XML           0.97-8
# RMySql        0.5-6
# DBI           0.1-9

# Analysis in LIMMA.
# Requirements = gal file, genepix gpr files,design/contrast matrix

library(limma)
files <- dir(pattern="gpr")
RG <- read.maimages(files,"genepix",wt.fun=wtflags(0.1))
RG <- backgroundCorrect(RG, method="subtract")
RG$genes <- readGAL()
RG$printer <- getLayout(RG$genes)
mm<-function(x) if(sum(x<400)>length(x)/3) rep(NA,length(x)) else x
for(i in 1:dim(RG$R)[1]) RG$R[i,]<-mm(RG$R[i,])
for(i in 1:dim(RG$G)[1]) RG$G[i,]<-mm(RG$G[i,])
RG$G[(RG$G<200)|(RG$G>65000)]<-NA
RG$R[(RG$R<200)|(RG$R>65000)]<-NA
RG$R[(RG$R>65000)]<-NA
RG$G[(RG$G>65000)]<-NA
#RG$weights[RG$R <400]   <- 0.1
#RG$weights[RG$G <400]   <- 0.1
#RG$weights[RG$R >65000]   <- 0.1
#RG$weights[RG$G >65000]   <- 0.1
MA <- normalizeWithinArrays(RG,RG$printer,method="loess",span=0.3,iterations=4)
MA <- normalizeBetweenArrays(MA,method="scale")
fit <- lmFit(MA, design=c(1,1))
fit <- eBayes(fit)
M <- fit$coef
A <- fit$Amean
ord <- order(fit$lods,decreasing=TRUE)
top30 <- ord[1:30]

x11()
pdf(file="MA_plot.pdf")
plot(A,M,pch=16,main="MA plot, lowess normalized, with scaling",cex=0.1)
text(A[top30],M[top30],labels=MA$genes$geneName[top30,"Name"],cex=0.8,col="blue")
abline(0,0,col="red")
dev.off()
x11()
pdf(file="Boxplot.pdf")
boxplot(MA$M~col(MA$M),names=colnames(MA$M),main="After median normalization and with scaling")
text(A[top30],M[top30],labels=MA$genes$geneName[top30,"Name"],cex=0.8,col="blue")
dev.off()

plot(fit$coef,fit$lods,pch=16,cex=0.2,xlab="Log Fold Change",ylab="Log Odds",main="default cutoffs")
abline(0,0,col="red")
abline(v=1,col="red")
abline(v=-1,col="red")

plot(A,M,pch=16,main="MA plot, loess normalized, without scaling",cex=0.1)
top <- topTable(fit,number=100,sort.by="B",adjust="fdr")
x <- topTable(fit, adjust="fdr", sort.by="M", number=2000)
y <- x[x$P.Value < 1 & (x$M > 0.5 | x$M < -0.5) & x$A > 9,];y; print("Number of genes in this list:"); length(y$ID)  # Same as above but with complex filter: P-value < 0.01 AND at least 2-fold change AND expression value A > 10.
write.table(top,file="top100B.txt",quote=FALSE,row.names=FALSE,sep="\t")


# Analysis in LIMMA based on referance design experiments
# Requirements = same as above, in addition design file

targets <- readTargets("design.txt")
design <- modelMatrix(targets,ref="ref")
MA <- normalizeWithinArrays(RG,RG$printer,method="loess",span=0.3,iterations=4)
MA <- normalizeBetweenArrays(MA,method="scale")
fit <- lmFit(MA, design)
contrast.matrix <- makeContrasts(diffv-diff,prov-pro,levels=design)		#define your contrasts
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
top <- topTable(fit2,number=100,sort.by="B",adjust="fdr")

# Import a txt file of ratios/log2 intenseties and build a linear model usin LIMMA
library(limma)
x<-as.matrix(read.table("genelist.txt",header=T,sep="\t",row.names=1))
design <- model.matrix(~ -1+factor(c(1,1,2,2,3,3,1,1,2,3,3,1,1,2,2,3,3)))
colnames(design) <- c("Ustim", "Asp", "LPS")
fit <- lmFit(x, design)
contrast.matrix <- makeContrasts(Asp-Ustim, LPS-Ustim, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
plot(fit2$coef[,1],fit2$lods[,1],xlab="Log2 Fold Change",ylab="Log Odds",pch=16,cex=0.2,main="test")
abline(0,0,col="red")
abline(v=-1,col="red")
abline(v=1,col="red")
asp <- topTable(fit2,coef=1, adjust="fdr", sort.by="B", number=5000)
test <- asp[asp$B >0 & (asp$M > 1 | asp$M < -1),]
write.table(test,"test.txt",sep="\t")

#Intensity dependent and IQR (interquartile range) filter for an eset data object. Thanks to Raffaele Calogero for providing the code

intensityFilter<-function(x, fraction,  threshold){
                     require(affy)
                     require(genefilter)
                     f1<-pOverA(fraction, threshold)
                     ff <- filterfun(f1)
                     which <- genefilter(x, ff)
                     par(mfrow=c(1,2))
                     text<-paste("Unfiltered probes= ", dim(exprs(x))[1], sep="")
                     hist(exprs(x), main=text, breaks=100)
                     text<-paste("Filtered probes= ", sum(which), sep="")
                     hist(exprs(x[which,]), main=text, breaks=100)
                     return(x[which,])
}

eset_intensity<-intensityFilter(eset, 0.33, log2(100))


# x is an exprSet, threshold is the IQR value to be applied
iqrFilter<-function(x, threshold){
                     require(affy)
                     require(genefilter)
                     iqr <- function(x){IQR(x)>threshold}
                     ff <- filterfun(iqr)
                     which <- genefilter(x, ff)
                     par(mfrow=c(1,2))
                     text<-paste("Unfiltered probes= ", dim(exprs(x))[1], sep="")
                     hist(exprs(x), main=text, breaks=100)
                     text<-paste("Filtered probes= ", sum(which), sep="")
                     hist(exprs(x[which,]), main=text, breaks=100)
                     return(x[which,])
}

eset_filtered<-iqrFilter(eset_intensity, 0.25)

# Present / Absent calls filter for Affymetrix data.
# Filtering is done on the RMA normalized eset object, using info from the mas5calls function from the affy package

data=ReadAffy(filenames=targets$FileNames)
PMA_calls=mas5calls(data)
eset=rma(data)
mascallsfilter <- function(cutoff = "p", number){function(x){sum(tolower(x) == tolower(cutoff)) >= number}}
f1 <- mascallsfilter("p", 6)							# minimum "P" in 6 arrays
filt <- filterfun(f1)
filter<-genefilter(PMA_calls,filt)
filtered_eset=eset[filter,]




########################################################################
# Venn Diagrams in LIMMA. Requires two lists of genes/probes
list1 <-x[,1]
list2 <-x[,2]
dooby.all <- union(list1,list2)
dooby.mat <- matrix(0,nrow=length(dooby.all),ncol=2)
colnames(dooby.mat) <- c("Aspergillus","LPS")
for (i in 1:length(dooby.all)) {
dooby.mat[i,1] <- dooby.all[i]%in%list1
dooby.mat[i,2] <- dooby.all[i]%in%list2}
dooby.ven <- vennCounts(dooby.mat)
vennDiagram(dooby.ven)


########################################################################
# Analysis in SAM (R inplementation of SAM)
# Requirements: valid LIMMA object MA

library(samr)
y<-paste(c("1Time1Start","1Time2","1Time3End","1Time1Start","1Time2","1Time3End","2Time1Start","2Time2","2Time3End","2Time1Start","2Time2","2Time3End"))
SAM<-list(x=MA$M,y=y,logged2=TRUE,geneid=paste(as.character(MA$genes$ID)),genenames=paste(as.character(MA$genes$GeneName),sep=""))
sam_r <- samr(SAM,resp.type="Two class unpaired timecourse",s0.perc=90, nperms=400)
delta.table <- samr.compute.delta.table(sam_r)
delta=0.7
x11()
pdf(file="SAM_plot.pdf")
samr.plot(sam_r,delta)
dev.off()
siggenes.table<-samr.compute.siggenes.table(sam_r,delta, SAM, delta.table)

#########################################################################
# Missing value estimation by kmeans nearest neighbour

library(EMV)
knn(RG$R,k=max(dim(m)[1]*0.01,2),na.rm=TRUE,correlation=FALSE, dist.bound=FALSE)
knn(RG$G,k=max(dim(m)[1]*0.01,2),na.rm=TRUE,correlation=FALSE, dist.bound=FALSE)

#########################################################################
# Analysis in MAANOVA
# Requirements = tab file of intensity in every channel(Big.txt)
# name,id,cy5(array1),cy3(array1),cy5(array2),cy3(array2).....
# description file,accession file,designfile

library(maanova)
x<-read.madata("Big.txt", designfile="design.txt",header=T,spotflag=T,pmt=3)
logx<-createData(x,log.trans=TRUE)
lowessx<-transform.madata(logx, method=c("glowess"),draw=c("off"))
model.maanova <- makeModel(data=lowessx, formula=~Dye+Array+Sample+vegf+Celltype) #define your formula based on design file
anova.maanova <- fitmaanova(lowessx, model.maanova)
varplot(anova.maanova)
resiplot(logx,anova.maanova)
test1000<-matest(lowessx,model.maanova,term="vegf", n.perm=1000)
testtmp1<-adjPval(test1000)
idx <- volcano(test1000, method=c(rep("unadj",4)),threshold=c(0.01,0.05,0.05,0.05),
highlight.flag=F,title="Volcanoplot with unadjusted P values after 1000 permutations.")$idx.Fs
ratio<-as.data.frame((anova.maanova$Sample[,1])-(anova.maanova$Sample[,2]))
pidx<-as.data.frame(test1000$Fs$Ptab)
dataidx<-as.data.frame(idx)
Ratio<-as.data.frame(ratio)
description<-read.table("description.txt", header=TRUE)
ACC<-read.table("ACC.txt", header=TRUE)
for (i in dataidx){
	Description<-paste(description[i,])
	Accession<-paste(ACC[i,])
	Ratio<-paste(ratio[i,])
	P<-paste(pidx[i,])
}
BIG<-data.frame(title1=Accession, title2=Description, title3=P, title4=Ratio)
ordBIG<-BIG[with(BIG,order(P, decreasing=FALSE)),]
write.table(ordBIG, "toplist_observed_p_1000perm.txt", row.names=FALSE, sep="\t")

#########################################################################
# Resourcerer builds an annotation package from available
# arrays from ftp://ftp.tigr.org/pub/data/tgi/Resourcerer

library(Resourcerer)
resourcerer2BioC("Agilent_MouseOligo.zip",organism="mouse",destDir=file.path(.path.package("Resourcerer"),"temp"),
pkgName=("AgilentMouseOligo"), srcUrls=getSrcUrl("all","Mus musculus"),
pkgPath=file.path(.path.package("Resourcerer"),"temp"),otherSrc=NULL,baseMapType="gbNRef",
version="1.0.0",fromWeb=True,baseUrl="ftp://ftp.tigr.org/pub/data/tgi/Resourcerer",
check=TRUE,author=list(author="Anonymous",maintainer="<morten.mattingsdal@medisin.uio.no>"))

# Or load in a gene list with probe ID & ref|gb and build annotation
# Working with R 2.4.0 and AnnBuilder 1.13.6

setwd("C:/hjemme/idag")
probe<-read.delim("new_operon.txt", header=F,sep="\t")
P<-as.matrix(probe)
pkgPath <- "C:/hjemme/idag"
write.table(P, file = file.path(myDir, "query"),sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
srcUrls=getSrcUrl('ALL',organism="Homo sapiens",xml=TRUE)
pkgName="OPERON"
organism="Homo sapiens"
version="1.0"
author=list(author="Anonymous",maintainer="<morten.mattingsdal@medisin.uio.no>")
baseMapType="ug"
ABPkgBuilder(baseName=file.path(myDir, "query"), srcUrls, baseMapType, otherSrc = NULL, pkgName, pkgPath, organism,version, author=list(author="Anonymous",maintainer="<morten.mattingsdal@medisin.uio.no>"), fromWeb = TRUE, lazyLoad = TRUE)

# Build several ogranism spesific datapackages based on sequence homology from Homologene with: (only replace author)
# the ID of Homologene seems to be wrong. Await a new release of AnnBuilder ?
library(AnnBuilder)
pkgPath<-tempdir()
version<-"1.1.0"
author<-list(author="Anonymous",maintainer="<morten.mattingsdal@medisin.uio.no>")
homoPkgBuilder(suffix = "homology", pkgPath, version, author, url =getSrcUrl("HG"))

# convert a list of pig acc nr Homologene ID
sschomology()
xx<-as.list(sschomologyACC2HGID)
pig<-read.delim("pigacc.txt", header=F,sep="\t")
PIG<-as.matrix(pig)
test<-intersect(names(xx),PIG)
jens=1
for (i in test){
ult[jens]<-paste(xx[i])
jens=jens+1
}

#########################################################################
# Builds a html table of LIMMA toplist. Requires a
# BioC metadata package. See Resourcerer or AnnBuilder
# Requirements = toplist and annotation library

library(annotate)
library(AgilentHuman1av2)
top<-read.delim("top100B.txt", header=T,blank.lines.skip = TRUE)		#read tab top list
top										#limma top object
probe<-as.matrix(top[,5])
description<-as.list(top[,10])
m<-top[,11]
a<-top[,12]
p<-top[,14]
b<-top[,15]
M<-round(m,2)
P<-round(p,5)
A<-round(a,2)
B<-round(b,2)
xx<-mget(probe,ifnotfound="NA",AgilentHuman1av2ACCNUM)
acc<-as.data.frame(cbind(xx))
ll <- getLL(probe,data="AgilentHuman1av2")
symbol <- getSYMBOL(probe,data="AgilentHuman1av2")
res<-data.frame(probe,cbind(unlist(symbol),unlist(acc),unlist(description),unlist(M),unlist(A),unlist(P),unlist(B)))
names(res)<-cbind("Probe","Symbol","Acc","Description","log2 Ratio","log2 Intensity","p value","B")
htmlpage(ll,filename="Results.html",title="Top 100 genes ranked by B(lods)",othernames=res,
table.head=c(("Entrez ID"),names(res)),table.center=TRUE)

#########################################################################
# Buildes a html page of LIMMA toplist of Affymetrix data

probe<-top[,1]
description<-aafDescription(as.character(probe),"hgu133plus2")
DE<-getText(description)
m<-top[,2]
p<-top[,5]
p_adj=top[,6]
b<-top[,7]
M<-round(m,2)
P<-round(p,5)
P_adj=round(p_adj,2)
B<-round(b,5)
xx<-mget(probe,ifnotfound="NA",hgu133plus2ACCNUM)
acc<-as.data.frame(cbind(xx))
ll <- getLL(probe,data="hgu133plus2")
symbol <- getSYMBOL(probe,data="hgu133plus2")
res<-data.frame(cbind(unlist(probe)),cbind(unlist(symbol),unlist(acc),unlist(DE),unlist(M),unlist(P),unlist(P_adj),unlist(B)))
names(res)<-cbind("Probe","Symbol","Acc","Description","log2 Ratio","p value","adjusted P value","B")
htmlpage(ll,filename="Results_lps_choe.html",title="Positive - Negative, sorted by M>1",othernames=res,
table.head=c(("Entrez ID"),names(res)),table.center=TRUE)

#########################################################################
# Given a file of Probe ID the following commands extracts Entrez ID
# from an annotation package and writes a file

library(mgug4121a)
id<-read.table("probeID.txt",header=T)
probe<-as.matrix(id)
ll<-getLL(probe,data="mgug4121a")
out<-as.data.frame(ll)
write.table(out,file="probe2ID.txt",quote=FALSE,row.names=FALSE)

#########################################################################
# GOstats
# Functions to explore GO ontologies with a file containing  probeIDs
# below are "Biological Process" with significance below 0.01
#
# Recuirements: file with differentially expressed genes (probe ID
# + annotation package for speficic arrays

library(GOstats)
id<-read.table("top100.txt",header=T)
probe<-as.matrix(id)
ll<-getLL(probe,data="AgilentHuman1av2")
go <- makeGOGraph(as.character(ll), "BP", removeRoot=FALSE)
hyp <- GOHyperG(unique(ll), lib="AgilentHuman1av2", what="BP")
gopval <- hyp$pvalues[nodes(go)]
gopval<- sort(gopval)
sig <- gopval[gopval < 0.01]
counts <- hyp$goCounts[names(sig)]
terms <- getGOTerm(names(sig))[["BP"]]
nch <- nchar(unlist(terms))
siggo<-matrix(c(names(terms), terms, round(sig, 3), counts),ncol=4, dimnames=list(1:length(sig),c("GO ID", "Term", "p-value", "# Genes")))
write.table(siggo,file="significant_BP_20t.txt",quote=FALSE,row.names=FALSE,sep="\t")

# Retrieve all genes on array with GO annotation "apoptosis (GO:0006915) and extract ratio values from the limma object "MA"

response to stress		GO:0006950
apoptosis   			GO:0006915
muscle development		GO:0007517
reponse to wounding		GO:0009611
development			GO:0007275
cell differention		GO:0030154
myoblast differentiation	GO:0045445
cell proliferation   		GO:0008283

md<-get("GO:0007517", GOBPCHILDREN)
myGO<-append(md,"GO:0007517")
for (i in myGO){
 probe<-unique(lookUp(i, "mgug4121a", "GO2ALLPROBES"))
}
res<-NULL
for (i in probe){
 res<-append(res,grep(i,MA$gene$ID))
}
myM<-matrix(nrow=(length(res)),ncol=12)
jens = 1
for(i in res){
 myM[jens,] = MA$M[i,]
 jens <- jens+1
}
name=NULL
for (i in res){
name<-append(name,MA$gene$GeneName[i])
jens<-jens+1
}
row.names(myM)<-name
write.table(myM,file="muscle development.txt",sep="\t",row.names=T)


### for Affymetrix raw data (only rma correction
library(annotate)
library(hgu133plus2)
library(affy)
library(limma)
library(convert)
library(GO)
targets <- readTargets("targets.txt", sep="")
Data <- ReadAffy(filenames=targets$filename)
RMAeset<-rma(Data)
exprs2excel(RMAeset, file="RMAwiggo.csv")
wiggo<-read.delim("RMAwiggo.txt", header=T,blank.lines.skip = TRUE)
wiggo1<-as.matrix(wiggo)
md<-get("GO:0006955", GOBPCHILDREN)
myGO<-append(md,"GO:0006955")
for (i in myGO){
 probe<-unique(lookUp(i, "hgu133plus2", "GO2ALLPROBES"))
}
res<-NULL
for (i in probe){
 res<-append(res,grep(i,wiggo1[,1]))
}
myM<-matrix(nrow=(length(res)),ncol=6)
jens = 1
for(i in res){
 myM[jens,] = wiggo1[i,]
 jens <- jens+1
}
write.table(myM,file="immune response.txt",sep="\t",row.names=T)


#Affymetrix analysis with "Preferred analysis methods for Affymetrix GeneChips revealed by a wholly defined control dataset"
data <- ReadAffy(filenames=targets$FileName)
eset <- expresso(data, bgcorrect.method="mas",normalize.method="gcrma",pmcorrect.method="mas",summary.method="medianpolish")
#do second loess normalization on eset object and export data.frame




#########################################################################
# Plots of different GO graphs. Requires Unix & Rgraphviz
# Uses objects from above

# Plot of GO terms
att <- list()
lab <- rep("", length(nodes(go)))
names(lab) <- nodes(go)
att$label <- lab
plot(go, nodeAttrs=att)

#Plot of GO terms with colored nodes
col <- ifelse(gopval < 0.01, ifelse(counts >= 5, "blue", "orange"), "white")
names(col) <- names(gopval)
att$fillcolor <- col
plot(go, nodeAttrs=att)

#########################################################################
# biomaRt functions. To extract information from vega / ensembl


# Very useful Biomart functions
library(biomaRt)
mart <- martConnect( biomarts = c("ensembl","vega"))
     #list available AffyArrays at mart
     getAffyArrays(mart = mart)
     #example using affy id from Ensembl
     g = getGene( id = "1939_at", array = "hg_u95av2", mart = mart)
     #example using locuslink id from vega
     g = getGene( id = 100, type = "entrezgene", species = "hsapiens", db = "vega", mart = mart)
     #list availabe species from ensembl & vega
     species =getSpecies(mart = mart)
     show(species)
     #retrives SNPs
     mart <- martConnect(biomarts=c("ensembl","snp"))
     snp = getSNP(chromosome = 1, start = 30065068, end = 30104323, species ="hsapiens", mart = mart)
     #OMIM
     omim = getOMIM( id = 672, type = "entrezgene", mart = mart) 		# unigene
     omim = getOMIM( id = "1939_at", array = "hg_u95av2", mart = mart)		# affy id
     # get Homologs
     #HUGO to Entrez Gene
     homolog = getHomolog(id = 1:20, from.species = 'hsapiens', to.species = 'mmusculus', from.type = 'entrezgene', to.type = 'refseq', mart = mart)
     show(homolog)
     #ensembl to ensembl
     homolog = getHomolog( id = "ENSG00000072778", from.species = "hsapiens", from.type = "ensembl", to.type="ensembl",to.species="mmusculus", mart = mart )
     show(homolog)
     #Affy to Affy
     homolog = getHomolog( id = "1424184_at", to.array = "hg_u95av2", from.array = "mouse430a_2", mart = mart )
     show(homolog)
     #Affy to ensembl
     homolog = getHomolog( id = "1424184_at", to.species = "hsapiens", to.type = "ensembl", from.array = "mouse430a_2", mart = mart )
     show(homolog)
     #example using affy id
     go = getGO( id = "1939_at", array = "hg_u95av2", mart = mart)
     #example using locuslink id
     go = getGO( id = 672, type = "entrezgene", species = "hsapiens", mart = mart)
     # get FASTA sequences
     mart <- martConnect()
     martConnect(biomarts='sequence',mart=mart,host='ensembldb.ensembl.org',user='anonymous',password='')
     seq<-getSequence(chromosome=c(2,2),start=c(100000,30000),end=c(100300,30500),species="hsapiens",mart=mart)
     exportFASTA(seq,file="test.fasta")

martDisconnect(mart = mart)

########################################################################
# Various clustering methods

# Heatmap
library(limma)
library(marray)
x<-read.table("matrix.txt",header=T,row.names=1)			#matrix with genes & log ratio values
y<-as.matrix(x)
my.pal <- maPalette(low="white", high="blue")
myhc <- function(n) hclust(n, method = "complete")
mydist <- function(m) dist(m, method = "euclidian")
heatmap(y,distfun = mydist,hclustfun = myhc,Rowv=NA,scale="none",col=my.pal,labRow=FALSE,cexCol=0.5)

#Hierarchical clustering
x<-read.table("matrix.txt",header=T,row.names=1)
distance<-dist(x,method="euclidian")
cluster<-hclust(distance)
par(cex=0.7)
plot(cluster)

#kmeans clustering
x<-read.table("matrix.txt",header=T,row.names=1)
y<-as.matrix(x)
kmean<-kmeans(y,5,iter.max=100,algorithm = c("Forgy"))
plot(y, col = kmean$cluster)
points(kmean$centers, col = 1:5, pch = 8)
text(y,labels=(row.names=y),cex=0.6,col="black")


########################################################################
# Some plotting options

plot.ts(Ratio,plot.type="single",col=4:2,xy.lines=T,axes=T,xlab="gene nr")
text(x1,labels=n,cex=0.6,col="black")
legend(.4,-1, "9 timer", pch=21, pt.bg="white", lty=1,bty="n", col = "blue")
legend(.4,-1.1, "24 timer", pch=21, pt.bg="white", lty=1, bty="n", col = "green")
legend(.4,-1.2, "48 timer", pch=21, pt.bg="white", lty=1, bty="n",col = "red")



data<-matrix(ncol=3,nrow=3)
data[,1]<-c(6780,6308,7402)
data[,2]<-c(4452,3919,3860)
data[,3]<-c(7574,7038,8691)
row.names(data)<-c("Basal level","Abstention period","Red wine period")
colnames(data)<-c("Leptin (all)","Leptin (males)","Leptin (females)")
par(cex=1.2,yaxt="l")
plot.ts(data,plot.type="single",col=1:3,xy.lines=T,ylab="pg/ml",lwd=4,xlab="Abstention period                       Red wine period")
legend(1,8600, "Leptin, all", pch=0, pt.bg="white", lwd=4,bty="n",pt.cex=FALSE, col = "black")
legend(1,8800, "Leptin, females", pch=21, pt.bg="white", bty="n",lwd=4,pt.cex=FALSE, col = "green")
legend(1,8400, "Leptin, males", pch=21, pt.bg="white", bty="n",lwd=4,pt.cex=FALSE, col = "red")


#########################################################################
# Basic text miner in R/CRAN. Count co-occuring terms in MedLine/PubMed
# Adopted from MedLineR


library(XML)
options("serviceUrl.entrez" = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/")

countApair<- function (
    term1, term2,
    termAdditional="",
    baseUrl=getOption("serviceUrl.entrez")
    ) {

    # QC: make sure the baseUrl is all right.
    if (is.null (baseUrl)) {
        stop ("Need to define the URL of the Pubmed service!")
    }
    query<- paste (baseUrl,
                   "esearch.fcgi?",
                   "db=pubmed&",					# replace pubmed with db of interest (omim,snp,unigene ++)
                   "rettype=count&",
                   "term=",
                     term1, "+AND+", term2, termAdditional,
                   sep="")
    result.xml<- try (xmlTreeParse(file=query, isURL=T))
    count<- as.numeric(xmlValue (xmlRoot (result.xml) [["Count"]]))
    return (count)
}

pauseBetweenQueries<- function (
   sleep.peak=15,                # pause (in seconds) during peak hours
   sleep.offpeak=3               # pause (in seconds) during off-peak
  ) {
 result.date<- unlist (strsplit(
   date(), split=" "))
 hour<- as.numeric(unlist (strsplit (result.date[4], split=':'))[1])
 if (
   (result.date[1]=="Sat") | (result.date[1]=="Sun") |
   (hour > 21) | (hour<5)
  ) {off.peak<-T} else {off.peak<-F}
 if (off.peak) {
 print("--Off hours at NCBI (faster)--")
  Sys.sleep (sleep.offpeak)
 } else {
 print("--Its peaktime at NCBI (slower)--")
  Sys.sleep (sleep.peak)
 }
}
termList=c("alcohol","benefit","gene","chromosome","income","norway")	# Define your searchterms
n.terms<- length (termList)
matrix (0, ncol=n.terms, nrow=n.terms)
n.terms<- length (termList)
result.matrix<-matrix (0, ncol=n.terms, nrow=n.terms)
 for (i in 1:n.terms){
  result.matrix[i,i]<- countApair (
    term1=termList[i],
    term2=termList[i])
  pauseBetweenQueries()
 }

for (i in 1:(n.terms-1)){
if (result.matrix [i,i]==0) {next}
 for (j in (i+1):n.terms) {
     n.counts <-countApair (
         term1=termList[i],
         term2=termList[j])
     pauseBetweenQueries()
     result.matrix[i,j]<- n.counts
     result.matrix[j,i]<- n.counts
   }
 }
row.names(result.matrix)<-termList
colnames(result.matrix)<-termList
result.matrix
dotchart(result.matrix,cex=0.7)

#########################################################################
# GEOquery. Retriveve public data from GEO, examples are from affymetrix & 2-channel data

GDS90<- getGEO("GDS90", destdir=".")		#Perou CM, PNAS, 1999
Meta(GDS90)$channel_count
MA<- GDS2MA(GDS90)				#2 channel data (creates a MA LIMMA object)

GDS360<- getGEO("GDS360", destdir=".")		#Chang, Lancet, 2003
Meta(GDS360)$channel_count
eset<-GDS2eSet(GDS360)				#Affymetrix data (creates an eset Affy object)




#########################################################################
# Power calculation for microarray data in Bioconductor/R (both dual and single channel)
# The package ssize requires a control experiment variable. The data below is derived from Affymetrix arrays
# NB package financed by Pfizer

library(ssize)
library(affy)
library(limma)
targets <- readTargets("Targets.txt")
data <- ReadAffy(filenames=targets$FileName)

eset <- expresso(data, normalize.method="quantiles", bg.correct=TRUE,bgcorrect.method="mas",
3,pmcorrect.method="mas", summary.method="medianpolish")# normalization according to choe et al., 2005

eset<-normalize.loess(eset,log.it=FALSE) 		# reccomended second loess by choe et al., 2005
eset=as.data.frame(exprs(eset))
design <- c(1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3)
controls <- eset[,design=="3"]
exp.sd <- apply(controls, 1, sd)			# make your matrix with control arrays

n <- 6							# parameters for sstats. The sig.level is Bonferonni adjusted below
fold.change <- 2
power <- 0.80
sig.level <- 0.05

hist(exp.sd, n = 20, col = "cyan", border = "blue", main = "",xlab = "Standard Deviation (for data on the log2 scale)")
dens <- density(exp.sd)
lines(dens$x, dens$y * par("usr")[4]/max(dens$y), col = "red",lwd = 2)
title("Histogram of Standard Deviations")

# the following commands asnswers: "What is the power for 6 patients pr group with log2(ratio)>1 and p<0.05"
all.power<-pow(sd=exp.sd, n=n, delta=log2(fold.change), sig.level=sig.level, alpha.correct = "Bonferonni")
power.plot(all.power, lwd = 2, col = "blue")
xmax <- par("usr")[2] - 0.05
ymax <- par("usr")[4] - 0.05
legend(x = xmax, y = ymax, legend = strsplit(paste("n=", n, ",","fold change=", fold.change, ",", "alpha=", sig.level, ",","# genes=", nobs(exp.sd), sep = ""), ",")[[1]], xjust = 1,yjust = 1, cex = 1)
title("Power to Detect 2-Fold Change")

# the following commands asnswers: "What is the size pr group to archive 80% power, log2(ratio)>1 ans p<0.05"
all.size <- ssize(sd = exp.sd, delta = log2(fold.change), sig.level = sig.level,power = power,alpha.correct = "Bonferonni")
ssize.plot(all.size, lwd = 2, col = "magenta", xlim = c(1, 20))
xmax <- par("usr")[2] - 1
ymin <- par("usr")[3] + 0.05
legend(x = xmax, y = ymin, legend = strsplit(paste("fold change=",fold.change, ",", "alpha=", sig.level, ",", "power=", power,",", "# genes=", nobs(exp.sd), sep = ""), ",")[[1]], xjust = 1,yjust = 0, cex = 1)
title("Sample Size to Detect 2-Fold Change")

# the following commands asnswers: "What is necessary fold change to achive 80% power with 3 patients pr group when log2(ratio)>1 and p<0.05"
all.delta <- delta(sd = exp.sd, power = power, n = n, sig.level = sig.level,alpha.correct = "Bonferonni")
delta.plot(all.delta, lwd = 2, col = "magenta", xlim = c(1, 10))
title("Fold Change to Achieve 80% Power")
xmax <- par("usr")[2] - 1
ymin <- par("usr")[3] + 0.05
legend(x = xmax, y = ymin, legend = strsplit(paste("n=", n, ",","alpha=", sig.level, ",", "power=", power, ",", "# genes=",nobs(exp.sd), sep = ""), ",")[[1]], xjust = 1, yjust = 0,cex = 1)
title("Fold Change to Achieve 80% Power")


# The following code preforms the exploratory "stam" algorithm
library(stam)
library(limma)
library(affy)
library(gcrma)
setwd("C:/prosjekter/trine/cel")
data=ReadAffy(widget=T)
eset=gcrma(data)
golubNorm.eval.explore <-
stam.evaluate(eset, "Covar1", testset = NULL,
chip="hgu133plus2", root="GO:0044237",
alpha=seq(0, 1, 0.1), ndelta=5)

# The following code does significant pathway analysis, accoring to the sigPathway package. "G" is the pathway gene sets
library(sigPathway)
data(MuscleExample)
phenotype=c("vector","vector","vector","vector","vector","PCSK9","PCSK9","PCSK9","PCSK9","PCSK9")
sig_pathways <- runSigPathway(G, 10, 500, new, phenotype, nsim = 1000,weightType = "constant", ngroups = 2, npath = 25, verbose = FALSE,allpathways = TRUE, annotpkg = "hgu133plus2", alwaysUseRandomPerm = FALSE)
writeSigPathway(sig_pathways)

#END