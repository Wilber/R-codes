library(limma)
library(Glimma)
library(edgeR)
#library(Mus.musculus)
library(Rsubread)
library(RColorBrewer)
#library('knitr')
library('limma')
library('reshape2')
#library('WGCNA')
library('gplots')
library(Rsubread)

#Adapted largely from: https://f1000research.com/articles/5-1408


#For phytopthera
setwd("/users/PAS0471/osu9657/RNA_Seq_Analysis/Data/chrX_data")
setwd("/users/PAS0107/osu6702/project/samples2")

#Step 1: Building an index
#Must provide a single FASTA file (eg. “genome.fa”)
#build index
buildindex(basename="PhytoRsubread_index",reference="genome/Phytophthora_infestans.ASM14294v1.30.dna.genome.fa", memory=128000)

#Load files: 
fastq.files.R1<-list.files(path = "samples2/", pattern = "SL.*_R1.fastq.gz$", full.names = TRUE)
fastq.files.R2<-list.files(path = "samples2/", pattern = "SL.*_R2.fastq.gz$", full.names = TRUE)


#Step 2: Aligning the reads
#Map paired-end reads:
align(index="PhytoRsubread_index",readfile1 = fastq.files.R1 ,readfile2 = fastq.files.R2 ,type = "rna", nthreads = 28)

#NOTE:Consider mapping with "Subjunc" function, enables exon-spanning and alternative splicing alignment 


#Check parameters used in alignment: 
args(align)

##Summary of proportion of read alignment: 
bam.files <- list.files(path = "samples", pattern = ".BAM$", full.names = TRUE)
props<-propmapped(bam.files, properlyPaired=TRUE)
write.table(props,"PhalignmentProportionsRsubread.txt", sep = "\t")


#Get bam files:
bam.files <- list.files(path = "samples2/", pattern = "SL.*.BAM$", full.names = TRUE) #Second 
bam.files <- list.files(path = "samples/", pattern = "SL.*.BAM$", full.names = TRUE)


##Remane to reduce the filenames, for aesthetics:
library(stringr)
bam.files<-word(bam.files, 6, sep = fixed('_')) #change the number, and the character  
bam.files<-word(bam.files, 1, sep = fixed('.'))

#Get feature counts 
fc <- featureCounts(bam.files, annot.ext = "genome/Phytophthora_infestans.ASM14294v1.30.gtf", 
                    isGTFAnnotationFile = TRUE, nthreads=28, isPairedEnd=TRUE, 
                    GTF.featureType = "gene")
# See what slots are stored in fc
names(fc)

## Take a look at the featurecounts stats
fc$stat
annotation<-(fc$annotation)

##Counts 
head(fc$counts)

#####DATA PRE-PROCESSING:
##Transformations from the raw-scale
##CPM normalizaton:
#raw counts are converted to CPM and log-CPM values using the cpm function
cpm <- cpm(fc)
lcpm <- cpm(fc, log=TRUE)

#####Removing genes that are lowly expressed
##OPTION 1:
# get # of genes in this dataset that have zero counts across all nine?? samples.
table(rowSums(fc$counts==0)==9) ##change 9 to the number of samples you have in your dataset
keep.exprs <- rowSums(cpm>1)>=3 #Filter by cpm > 1 in 3? or more samples 
fc <- fc[keep.exprs,, keep.lib.sizes=FALSE]
dim(fc)
##OPTION 2:
#Use the "filterByExpr" edgeR function. By default, the fxn keeps genes with about 10 read counts or more in a 
#minimum number of samples, where the number of samples is chosen according to the minimum group sample size. 
#The actual filtering uses CPM values rather than counts in order to avoid giving preference to samples with large library sizes. 
#For example, if the median library size is about 51 million and 10/51 ≈ 0.2, so the filterByExpr function keeps genes that have a CPM of 0.2 or
#more in at least three samples. 

keep.exprs <- filterByExpr(fc, group=group)
fc <- fc[keep.exprs,, keep.lib.sizes=FALSE]
dim(fc)



library(RColorBrewer)
nsamples <- ncol(fc)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
samplenames<-colnames(fc)
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(fc, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
     main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")


####Normalising gene expression distributions
fcUnnormalized<-fc

#TMM normalization
#Read more here: http://www.rna-seqblog.com/which-method-should-you-use-for-normalization-of-rna-seq-data/



lcpm <- cpm(fcUnnormalized, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="Unnormalised data",ylab="Log-cpm")

fc <- calcNormFactors(fc, method = "TMM")

lcpm <- cpm(fc, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="Normalised data",ylab="Log-cpm")



fc2 <- calcNormFactors(fcUnnormalized, method = "TMM")

lcpm2 <- cpm(fc2, log=TRUE)
boxplot(lcpm2, las=2, col=col, main="")
title(main="Normalized data",ylab="Log-cpm")

##OR IF NOT REMOVING THEN:
#fc2<-fc

####Unsupervised clustering of samples
lcpmUnnormalized <- cpm(fcUnnormalized, log=TRUE)


plotMDS(lcpmUnnormalized, top = 1000,labels=fc2$samples$strain, col=as.numeric(fc2$samples$strain)) #slight improvement in clustering 
plotMDS(lcpmUnnormalized, top = nrow(lcpmUnnormalized), labels=fc2$samples$strain, col=as.numeric(fc2$samples$strain))

##########
#OR
#PCA Adapted from: http://monashbioinformaticsplatform.github.io/RNAseq-DE-analysis-with-R/RNAseq_DE_analysis_with_R.html#install-and-load-packages
#########

# select data for the 1000 most highly expressed genes
select = order(rowMeans(fcUnnormalized$counts), decreasing=TRUE)[1:1000]
highexprgenes_counts <- fcUnnormalized$counts[select,]


# transpose the data to have variables (genes) as columns
data_for_PCA <- t(highexprgenes_counts)
dim(data_for_PCA)



## calculate MDS (matrix of dissimilarities)
mds <- cmdscale(dist(data_for_PCA), k=3, eig=TRUE)  
# k = the maximum dimension of the space which the data are to be represented in
# eig = indicates whether eigenvalues should be returned
mds$eig

##How many components can explain the variabilty?
# transform the Eigen values into percentage
eig_pc <- mds$eig * 100 / sum(mds$eig)
# plot the PCA
#png(file="~/PCA_PropExplainedVariance.png")
barplot(eig_pc,
        las=1,
        xlab="Dimensions", 
        ylab="Proportion of explained variance (%)", y.axis=NULL,
        col="darkgrey")
#dev.off()
## calculate MDS
mds <- cmdscale(dist(data_for_PCA)) # Performs MDS analysis 

#Samples representation
#png(file="~/PCA_Dim1vsDim2.png")
plot(mds[,1], -mds[,2], type="n", xlab="Dimension 1", ylab="Dimension 2", main="")
text(mds[,1], -mds[,2], rownames(mds), cex=0.8) #Pairwise comparisons between PAMA3qg and the rest???


#dev.off()

##FITTING MODELS TO DATA:

# create design matrix for differential expression analysis;
# if you wanted to account for batch here, you could simply include a batch
# term in the linear model at this step, e.g.:

mod <- model.matrix(~0+fc2$samples$strain )


# generate a list of all possible pairwise contrasts
#condition_pairs <- t(combn(levels(fc2$samples$strain), 2))

# generate a list of possible pairwise contrasts between PAMA3_qg and the rest:
#Below, a +ve log FC denotes gene upregulated in strain relative to PAMA3_qc baseline 
contr.matrix <- makeContrasts(
  PAMA15vsPAMA3_qg = PAMA15 - PAMA3_qg, 
  PAMA18vsPAMA3_qg = PAMA18 - PAMA3_qg, 
  PAN27vsPAMA3_qg = PAN27 - PAMA3_qg,
  PAN39vsPAMA3_qg = PAN39 - PAMA3_qg,
  PCA23vsPAMA3_qg = PCA23 - PAMA3_qg,
  PCA79vsPAMA3_qg = PCA79 - PAMA3_qg,
  PLL66vsPAMA3_qg = PLL66 - PAMA3_qg,
  PLL67vsPAMA3_qg = PLL67 - PAMA3_qg,
  PPU103vsPAMA3_qg = PPU103 - PAMA3_qg,
  PAMA3vsPAMA3_qg = PAMA3 - PAMA3_qg,
  levels = colnames(mod))

#Estimate dispersion: 
fc3 <- estimateDisp(fc2, design = mod )
plotBCV(fc3)


#Fit GLM QL:
fit <- glmQLFit(fc3, mod)

sampleNames<-sampleNames[-remove, ]
rownames(sampleNames)<-1:30
colnames(fc3$counts)<-sampleNames$sample

head(fit$coefficients) # the estimated values of the GLM coefficients for each gene
#QL F-tests instead of the more usual LRT as they give stricter error rate  ...
#...control by accounting for the uncertainty in dispersion estimation:
PAMA15vsPAMA3_qg.qlf <- glmQLFTest(fit, contrast=contr.matrix[,"PAMA15vsPAMA3_qg"])
PAMA18vsPAMA3_qg.qlf <- glmQLFTest(fit, contrast=contr.matrix[,"PAMA18vsPAMA3_qg"])
PAN27vsPAMA3_qg.qlf <- glmQLFTest(fit, contrast=contr.matrix[,"PAN27vsPAMA3_qg"])
PAN39vsPAMA3_qg.qlf <- glmQLFTest(fit, contrast=contr.matrix[,"PAN39vsPAMA3_qg"])
PCA23vsPAMA3_qg.qlf <- glmQLFTest(fit, contrast=contr.matrix[,"PCA23vsPAMA3_qg"])
PCA79vsPAMA3_qg.qlf <- glmQLFTest(fit, contrast=contr.matrix[,"PCA79vsPAMA3_qg"])
PLL66vsPAMA3_qg.qlf <- glmQLFTest(fit, contrast=contr.matrix[,"PLL66vsPAMA3_qg"])
PLL67vsPAMA3_qg.qlf <- glmQLFTest(fit, contrast=contr.matrix[,"PLL67vsPAMA3_qg"])
PPU103vsPAMA3_qg.qlf <- glmQLFTest(fit, contrast=contr.matrix[,"PPU103vsPAMA3_qg"])
PAMA3vsPAMA3_qg.qlf <- glmQLFTest(fit, contrast=contr.matrix[,"PAMA3vsPAMA3_qg"])

topTags(PAMA15vsPAMA3_qg.qlf)


DEgenesPAMA15<-topTags(PAMA15vsPAMA3_qg.qlf, n = Inf, adjust.method = "BH", sort.by = "PValue", p.value = 0.01)
summary(decideTestsDGE(PAMA15vsPAMA3_qg.qlf, adjust.method = "BH", p.value = 0.01, lfc = 1))
write.table(DEgenesPAMA15$table[abs(DEgenesPAMA15$table$logFC)>=1,], "DEgenesPAMA15.txt", sep = "\t", quote = FALSE)

DEgenesPAMA18<-topTags(PAMA18vsPAMA3_qg.qlf, n = Inf, adjust.method = "BH", sort.by = "PValue", p.value = 0.01)
summary(decideTestsDGE(PAMA18vsPAMA3_qg.qlf, adjust.method = "BH", p.value = 0.01, lfc = 1))
write.table(DEgenesPAMA18$table[abs(DEgenesPAMA18$table$logFC)>=1,], "DEgenesPAMA18.txt", sep = "\t", quote = FALSE)

DEgenesPAN27<-topTags(PAN27vsPAMA3_qg.qlf, n = Inf, adjust.method = "BH", sort.by = "PValue", p.value = 0.01)
summary(decideTestsDGE(PAN27vsPAMA3_qg.qlf, adjust.method = "BH", p.value = 0.01, lfc = 1 ))
write.table(DEgenesPAN27$table[abs(DEgenesPAN27$table$logFC)>=1,], "DEgenesPAN27.txt", sep = "\t", quote = FALSE)

DEgenesPAN39<-topTags(PAN39vsPAMA3_qg.qlf, n = Inf, adjust.method = "BH", sort.by = "PValue", p.value = 0.01)
summary(decideTestsDGE(PAN39vsPAMA3_qg.qlf, adjust.method = "BH", p.value = 0.01, lfc = 1))
write.table(DEgenesPAN39$table[abs(DEgenesPAN39$table$logFC)>=1,], "DEgenesPAN39.txt", sep = "\t", quote = FALSE)

DEgenesPCA23<-topTags(PCA23vsPAMA3_qg.qlf, n = Inf, adjust.method = "BH", sort.by = "PValue", p.value = 0.01)
summary(decideTestsDGE(PCA23vsPAMA3_qg.qlf, adjust.method = "BH", p.value = 0.01, lfc = 1))
write.table(DEgenesPCA23$table[abs(DEgenesPCA23$table$logFC)>=1,], "DEgenesPCA23.txt", sep = "\t", quote = FALSE)

DEgenesPCA79<-topTags(PCA79vsPAMA3_qg.qlf, n = Inf, adjust.method = "BH", sort.by = "PValue", p.value = 0.01)
summary(decideTestsDGE(PCA79vsPAMA3_qg.qlf, adjust.method = "BH", p.value = 0.01, lfc = 1))
write.table(DEgenesPCA79$table[abs(DEgenesPCA79$table$logFC)>=1,], "DEgenesPCA79.txt", sep = "\t", quote = FALSE)

DEgenesPLL66<-topTags(PLL66vsPAMA3_qg.qlf, n = Inf, adjust.method = "BH", sort.by = "PValue", p.value = 0.01)
summary(decideTestsDGE(PLL66vsPAMA3_qg.qlf, adjust.method = "BH", p.value = 0.01, lfc = 1))
write.table(DEgenesPLL66$table[abs(DEgenesPLL66$table$logFC)>=1,], "DEgenesPLL66.txt", sep = "\t", quote = FALSE)

DEgenesPLL67<-topTags(PLL67vsPAMA3_qg.qlf, n = Inf, adjust.method = "BH", sort.by = "PValue", p.value = 0.01)
summary(decideTestsDGE(PLL67vsPAMA3_qg.qlf,adjust.method = "BH", p.value = 0.01, lfc = 1))
write.table(DEgenesPLL67$table[abs(DEgenesPLL67$table$logFC)>=1,], "DEgenesPLL67.txt", sep = "\t", quote = FALSE)

DEgenesPPU103<-topTags(PPU103vsPAMA3_qg.qlf, n = Inf, adjust.method = "BH", sort.by = "PValue", p.value = 0.01)
summary(decideTestsDGE(PPU103vsPAMA3_qg.qlf, adjust.method = "BH", p.value = 0.01, lfc = 1))
write.table(DEgenesPPU103$table[abs(DEgenesPPU103$table$logFC)>=1,], "DEgenesPPU103.txt", sep = "\t", quote = FALSE)

DEgenesPAMA3<-topTags(PAMA3vsPAMA3_qg.qlf, n = Inf, adjust.method = "BH", sort.by = "PValue", p.value = 0.01)
summary(decideTestsDGE(PAMA3vsPAMA3_qg.qlf,adjust.method = "BH", p.value = 0.01, lfc = 1))
write.table(DEgenesPAMA3$table[abs(DEgenesPAMA3$table$logFC)>=1,], "DEgenesPAMA3.txt", sep = "\t", quote = FALSE)
#DEgenesPAMA3<-topTags(qlf, n = 400, adjust.method = "BH", sort.by = "logFC", p.value = 0.05)
#head(DEgenesPAMA15$table)
#summary(decideTestsDGE(PAMA15vsPAMA3_qg.qlf))


######Annotation of DE genes:
library("biomaRt")
#Load protist database
ensembl<-useMart(biomart = "protists_mart", host = "protists.ensembl.org")
#List datasets:
listDatasets(ensembl)
#Select dataset to use:
ensembl =  useDataset("pinfestans_eg_gene", mart = ensembl)
getBM(attributes = c("broad_p_infestans","description"), filters = "broad_p_infestans",values = rownames(DEgenesPAMA3$table), mart = ensembl)
getBM(attributes = c("broad_p_infestans","description"), filters = "broad_p_infestans",values = , mart = ensembl)

##For Arabidopsis:

##Adapted 
######Annotation of DE genes:
install.packages("rlang")
library("biomaRt")
listMarts() #list available databases

#Load protist database
ensembl <- useMart(biomart = "plants_mart",
                   dataset = "athaliana_eg_gene",
                   host = "plants.ensembl.org")
#List datasets:
listDatasets(ensembl)
#Select dataset to use:
#ensembl =  useDataset("pinfestans_eg_gene", mart = ensembl)
#getBM(attributes = c("broad_p_infestans","description"), filters = "broad_p_infestans",values = rownames(DEgenesPAMA3$table), mart = ensembl)
features.RA <- getBM(attributes = c("ensembl_gene_id","go_id",
                                "description"),
                 filters = c("ensembl_gene_id"),
                 values = rownames(toptags.RA),
                 mart = ensembl)
#examine result
head (features.RA)
#Remove blank entries
features.RA <- features.RA[features.RA$go_id != '',]
# convert from table format to list format
geneID2GO <- by(features.RA$go_id,
                features.RA$ensembl_gene_id,
                function(x) as.character(x))
#examine result
head (geneID2GO)

#create 'topGOdata' object by calling
go.obj = new("topGOdata", ontology='BP'
, allGenes = int.genes
, annot = annFUN.gene2GO
, gene2GO = geneID2GO)
        

