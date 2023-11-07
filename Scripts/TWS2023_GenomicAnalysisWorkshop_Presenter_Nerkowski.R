##### TWS 2023 Genomic Analysis Workshop #####

##Set working directory
##Use the path to where you downloaded the data ("batch_1.vcf")
setwd("~/Desktop/MEWG/TWS R Workshop 2023/")

## Load Packages ##
install.packages("devtools")
install.packages("adegenet")
install.packages("vcfR")
install.packages("tidyverse")
install.packages("hierfstat")
install.packages("ggplot2")
install.packages("pcadapt") #this one might need to be downloaded from github (instructions below)
#website for pcadapt https://bcm-uga.github.io/pcadapt/articles/pcadapt.html

## Install LEA ##
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("LEA")

## Install ggman and qvalue ##
library(devtools)
install_github("drveera/ggman")
install_github("jdstorey/qvalue")
#here is an alternative way to get pcadapt if the CRAN method doesn't work
install_github("bcm-uga/pcadapt")
devtools::install_github("pievos101/PopGenome")


#Load Packages ##
library(LEA)
library(vcfR)
library(adegenet)
library(tidyverse)
library(hierfstat)
library(ggman)
library(ggplot2)
library(pcadapt)
library(qvalue)
library(PopGenome)

##Load SambaR 
source("https://github.com/mennodejong1986/SambaR/raw/master/SAMBAR_v1.09.txt")
getpackages()
# This will take some time and install many different genomic packaeges that are useful
# You may get some prompts to "install from sources the package which needs compilation" messages
# I typically just type 'yes'
# If you run into errors, please see the manual and the github page
# https://github.com/mennodejong1986/SambaR/tree/master
# There is a ton of useful information there!

## Load prairie dog dataset ##
vcf <- read.vcfR("batch_1.vcf")
vcf

# ***** Object of Class vcfR *****
#   237 samples
# 1826 CHROMs
# 3,549 variants
# Object size: 8.2 Mb
# 0 percent missing data
# *****        *****         *****

vcf@meta

vcf@fix

vcf@gt

pdogGen <- vcfR2genlight(vcf)
pdogGen
#' /// GENLIGHT OBJECT /////////
#'   
#'   // 237 genotypes,  3,549 binary SNPs, size: 1.1 Mb
#' 39695 (4.72 %) missing data
#' 
#' // Basic content
#' @gen: list of 237 SNPbin
#' 
#' // Optional content
#' @ind.names:  237 individual labels
#' @loc.names:  3549 locus labels
#' @chromosome: factor storing chromosomes of the SNPs
#' @position: integer storing positions of the SNPs
#' @other: a list containing: elements without names 

#include popIDs
pop(pdogGen)<- sapply(strsplit(indNames(pdogGen),"_"), function(pdogGen){pdogGen[1]})
levels(pop(pdogGen))<-c(paste0("CC",1:3),paste0("HE",1:4))

table(pop(pdogGen))
# CC1 CC2 CC3 HE1 HE2 HE3 HE4 
# 40  59  21  44  24  15  34 

#Let's start with a brief look at our data using the sambaR package
genlight2sambar("pdogGen",do_confirm=TRUE)

#filter data
filterdata(indmiss=0.25,snpmiss=0.1,min_mac=2,dohefilter=TRUE,snpdepthfilter=TRUE, min_spacing=500, nchroms=NULL, silent=TRUE,maxprop_hefilter = 0.06)

#calculate kinship
calckinship()

#genetic diversity
calcdiversity() 

###### Now let's identify loci under selection ###############

# https://https://academic.oup.com/mbe/article/37/7/2153/5826356
#source("https://https://github.com/bcm-uga/pcadapt/blob/master/README.md")

## Load prairie dog dataset ##
#this vcf has been filtered using vcftools
pdog.vcf<-"./raw.dplm.recode.vcf"
p.dog <- read.pcadapt("./raw.dplm.recode.vcf", type ="vcf")
pdog2=read.vcfR(pdog.vcf)
pdog2=vcfR2genlight(pdog2)
pop(pdog2)<- sapply(strsplit(indNames(pdog2),"_"), function(pdog2){pdog2[1]})
levels(pop(pdog2))<-c(paste0("CC",1:3),paste0("HE",1:4))
#No variant got discarded.
#Summary:

#- input file:				./raw.dplm.recode.vcf
#- output file:				/var/folders/gn/lk3nlh_n0qz4ccn0qn01bc8r0000gn/T//RtmpE9YjUG/file153b91b45280.pcadapt

#- #number of individuals detected:	233
#- number of loci detected:		2954

#2955 lines detected.
#233 columns detected.

#Choosing the number K of Principal Components
pdog.pca<-pcadapt(input = p.dog, K=20)

#Scree plot - up to K=20
plot(pdog.pca, option = "screeplot", plt.pkg = "ggplot")
#Focus in on K=10
plot(pdog.pca, option="screeplot", K=10, plt.pkg = "ggplot")

#PCA
#include popIDs


PC12=plot(pdog.pca, option = "scores", i = 1, j = 2, plt.pkg = "ggplot", pop=pdog2$pop)
PC23=plot(pdog.pca, option = "scores", i = 2, j = 3, plt.pkg = "ggplot", pop=pdog2$pop)
PC34=plot(pdog.pca, option = "scores", i =3, j = 4, plt.pkg = "ggplot", pop=pdog2$pop)
PC45=plot(pdog.pca, option = "scores", i = 4, j = 5, plt.pkg = "ggplot", pop=pdog2$pop)
PC56=plot(pdog.pca, option = "scores", i = 5, j = 6, plt.pkg = "ggplot", pop=pdog2$pop)

#calculate variance on each axis
EV <- (pdog.pca$singular.values^2)
EV2 <-EV*100
EV2

#Using ggplot to plot your PCA
PC12+theme_classic()+theme(plot.title = element_blank()) +
  xlab("PC1 (47.8%)")+
  ylab("PC2 (14.2%)")+
  scale_x_continuous(limits = c(-0.2, 0.2), 
                     breaks = c(-0.2,-0.1,0,0.1,0.2))+
  scale_y_continuous(limits = c(-0.2, 0.2), breaks = c(-0.2,-0.1,0,0.1,0.2))
ggsave(".PC12.pdf", dpi=600)

#Set your K value
pdog_k5 <- pcadapt(p.dog, K = 5)
summary(pdog_k5)

#Length Class  Mode   
#scores           1165  -none- numeric
#singular.values     5  -none- numeric
#loadings        14770  -none- numeric
#zscores         14770  -none- numeric
#af               2954  -none- numeric
#maf              2954  -none- numeric
#chi2.stat        2954  -none- numeric
#stat             2954  -none- numeric
#gif                 1  -none- numeric
#pvalues          2954  -none- numeric
#pass             2940  -none- numeric

#Basic Manhattan Plot
plot(pdog_k5, option="manhattan")

#Q-Q Plots
plot(pdog_k5, option="qqplot")

#Histograms of the test statistic and of the p-values
hist(pdog_k5$pvalues, xlab = "p-value",main = NULL, breaks = 50, col = "#2471A3")

plot(pdog_k5, option = "stat.distribution")

## q-values

#Choosing a cutoff for outlier detection
#To provide a list of outliers and choose a cutoff for outlier detection, there are several methods 
#that are listed below from the less conservative one to the more conservative one.

#For a given α (real valued number between 0 and 1), SNPs with q-values less than α will be considered as 
#outliers with an expected false discovery rate bounded by α. The false discovery rate is defined as the percentage 
#of false discoveries among the list of candidate SNPs. Here is an example of how to provide a list of candidate SNPs 
#for the geno3pops data, for an expected false discovery rate lower than 10%:
loci=pdog2$loc.names
qval <- qvalue(pdog_k5$pvalues)$qvalues
alpha10 <- 0.1
outliers_qvalue10 <- which(qval < alpha10)
length(outliers_qvalue10)
alpha5 <- 0.05
outliers_qvalue5.ints <- which(qval < alpha5)
length(outliers_qvalue5.ints)
outliers_qvalue5 <- as.character(loci)[outliers_qvalue5.ints]
outliers_qvalue5
#write.csv(outliers_qvalue5, "./827_adaptive_snps_q.csv")

#Benjamini-Hochberg Procedure
padj_BH <- p.adjust(pdog_k5$pvalues,method="BH")
alpha10 <- 0.1
outliers_BH10 <- which(padj_BH < alpha10)
length(outliers_BH10)
#outliers_BH10

alpha5 <- 0.05
outliers_BH5.ints <- which(padj_BH < alpha5)
length(outliers_BH5.ints)
outliers_BH5 <- as.character(loci)[outliers_BH5.ints]
outliers_BH5


alpha01 <- 0.01
outliers_BH1.ints <- which(padj_BH < alpha01)
length(outliers_BH1.ints)
outliers_BH1 <- as.character(loci)[outliers_BH1.ints]
outliers_BH1

#Bonferroni correction
padj_BC <- p.adjust(pdog_k5$pvalues,method="bonferroni")
alpha10 <- 0.1
outliers_BH10 <- which(padj_BC < alpha10)
length(outliers_BH10)
#outliers_BH10

alpha5 <- 0.05
outliers_BC5.ints <- which(padj_BC < alpha5)
length(outliers_BC5.ints)
outliers_BC5 <- as.character(loci)[outliers_BC5.ints]
outliers_BC5

#If you want the pvalues for each SNP
snps_pvalues <- cbind(loci, pdog_k5$pvalues)
snps_pvalues_no_na <- na.omit(snps_pvalues)
write.table(snps_pvalues, "All_Pvalues.txt", sep="\t", quote=FALSE)

#Population Structure

#Convert vcf file to lfmm and geno file type
pdog_lfmm<- vcf2lfmm('./raw.dplm.recode.vcf', output.file = "./pdog.lfmm", force = TRUE)
#- number of detected individuals:	233
#- number of detected loci:		2955
pdog_genos<-vcf2geno("./raw.dplm.recode.vcf", output.file = "./pdog.geno", force=TRUE)
#For SNP info, please check ./raw.dplm.recode.vcfsnp.

#0 line(s) were removed because these are not SNPs.
#Please, check ./raw.dplm.recode.removed file, for more informations.

# number of detected individuals:	233
# number of detected loci:		2955

#run of pca # Available options, K (the number of PCs), 
# center and scale. # Create files: genotypes.eigenvalues- eigenvalues, 
# # # genotypes.eigenvectors- eigenvectors, genotypes.sdev- standard deviations, genotypes.projections- projections, 
# Create a pcaProject object: 

# main options # K = number of ancestral populations 
# entropy = TRUE: computes the cross-entropy criterion, 
# CPU = 4 the number of CPUs. project = NULL project = 
pdog_project<-snmf("pdog.geno", K = 1:10, entropy = TRUE, repetitions = 10, project = "new")

# plot cross-entropy criterion for all runs in the snmf project 
pdog_prj_plot=plot(pdog_project, col="midnightblue",cex=1.5,pch=19)
#Changing format of data to plot CE with error bars

library(tidyverse)
PD <- summary(pdog_project)$crossEntropy %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("temp") %>%
  mutate(K = as.numeric(str_extract(temp, "(\\d)+"))) %>%
  select(-temp)

library(ggpubr)
PD5=ggplot(PD, aes(x = K, y = mean)) +
  geom_line(color = "black", size = 0.25 ) +
  geom_segment(aes(x = K, y = min, xend = K, yend = max)) +
  geom_point(shape = 21, size = 3, color = "black", fill = "slateblue") +
  scale_x_continuous(breaks = seq(0, 20, by = 1)) +
  labs(x = "Number of Ancestral Populations", y = "Cross-entropy Criterion") +
  theme_classic()+
  ggtitle("Praire Dog Cross Entropy")
PD5
  
# select the best run for K=3 through K = 5
#K=3
best3 = which.min(cross.entropy(pdog_project, K = 3)) 
my.colors <- c("darkslateblue", "darkturquoise", "burlywood", "darkmagenta","gold") 
barchart(pdog_project, K = 3, run = best2, 
         border = NA, space = 0, col = 
           my.colors, xlab = "Individuals", 
         ylab = "Ancestry proportions", main = "Ancestry matrix")-> bp 
axis(1, at = 1:length(bp$order), labels = bp$order, las=1, cex.axis = .4)

#K=4
best4 = which.min(cross.entropy(pdog_project, K = 4)) 
barchart(pdog_project, K = 4, run = best2, 
         border = NA, space = 0, col = 
           my.colors, xlab = "Individuals", 
         ylab = "Ancestry proportions", main = "Ancestry matrix")-> bp 
axis(1, at = 1:length(bp$order), labels = bp$order, las=1, cex.axis = .4)

#K=5
best5 = which.min(cross.entropy(pdog_project, K = 5)) 
barchart(pdog_project, K = 5, run = best, 
         border = NA, space = 0, col = 
           my.colors, xlab = "Individuals", 
         ylab = "Ancestry proportions", main = "Ancestry matrix")-> bp 


##### We can also look at outlier loci using FST and a manhattan plot

## Let's start off by reading in the VCF file
file<-read.vcfR("raw.dplm.recode.vcf",verbose=FALSE)
indiv_ids<-as.data.frame(colnames(file@gt))
indiv_ids<-tail(indiv_ids,-1)
colnames(indiv_ids)[1]<-"Individual_IDs"
test<-str_split_fixed(indiv_ids$Individual_IDs,'_', 2)
pop_names<-data.frame(pops=test[,1],Numeric_IDs=test[,2])
pop_names$id<-paste0(pop_names$pops,'_',pop_names$Numeric_IDs)
pop_names<-pop_names[c("pops","id")]
# How many populations do we have?
length(unique(pop_names$pops))
pop <- pop_names[order(pop_names$id),]

## Now that we have a population file, we also need to do a few things
# Let's start by making a genind file
my_genind <- vcfR2genind(file)
my_genind@pop<-as.factor(pop$pops)

## Let's take a quick look at this object to make sure it looks as we expect
my_genind

##/// GENIND OBJECT /////////
#  // 233 individuals; 2,955 loci; 5,910 alleles; size: 6.8 Mb

# // Basic content
# @tab:  233 x 5910 matrix of allele counts
# @loc.n.all: number of alleles per locus (range: 2-2)
# @loc.fac: locus factor for the 5910 columns of @tab
# @all.names: list of allele names for each locus
# @ploidy: ploidy of each individual  (range: 2-2)
# @type:  codom
# @call: adegenet::df2genind(X = t(x), sep = sep)

# // Optional content
# @pop: population of each individual (group size range: 13-58)

## Now, let's calculate pariwise FSTs according to Weir and Cockerham
# This will calculate pariwise FST for each population pair that we encoded into our data
FST<-pairwise.WCfst(dat=my_genind,diploid=TRUE)
stats<-basic.stats(data=my_genind)
locus_fst<-data.frame(ID=rownames(stats$perloc),FST=stats$perloc$Fst,P=stats$perloc$Fstp)
hist(locus_fst$P)

## So one thing to keep in mind, this is looking at per locus FST across all pairwise population comparisons
# Let's just look at two population pairs: 14CCUT3 and 14HEUT4
# Let's do a
subset<-my_genind[pop=c("14HEUT3","14HEUT4")]
# We'll quickly check to make sure it has roughly the same number of loci/individuals
subset
FST<-pairwise.WCfst(dat=subset,diploid=TRUE)
FST
sub_stats<-basic.stats(data=subset)
sub_locus_fst<-data.frame(ID=rownames(sub_stats$perloc),FST=sub_stats$perloc$Fst,P=sub_stats$perloc$Fstp)
hist(sub_locus_fst$P)

## We're now going to make a reference file for the chromosome/contig IDs and positions
# These genomic data are in contigs, not in chromosomes, so we're going to generate some random bins
contigs<-data.frame(file@fix)
contigs$CONTIG<-sub('_pilon','',contigs$CHROM)
snp_position<-contigs[c("POS","ID","REF","ALT","CONTIG")]
snp_position$Physical_SNP_Position<-paste(snp_position$CONTIG,snp_position$POS,sep="_")

## Now we'll merge our FST dataframe with our SNP loci dataframe
sub_locus_fst$ID<-sub('X','',sub_locus_fst$ID)
df<-merge(sub_locus_fst,snp_position,by="ID")
final_df<-df[!df$FST=="NaN",]
final_df$FST[final_df$FST < 0] <- "0"
final_df$FST<-as.numeric(final_df$FST)
final_df$POS<-as.numeric(final_df$POS)
df<-final_df[order(final_df$POS),]

## Time for some fun figures- let's generate a manhattan plot for each chromosome for these measures
# Very basic plot- but shows you what you can do
p1<-ggman(df,snp="Physical_SNP_Position",chrom="CONTIG",pvalue="FST",bp="POS",ymin=0,ymax=.8,
          pointSize=1,logTransform=FALSE)
p1<-p1+scale_color_manual(values=c("midnightblue","cyan3"))
p1<-p1+theme(axis.text.y=element_text(size=10),axis.title.y=element_text(size=15,face="bold"))
p1<-p1+theme(axis.text.x=element_text(size=10),axis.title.x=element_text(size=15,face="bold"))
p1<-p1+xlab(label="Contigs")
p1<-p1+ylab(label="Paiwsie FST")
p1<-p1+geom_hline(yintercept=0.12,linetype="dashed",color="red")
p1<-p1+labs(title="Pairwise FST 14HEUT3 and 14HEUT4")
p1<-p1+theme(strip.text=element_text(face="bold",size=20),strip.background=element_rect(fill="white",color="white"))
p1<-p1+theme(legend.text=element_text(size=rel(0.8)))
p1<-p1+theme(legend.title=element_text(face="bold",size=rel(1.1)))
p1
