##### TWS 2024 Genomic Analysis Workshop #####

##Set working directory
##Use the path to where you downloaded the data ("batch_1.vcf")
setwd("C:/Users/tbost/OneDrive/Documents/TAMUK/SPEC Lab Materials/GitHub/TWS_GenomicWS_24/Data")

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

#### Part 1: SambaR ####
### --- in samabaR; creates many, many figures; those figures not uploaded --- ###
### those figures can be tailored by using arguments found in the handbook ###
### this section of the code goes over basic samba codes that gives you overview of your SNP data
## not necessarily analyses, but good for prelim

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
#Prelim Filters:1 SNP per contig (keep loci independent), 
#loci represented in at least 80% of indiv, 
#removed loci w/ het over .7 (filter to account for assembly errors)\
#Minor allele frequency (MAF) removed if less than 5%
#done in Stacks pipeline v1.48
#vcf file is the output from the pipeline -- not the raw
vcf

# ***** Object of Class vcfR *****
#   237 samples; check to make sure this is the expected number
# 1826 CHROMs
# 3,549 variants
# Object size: 8.2 Mb
# 0 percent missing data
# *****        *****         *****

vcf@meta

vcf@fix

vcf@gt

pdogGen <- vcfR2genlight(vcf) #making vcf file into genlight file
pdogGen
#' /// GENLIGHT OBJECT ///////// : genlight is a type of input file, some R packages want other types such as Genid and genpop
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

indNames(pdogGen) #shows population/name ID's

#include popIDs
pop(pdogGen)<- sapply(strsplit(indNames(pdogGen),"_"), function(pdogGen){pdogGen[1]})
#strsplit seperates ID's by _ and then function keeps the first half to create the population ID
levels(pop(pdogGen))<-c(paste0("CC",1:3),paste0("HE",1:4))
#levels: categorical, code above is renaming the pop ID's to be more workable
#make sure that when renaming they are in the right order
levels(pop(pdogGen))

table(pop(pdogGen)) #count of how many individ in each pop
# CC1 CC2 CC3 HE1 HE2 HE3 HE4 : the result of renaming the pop Id's 
# 40  59  21  44  24  15  34 

#Let's start with a brief look at our data using the sambaR package
#samabaR -- look at assumptions made by package, this package good for surface level analyses; a collection of Rpackages for genomics
#samabaR -- a wrapper for multiple packages but can use each package individually; automatically generates figures for quick peeks
#default versus alternate arguments and filtering options
#make multiple datasets with different filtering to answer different questions
genlight2sambar("pdogGen",do_confirm=TRUE)

#filter data : can only work on one project at a time in sambaR
filterdata(indmiss=0.25,snpmiss=0.1,min_mac=2,dohefilter=TRUE,snpdepthfilter=TRUE, min_spacing=500, nchroms=NULL, silent=TRUE,maxprop_hefilter = 0.06)
#filter conditions: indmiss (level of missing data per indivd), min_mac (min minor allele count), dohefilter:hetero filter, min_spacing(minimum space between snps)
#in sambaR_outputs file folder: automatically generating figures 

#calculate kinship
calckinship()
#outputs at sambaR Kinship
#inbreeding coefficents: outliers colored in, normal open circle; each dot an individual
#heterozygosities lower for SNPs than for microsatts -- on average half, so healthy pop average .1-.3 range; will always have lower levels

#genetic diversity
calcdiversity() 

#calculate genetic distance
calcdistance()

#pop structure
findstructure(Kmax=6)
#pick high k initially
#cross-entropy criterion: helps select correct k -- when the line dips first (ex 5)

##### Part 2: Identify loci under selection #####

## go to presenter script for part 2 ##

# https://https://academic.oup.com/mbe/article/37/7/2153/5826356
#source("https://https://github.com/bcm-uga/pcadapt/blob/master/README.md")

## Load prairie dog dataset ##
#this vcf has been filtered using vcftools
pdog.vcf<-"./raw.dplm.recode.vcf"
p.dog <- read.pcadapt("./raw.dplm.recode.vcf", type ="vcf")



#Choosing the number K of Principal Components


#Scree plot - up to K=20

#Focus in on K=10


#PCA


#Set your K value


#Basic Manhattan Plot


#Q-Q Plots


#Histograms of the test statistic and of the p-values


## q-values

#Choosing a cutoff for outlier detection
#To provide a list of outliers and choose a cutoff for outlier detection, there are several methods 
#that are listed below from the less conservative one to the more conservative one.

#For a given α (real valued number between 0 and 1), SNPs with q-values less than α will be considered as 
#outliers with an expected false discovery rate bounded by α. The false discovery rate is defined as the percentage 
#of false discoveries among the list of candidate SNPs. Here is an example of how to provide a list of candidate SNPs 
#for the geno3pops data, for an expected false discovery rate lower than 10%:




#Benjamini-Hochberg Procedure


#Bonferroni correction


#Write out results to txt file
write.table(outliers_BH10, file = "./outlier_BH10.txt", sep = "\t", row.names = FALSE)
write.table(outliers_BH5, file = "./outlier_BH5.txt", sep = "\t", row.names = FALSE)
write.table(outliers_BC10, file = "./outlier_BC10.txt", sep = "\t", row.names = FALSE)
write.table(outliers_BC5, file = "./Aoutlier_BC5.txt", sep = "\t", row.names = FALSE)
write.table(padj_BC, file = "./padj_BC.txt",sep = "\t", row.names = FALSE)
write.table(qval, file = "./qvalues.txt",sep = "\t", row.names = FALSE)

#Population Structure

#Convert vcf file to lfmm and geno file type

#For SNP info, please check ./raw.dplm.recode.vcfsnp.





#run of pca # Available options, K (the number of PCs), 
# center and scale. # Create files: genotypes.eigenvalues- eigenvalues, 
# # # genotypes.eigenvectors- eigenvectors, genotypes.sdev- standard deviations, genotypes.projections- projections, 
# Create a pcaProject object: 

# Perfom Tracy-Widom tests on all eigenvalues. # create file: tuto.tracyWidom- tracy-widom test information. 


# display p-values for the Tracy-Widom tests (first 5 pcs). 


# plot the percentage of variance explained by each component 


# main options # K = number of ancestral populations 
# entropy = TRUE: computes the cross-entropy criterion, 
# CPU = 4 the number of CPUs. project = NULL project = 

# plot cross-entropy criterion for all runs in the snmf project 

# select the best run for K=3 through K = 5
#K=3

#K=4

#K=5



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
p1<-p1+scale_color_manual(values=c("grey19","red"))
p1<-p1+theme(axis.text.y=element_text(size=10),axis.title.y=element_text(size=15,face="bold"))
p1<-p1+theme(axis.text.x=element_text(size=10),axis.title.x=element_text(size=15,face="bold"))
p1<-p1+xlab(label="Contigs")
p1<-p1+ylab(label="Paiwsie FST")
p1<-p1+geom_hline(yintercept=0.12,linetype="dashed",color="darkblue")
p1<-p1+labs(title="Pairwise FST 14HEUT3 and 14HEUT4")
p1<-p1+theme(strip.text=element_text(face="bold",size=20),strip.background=element_rect(fill="white",color="white"))
p1<-p1+theme(legend.text=element_text(size=rel(0.8)))
p1<-p1+theme(legend.title=element_text(face="bold",size=rel(1.1)))
p1


