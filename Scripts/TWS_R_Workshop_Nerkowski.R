#### TWS 2022 Genomic Analysis Workshop - pcadpt #####
## Load Packages ##
install.packages("vcfR")
#install.packages("sambaR")

library(vcfR)
# https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13339
source("https://github.com/mennodejong1986/SambaR/raw/master/SAMBAR_v1.08.txt")
getpackages(mylib=NULL)

## Load prairie dog dataset ##
setwd("../TWS R Workshop/")

p.dog_vcf <- read.vcfR("raw.dplm.recode.vcf")

# ***** Object of Class vcfR *****
# 60 samples
# 7911 CHROMs
# 773,157 variants
# Object size: 1961.2 Mb
# 0 percent missing data
# *****        *****         *****

pdogGen <- vcfR2genlight(pdog_vcf)
#/// GENLIGHT OBJECT /////////

#'   // 233 genotypes,  2,955 binary SNPs, size: 995.8 Kb
#' 28919 (4.2 %) missing data
#' 
#' // Basic content
#' @gen: list of 233 SNPbin
#' 
#' // Optional content
#' @ind.names:  233 individual labels
#' @loc.names:  2955 locus labels
#' @chromosome: factor storing chromosomes of the SNPs
#' @position: integer storing positions of the SNPs
#' @other: a list containing: elements without names 
#' 
#' 

# Create some code here to use indNames and extract the 3-7 characters for pop name
# Maybe include geographic coordinates

genlight2sambar("pdogGen",do_confirm=TRUE)

#filter data
library(gplots) #needed for heatmap.2 in filterdata
filterdata(indmiss=0.25,snpmiss=0.1,min_mac=2,dohefilter=TRUE,snpdepthfilter=TRUE, min_spacing=500, nchroms=NULL, silent=TRUE)
library(ape)
library(poppr)
install.packages("BiocManager")
library(BiocManager)
install.packages("StAMPP")
library(StAMPP)
BiocManager::install("LEA")

library(LEA)
findstructure(Kmax=6,add_legend=TRUE,legend_pos='bottomright',legend_cex=3,symbol_size=3)

##Population Structure and Outlier Loci
## Load Packages ##
install.packages("pcadapt")

library(pcadapt)
# https://https://academic.oup.com/mbe/article/37/7/2153/5826356
#source("https://https://github.com/bcm-uga/pcadapt/blob/master/README.md")

## Load prairie dog dataset ##
setwd("~/Desktop/MEWG/TWS R Workshop/")
pdog.vcf<-"./raw.dplm.recode.vcf"
p.dog <- read.pcadapt("./raw.dplm.recode.vcf", type ="vcf")
#No variant got discarded.
#Summary:
  
#- input file:				./raw.dplm.recode.vcf
#- output file:				/var/folders/gn/lk3nlh_n0qz4ccn0qn01bc8r0000gn/T//RtmpE9YjUG/file153b91b45280.pcadapt

#- #number of individuals detected:	233
#- number of loci detected:		2954

#2954 lines detected.
#233 columns detected.

#Choosing the number K of Principal Components
pdog.pca<-pcadapt(input = p.dog, K=20)

#Scree plot - up to K=20
plot(pdog.pca, option = "screeplot")
#Focus in on K=10
plot(pdog.pca, option="screeplot", K=10)

#PCA
plot(pdog.pca, option = "scores", i = 1, j = 2)
plot(pdog.pca, option = "scores", i = 2, j = 3)
plot(pdog.pca, option = "scores", i =3, j = 4)
plot(pdog.pca, option = "scores", i = 4, j = 5)
plot(pdog.pca, option = "scores", i = 5, j = 6)

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
#The R package qvalue, transforms p-values into q-values. To install and load the package, 
#type the following command lines:
install.packages("devtools")
library("devtools")
install_github("jdstorey/qvalue")

#Choosing a cutoff for outlier detection
#To provide a list of outliers and choose a cutoff for outlier detection, there are several methods 
#that are listed below from the less conservative one to the more conservative one.

#For a given α (real valued number between 0 and 1), SNPs with q-values less than α will be considered as 
#outliers with an expected false discovery rate bounded by α. The false discovery rate is defined as the percentage 
#of false discoveries among the list of candidate SNPs. Here is an example of how to provide a list of candidate SNPs 
#for the geno3pops data, for an expected false discovery rate lower than 10%:

library(qvalue)
qval <- qvalue(pdog_k5$pvalues)$qvalues
alpha10 <- 0.1
outliers_qvalue10 <- which(qval < alpha10)
length(outliers_qvalue10)
alpha5 <- 0.05
outliers_qvalue5 <- which(qval < alpha5)
length(outliers_qvalue5)


#Benjamini-Hochberg Procedure
padj_BH <- p.adjust(pdog_k5$pvalues,method="BH")
alpha10 <- 0.1
outliers_BH10 <- which(padj_BH < alpha10)
length(outliers_BH10)
#outliers_BH10

alpha5 <- 0.05
outliers_BH5 <- which(padj_BH < alpha5)
length(outliers_BH5)
#outliers_BH5

alpha01 <- 0.01
outliers_BH01 <- which(padj_BH < alpha01)
length(outliers_BH01)
#outliers_BH01

#Bonferroni correction
padj_BC <- p.adjust(pdog_k5$pvalues,method="bonferroni")
alpha10 <- 0.1
outliers_BC10 <- which(padj_BC < alpha10)
length(outliers_BC10)
#outliers_BC10

alpha10 <- 0.05
outliers_BC5 <- which(padj_BC < alpha5)
length(outliers_BC5)
#outliers_BC5

write.table(outliers_BH10, file = "./outlier_BH10.txt", sep = "\t", row.names = FALSE)
write.table(outliers_BH5, file = "./outlier_BH5.txt", sep = "\t", row.names = FALSE)
write.table(outliers_BC10, file = "./outlier_BC10.txt", sep = "\t", row.names = FALSE)
write.table(outliers_BC5, file = "./Aoutlier_BC5.txt", sep = "\t", row.names = FALSE)
write.table(padj_BC, file = "./padj_BC.txt",sep = "\t", row.names = FALSE)
write.table(qval, file = "./qvalues.txt",sep = "\t", row.names = FALSE)

#Population Structure
#Install packages if needed
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("LEA")

#Load library
library(LEA)

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
pdog_pc = pca(pdog_lfmm, scale = TRUE)
"******************************"
#[1] " Principal Component Analysis "
#[1] "******************************"
#summary of the options:
  
#  -n (number of individuals)          233
#-L (number of loci)                 2955
#-K (number of principal components) 233
#-x (genotype file)                  /Users/Gaea/Desktop/MEWG/TWS R Workshop/raw.dplm.recode.lfmm
#-a (eigenvalue file)                /Users/Gaea/Desktop/MEWG/TWS R Workshop/raw.dplm.recode.pca/raw.dplm.recode.eigenvalues
#-e (eigenvector file)               /Users/Gaea/Desktop/MEWG/TWS R Workshop/raw.dplm.recode.pca/raw.dplm.recode.eigenvectors
#-d (standard deviation file)        /Users/Gaea/Desktop/MEWG/TWS R Workshop/raw.dplm.recode.pca/raw.dplm.recode.sdev
#-p (projection file)                /Users/Gaea/Desktop/MEWG/TWS R Workshop/raw.dplm.recode.pca/raw.dplm.recode.projections
#-s data centered and scaled 

# Perfom Tracy-Widom tests on all eigenvalues. # create file: tuto.tracyWidom- tracy-widom test information. 
tw = tracy.widom(pdog_pc)

# display p-values for the Tracy-Widom tests (first 5 pcs). 
tw$pvalues[1:20]

# plot the percentage of variance explained by each component 
plot(tw$percentage)

# main options # K = number of ancestral populations 
# entropy = TRUE: computes the cross-entropy criterion, 
# CPU = 4 the number of CPUs. project = NULL project = 
pdog_project<-snmf("pdog.geno", K = 1:10, entropy = TRUE, repetitions = 10, project = "new")

# plot cross-entropy criterion for all runs in the snmf project 
plot(pdog_project, col = "blue", pch = 19, cex = 1.2)

# select the best run for K=3 through K = 5
#K=3
best2 = which.min(cross.entropy(pdog_project, K = 3)) 
my.colors <- c("darkslateblue", "darkturquoise", "burlywood", "darkmagenta","gold") 
barchart(pdog_project, K = 3, run = best2, 
          border = NA, space = 0, col = 
           my.colors, xlab = "Individuals", 
         ylab = "Ancestry proportions", main = "Ancestry matrix")-> bp 
axis(1, at = 1:length(bp$order), labels = bp$order, las=1, cex.axis = .4)

#K=4
best2 = which.min(cross.entropy(pdog_project, K = 4)) 
barchart(pdog_project, K = 4, run = best2, 
         border = NA, space = 0, col = 
           my.colors, xlab = "Individuals", 
         ylab = "Ancestry proportions", main = "Ancestry matrix")-> bp 
axis(1, at = 1:length(bp$order), labels = bp$order, las=1, cex.axis = .4)

#K=5
best = which.min(cross.entropy(pdog_project, K = 5)) 
barchart(pdog_project, K = 5, run = best, 
         border = NA, space = 0, col = 
           my.colors, xlab = "Individuals", 
         ylab = "Ancestry proportions", main = "Ancestry matrix")-> bp 

