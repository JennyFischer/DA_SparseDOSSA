library(microbiomeMarker)
library(microbiome)
library(compositions)
library(phyloseq)
library(ALDEx2)
library(tibble)
library(dplyr)

#FIRST GET ALL INPUT DATA
#whole otu file from output of sparseDOSSA
otu_t <- read.delim("./first-dataset/otu_table.tsv")

#OTU file for the null (i.e., not spiked-in) absolute abundances
#first seperate python script and manual deletion of first two columns
lognormal <- read.delim("./first-dataset/lognorm.tsv", header=FALSE, row.names=NULL)
rownames(lognormal) <- paste0("OTU", 1:nrow(lognormal))
colnames(lognormal) <- paste0("Sample", 1:ncol(lognormal))
OTU_lognormal = otu_table(lognormal, taxa_are_rows = TRUE)

#OTU file for spiked-in absolute abundancess
#first seperate python script and manual deletion of first two columns
logspike <- read.delim("./first-dataset/spike.tsv", header=FALSE, row.names=NULL)
rownames(logspike) <- paste0("OTU", 1:nrow(logspike))
colnames(logspike) <- paste0("Sample", 1:ncol(logspike))
OTU_logspike = otu_table(logspike, taxa_are_rows = TRUE)


n = ncol(OTU_logspike)
groups = sample(c(0,1), replace=TRUE, size=n)
#dont need normalisation beacuse aldex does by itself

#SECOND run aldex for both datasets
#ABS
x <- aldex.clr(OTU_lognormal, groups)
x.tt <- aldex.ttest(x)
abs =select(x.tt, we.ep, wi.ep )
#REL
y <- aldex.clr(OTU_logspike, groups)
y.tt <- aldex.ttest(y)
rel = select(y.tt, we.ep, wi.ep )

#Get overlap of both for wilcoxon 
sign_abs_wilc <- filter(abs[2], abs[2] < 0.05)
sign_rel_wilc <- filter(rel[2], rel[2] < 0.05)
intersect(sign_rel_wilc, sign_rel_wilc) 


#Get overlap of both for welch ttest 
sign_abs_welch <- filter(abs[1], abs[1] <0.05)
sign_rel_welch <- filter(rel[1], rel[1] <0.05)
intersect(sign_abs_welch, sign_rel_welch) 

