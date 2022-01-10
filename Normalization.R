#need to be installed 
#BiocManager::install("microbiome")

library(microbiomeMarker)
library(microbiome)
library(compositions)
library(phyloseq)
library(ALDEx2)

#whole otu file from output of sparseDOSSA
otu_t <- read.delim("~/Desktop/otu_table.tsv")

#OTU file for the null (i.e., not spiked-in) absolute abundances
#first seperate python script an manual deletion of first two columns
lognormal <- read.delim("~/Desktop/lognorm.tsv", header=FALSE, row.names=NULL)
rownames(lognormal) <- paste0("OTU", 1:nrow(lognormal))
colnames(lognormal) <- paste0("Sample", 1:ncol(lognormal))
OTU_lognorm = otu_table(lognormal, taxa_are_rows = TRUE)

#OTU file for spiked-in absolute abundancess
#first seperate python script an manual deletion of first two columns
logspike <- read.delim("~/Desktop/spike.tsv", header=FALSE, row.names=NULL)
rownames(logspike) <- paste0("OTU", 1:nrow(logspike))
colnames(logspike) <- paste0("Sample", 1:ncol(logspike))
OTU_logspike = otu_table(logspike, taxa_are_rows = TRUE)

############ NORMALIZATION ####################################################
#for the null (i.e., not spiked-in) absolute abundances
#clr normalization
norm_clr_lognorm =norm_clr(OTU_lognorm)
norm_clr_lognorm 
#microbiome::transform(OTU, 'clr') # gibt negative Werte
write.csv(norm_clr_lognorm, file = "Normalisation_output/normalise_lognorm_clr.csv", row.names = TRUE)

#css normalization
norm_css_lognorm = normalize(lognormal, "CSS")
norm_css_lognorm
write.csv(norm_css_lognorm, file = "Normalisation_output/normalise_lognorm_css.csv", row.names = TRUE)

#tss normalization
norm_tss_lognorm = normalize(lognormal, "TSS")#norm_tss(OTU)
norm_tss_lognorm
write.csv(norm_tss_lognorm, file = "Normalisation_output/normalise_lognorm_tss.csv", row.names = TRUE)

#relative abundance normalization
norm_relab_lognorm = transform(lognormal, 'compositional')
norm_relab_lognorm
write.csv(norm_relab_lognorm, file = "Normalisation_output/normalise_lognorm_relab.csv", row.names = TRUE)


##############################################################################
#for OTU file for spiked-in absolute abundancess
#clr normalization
norm_clr_logspike =norm_clr(OTU_logspike)
norm_clr_logspike 
#microbiome::transform(OTU, 'clr') # gibt negative Werte
write.csv(norm_clr_logspike, file = "Normalisation_output/normalise_logspike_clr.csv", row.names = TRUE)

#css normalization
norm_css_logspike = normalize(lognormal, "CSS")
norm_css_logspike
write.csv(norm_css_logspike, file = "Normalisation_output/normalise_logspike_css.csv", row.names = TRUE)

#tss normalization
norm_tss_logspike = normalize(lognormal, "TSS")#norm_tss(OTU)
norm_tss_logspike
write.csv(norm_tss_logspike, file = "Normalisation_output/normalise_logspike_tss.csv", row.names = TRUE)

#relative abundance normalization
norm_relab_logspike = transform(lognormal, 'compositional')
norm_relab_logspike
write.csv(norm_relab_logspike, file = "Normalisation_output/normalise_logspike_relab.csv", row.names = TRUE)

