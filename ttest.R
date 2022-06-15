library(microbiomeMarker)
library(microbiome)
library(compositions)
library(phyloseq)
library(ALDEx2)
library(ggpubr)
library(rstatix)
library(dplyr)

#whole otu file from output of sparseDOSSA
otu_t <- read.delim("./first_dataset/otu_table.tsv")

#OTU file for the null (i.e., not spiked-in) absolute abundances
#first seperate python script an manual deletion of first two columns
lognormal <- read.delim("./first_dataset/lognorm.tsv", header=FALSE, row.names=NULL)
rownames(lognormal) <- paste0("OTU", 1:nrow(lognormal))
colnames(lognormal) <- paste0("Sample", 1:ncol(lognormal))
OTU_lognorm = otu_table(lognormal, taxa_are_rows = TRUE)
OTU_lognorm
write.csv(OTU_lognorm, file = "./first_dataset/OTU_lognorm.csv", row.names = TRUE)

#OTU file for spiked-in absolute abundancess
#first seperate python script an manual deletion of first two columns
logspike <- read.delim("./first_dataset/spike.tsv", header=FALSE, row.names=NULL)
rownames(logspike) <- paste0("OTU", 1:nrow(logspike))
colnames(logspike) <- paste0("Sample", 1:ncol(logspike))
OTU_logspike = otu_table(logspike, taxa_are_rows = TRUE)
write.csv(OTU_logspike, file = "./first_dataset/OTU_logspike.csv", row.names = TRUE)
#write.table(OTU_logspike, file = "/Users/jf/Desktop/OTU_logspike.tsv", row.names=TRUE, sep="\t")

############ NORMALIZATION ####################################################
###############################################################################
#RELATIVE ABUNDANCE
###############################################################################
#clr normalization
norm_clr_logspike =norm_clr(OTU_logspike)
norm_clr_logspike 
#microbiome::transform(OTU, 'clr') # gibt negative Werte
write.csv(norm_clr_logspike, file = "./first_dataset/Normalisation_output/normalise_logspike_clr.csv", row.names = TRUE)

#css normalization
norm_css_logspike = normalize(logspike, "CSS")
norm_css_logspike
write.csv(norm_css_logspike, file = "./first_dataset/Normalisation_output/normalise_logspike_css.csv", row.names = TRUE)

#tss normalization
norm_tss_logspike = normalize(logspike, "TSS")#norm_tss(OTU)
norm_tss_logspike
write.csv(norm_tss_logspike, file = "./first_dataset/Normalisation_output/normalise_logspike_tss.csv", row.names = TRUE)

#relative abundance normalization
norm_relab_logspike = transform(logspike, 'compositional')
norm_relab_logspike
write.csv(norm_relab_logspike, file = "./first_dataset/Normalisation_output/normalise_logspike_relab.csv", row.names = TRUE)





###############################################################################
############  ABSOLUTE #################
#build TAX file 
OTU_lognormal = otu_table(lognormal, taxa_are_rows = TRUE)
taxmat_lognormal = matrix(sample(letters, 100, replace = TRUE), nrow = nrow(OTU_lognormal), ncol = 7)
rownames(taxmat_lognormal) <- rownames(OTU_lognormal)
colnames(taxmat_lognormal) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
TAX_lognormal = tax_table(taxmat_lognormal)
physeq_lognormal = phyloseq(OTU_lognormal, TAX_lognormal)
#write.table(OTU_lognormal_norm_relab, file = "/Users/jf/Desktop/OTU_lognormal.tsv", row.names=TRUE, sep="\t")

#Sample file with Metadata 
conds_lognormal = c(1,	0,	1,	1,	1,	1,	0,	1,	1,	1,	0,	0,	1,	1,	0,	1,	0,	1,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	1,	1,	1,	1,	0,	0,	0,	1,	0,	1,	0,	0,	1,	1,	0,	1,	1,	1,	1,	0,	0,	0,	1,	1,	0,	0,	1,	0,	0,	0,	0,	1,	0,	0,	1,	0,	0,	0,	0,	0,	0,	1,	0,	0,	1,	1,	0,	1,	1,	0,	0,	1,	0,	0,	1,	0,	1,	1,	1,	0,	1,	1,	0,	1,	1,	1,	0,	0)
sample_data_lognormal = sample_data(data.frame(Group = sample(conds_lognormal, size=nsamples(physeq_lognormal), 
                                                                         replace=TRUE),row.names=sample_names(physeq_lognormal),stringsAsFactors=FALSE))
#write.csv(sample_data, file = "metadata.csv", row.names = TRUE)

#physeq class
physeq_lognormal = phyloseq(OTU_lognormal, TAX_lognormal, sample_data_lognormal)
physeq_lognormal



smp_lognormal = row.names=sample_names(physeq_lognormal)

class.df<- data.frame(conds_lognormal, smp_lognormal) 

data_grp0_lognormal = subset(class.df, conds_lognormal == 0) # only metadata 0
data_grp1_lognormal = subset(class.df, conds_lognormal == 1) # only metadata 1

test_ttest_lognormal = lognormal
colnames_ttest_lognormal = ifelse(conds_lognormal == 0,"cond_zero","cond_one")
colnames_ttest_lognormal
colnames(test_ttest_lognormal) <- paste0(colnames_ttest_lognormal, 1:ncol(test_ttest_lognormal))
test_ttest_lognormal
write.csv(test_ttest_lognormal, file = "./first_dataset/ttest_output/test_ttest_lognormal.csv", row.names = TRUE)

#!!!!!!!!!!!!   einfugen von "gene_name"
test_lognormal <- read.csv("./first_dataset/ttest_output/test_ttest_lognormal.csv")


dlong_lognormal <- as.tbl(test_lognormal) %>% 
  gather(key, value, -gene_name) %>% 
  mutate(group=ifelse(grepl("cond_zero",key), "control", "case"))

tt_pvalues_all_lognormal  = dlong_lognormal %>% 
  group_by(gene_name) %>% 
  t_test(value~group)%>%
  adjust_pvalue(method = "BH")%>%
  add_significance()
tt_pvalues_all_lognormal
write.csv(tt_pvalues_all_lognormal , file = "./first_dataset/ttest_output/tt_pvalues_lognormal_all.csv", row.names = FALSE)


tt_pvalues_lognormal  = dlong_lognormal %>% 
  group_by(gene_name) %>% 
  summarise(p=t.test(value~group)$p.value)
tt_pvalues_lognormal
write.csv(tt_pvalues_lognormal , file = "./first_dataset/ttest_output/tt_pvalues_lognormal.csv", row.names = FALSE)






############  RELATVE #################
############  CLR #################
#build TAX file 
OTU_logspike_norm_clr = otu_table(norm_clr_logspike, taxa_are_rows = TRUE)
taxmat_logspike_norm_clr = matrix(sample(letters, 100, replace = TRUE), nrow = nrow(OTU_logspike_norm_clr), ncol = 7)
rownames(taxmat_logspike_norm_clr) <- rownames(OTU_logspike_norm_clr)
colnames(taxmat_logspike_norm_clr) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
TAX_logspike_norm_clr = tax_table(taxmat_logspike_norm_clr)
physeq_logspike_norm_clr = phyloseq(OTU_logspike_norm_clr, TAX_logspike_norm_clr)

#Sample file with Metadata 
conds_logspike_norm_clr = c(1,	0,	1,	1,	1,	1,	0,	1,	1,	1,	0,	0,	1,	1,	0,	1,	0,	1,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	1,	1,	1,	1,	0,	0,	0,	1,	0,	1,	0,	0,	1,	1,	0,	1,	1,	1,	1,	0,	0,	0,	1,	1,	0,	0,	1,	0,	0,	0,	0,	1,	0,	0,	1,	0,	0,	0,	0,	0,	0,	1,	0,	0,	1,	1,	0,	1,	1,	0,	0,	1,	0,	0,	1,	0,	1,	1,	1,	0,	1,	1,	0,	1,	1,	1,	0,	0)
sample_data_logspike_norm_clr = sample_data(data.frame(Group = sample(conds_logspike_norm_clr, size=nsamples(physeq_logspike_norm_clr), 
                                                                      replace=TRUE),row.names=sample_names(physeq_logspike_norm_clr),stringsAsFactors=FALSE))
#write.csv(sample_data, file = "metadata.csv", row.names = TRUE)

#physeq class
physeq_logspike_norm_clr = phyloseq(OTU_logspike_norm_clr, TAX_logspike_norm_clr, sample_data_logspike_norm_clr)
physeq_logspike_norm_clr



smp_logspike_norm_clr = row.names=sample_names(physeq_logspike_norm_clr)
#o = c("OTU1","OTU2","OTU3","OTU4","OTU5","OTU6","OTU7","OTU8","OTU9","OTU10")
#o
class.df<- data.frame(conds_logspike_norm_clr, smp_logspike_norm_clr) 

data_grp0_logspike_norm_clr = subset(class.df, conds_logspike_norm_clr == 0) # only metadata 0
data_grp1_logspike_norm_clr = subset(class.df, conds_logspike_norm_clr == 1) # only metadata 1

test_ttest_logspike_norm_clr = norm_css_logspike
colnames_ttest_logspike_norm_clr = ifelse(conds_logspike_norm_clr == 0,"cond_zero","cond_one")
colnames_ttest_logspike_norm_clr
colnames(test_ttest_logspike_norm_clr) <- paste0(colnames_ttest_logspike_norm_clr, 1:ncol(test_ttest_logspike_norm_clr))
test_ttest_logspike_norm_clr
write.csv(test_ttest_logspike_norm_clr, file = "./first_dataset/ttest_output/test_ttest_logspike_norm_clr.csv", row.names = TRUE)

#!!!!!!!!!!!!   einfugen von "gene_name"
test_logspike_norm_clr <- read.csv("./first_dataset/ttest_output/test_ttest_logspike_norm_clr.csv")


dlong_logspike_norm_clr <- as.tbl(test_logspike_norm_clr) %>% 
  gather(key, value, -gene_name) %>% 
  mutate(group=ifelse(grepl("cond_zero",key), "control", "case"))

tt_pvalues_logspike_all_norm_clr = dlong_logspike_norm_clr %>% 
  group_by(gene_name) %>% 
  t_test(value~group)%>%
  adjust_pvalue(method = "BH")%>%
  add_significance()
tt_pvalues_logspike_all_norm_clr
write.csv(tt_pvalues_logspike_all_norm_clr, file = "./first_dataset//ttest_output/tt_pvalues_logspike_norm_clr_all.csv", row.names = FALSE)


tt_pvalues_logspike_norm_clr = dlong_logspike_norm_clr %>% 
  group_by(gene_name) %>% 
  summarise(p=t.test(value~group)$p.value)
tt_pvalues_logspike_norm_clr
write.csv(tt_pvalues_logspike_norm_clr, file = "./first_dataset/ttest_output/tt_pvalues_logspike_norm_clr.csv", row.names = FALSE)



############  CSS #################
#build TAX file 
OTU_logspike_norm_css = otu_table(norm_css_logspike, taxa_are_rows = TRUE)
taxmat_logspike_norm_css = matrix(sample(letters, 100, replace = TRUE), nrow = nrow(OTU_logspike_norm_css), ncol = 7)
rownames(taxmat_logspike_norm_css) <- rownames(OTU_logspike_norm_css)
colnames(taxmat_logspike_norm_css) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
TAX_logspike_norm_css = tax_table(taxmat_logspike_norm_css)
physeq_logspike_norm_css = phyloseq(OTU_logspike_norm_css, TAX_logspike_norm_css)

#Sample file with Metadata 
conds_logspike_norm_css = c(1,	0,	1,	1,	1,	1,	0,	1,	1,	1,	0,	0,	1,	1,	0,	1,	0,	1,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	1,	1,	1,	1,	0,	0,	0,	1,	0,	1,	0,	0,	1,	1,	0,	1,	1,	1,	1,	0,	0,	0,	1,	1,	0,	0,	1,	0,	0,	0,	0,	1,	0,	0,	1,	0,	0,	0,	0,	0,	0,	1,	0,	0,	1,	1,	0,	1,	1,	0,	0,	1,	0,	0,	1,	0,	1,	1,	1,	0,	1,	1,	0,	1,	1,	1,	0,	0)
sample_data_logspike_norm_css = sample_data(data.frame(Group = sample(conds_logspike_norm_css, size=nsamples(physeq_logspike_norm_css), 
                                                                      replace=TRUE),row.names=sample_names(physeq_logspike_norm_css),stringsAsFactors=FALSE))
#write.csv(sample_data, file = "metadata.csv", row.names = TRUE)

#physeq class
physeq_logspike_norm_css = phyloseq(OTU_logspike_norm_css, TAX_logspike_norm_css, sample_data_logspike_norm_css)
physeq_logspike_norm_css



smp_logspike_norm_css = row.names=sample_names(physeq_logspike_norm_css)

class.df<- data.frame(conds_logspike_norm_css, smp_logspike_norm_css) 

data_grp0_logspike_norm_css = subset(class.df, conds_logspike_norm_css == 0) # only metadata 0
data_grp1_logspike_norm_css = subset(class.df, conds_logspike_norm_css == 1) # only metadata 1

test_ttest_logspike_norm_css = norm_css_logspike
colnames_ttest_logspike_norm_css = ifelse(conds_logspike_norm_css == 0,"cond_zero","cond_one")
colnames_ttest_logspike_norm_css
colnames(test_ttest_logspike_norm_css) <- paste0(colnames_ttest_logspike_norm_css, 1:ncol(test_ttest_logspike_norm_css))
test_ttest_logspike_norm_css
write.csv(test_ttest_logspike_norm_css, file = "./first_dataset/ttest_output/test_ttest_logspike_norm_css.csv", row.names = TRUE)

#!!!!!!!!!!!!   einfugen von "gene_name"
test_logspike_norm_css <- read.csv("./first_dataset/ttest_output/test_ttest_logspike_norm_css.csv")


dlong_logspike_norm_css <- as.tbl(test_logspike_norm_css) %>% 
  gather(key, value, -gene_name) %>% 
  mutate(group=ifelse(grepl("cond_zero",key), "control", "case"))


tt_pvalues_logspike_all_norm_css = dlong_logspike_norm_css %>% 
  group_by(gene_name) %>% 
  t_test(value~group)%>%
  adjust_pvalue(method = "BH")%>%
  add_significance()
tt_pvalues_logspike_all_norm_css
write.csv(tt_pvalues_logspike_all_norm_css, file = "./first_dataset/ttest_output/tt_pvalues_logspike_norm_css_all.csv", row.names = FALSE)


tt_pvalues_logspike_norm_css = dlong_logspike_norm_css %>% 
  group_by(gene_name) %>% 
  summarise(p=t.test(value~group)$p.value)
tt_pvalues_logspike_norm_css
write.csv(tt_pvalues_logspike_norm_css, file = "./first_dataset/ttest_output/tt_pvalues_logspike_norm_css.csv", row.names = FALSE)




############  TSS  #################
#build TAX file 
OTU_logspike_norm_tss = otu_table(norm_tss_logspike, taxa_are_rows = TRUE)
taxmat_logspike_norm_tss = matrix(sample(letters, 100, replace = TRUE), nrow = nrow(OTU_logspike_norm_tss), ncol = 7)
rownames(taxmat_logspike_norm_tss) <- rownames(OTU_logspike_norm_tss)
colnames(taxmat_logspike_norm_tss) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
TAX_logspike_norm_tss = tax_table(taxmat_logspike_norm_tss)
physeq_logspike_norm_tss = phyloseq(OTU_logspike_norm_tss, TAX_logspike_norm_tss)

#Sample file with Metadata 
conds_logspike_norm_tss = c(1,	0,	1,	1,	1,	1,	0,	1,	1,	1,	0,	0,	1,	1,	0,	1,	0,	1,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	1,	1,	1,	1,	0,	0,	0,	1,	0,	1,	0,	0,	1,	1,	0,	1,	1,	1,	1,	0,	0,	0,	1,	1,	0,	0,	1,	0,	0,	0,	0,	1,	0,	0,	1,	0,	0,	0,	0,	0,	0,	1,	0,	0,	1,	1,	0,	1,	1,	0,	0,	1,	0,	0,	1,	0,	1,	1,	1,	0,	1,	1,	0,	1,	1,	1,	0,	0)
sample_data_logspike_norm_tss = sample_data(data.frame(Group = sample(conds_logspike_norm_tss, size=nsamples(physeq_logspike_norm_tss), 
                                                                      replace=TRUE),row.names=sample_names(physeq_logspike_norm_tss),stringsAsFactors=FALSE))
#write.csv(sample_data, file = "metadata.csv", row.names = TRUE)

#physeq class
physeq_logspike_norm_tss = phyloseq(OTU_logspike_norm_tss, TAX_logspike_norm_tss, sample_data_logspike_norm_tss)
physeq_logspike_norm_tss



smp_logspike_norm_tss = row.names=sample_names(physeq_logspike_norm_tss)
#o = c("OTU1","OTU2","OTU3","OTU4","OTU5","OTU6","OTU7","OTU8","OTU9","OTU10")
#o
class.df<- data.frame(conds_logspike_norm_tss, smp_logspike_norm_tss) 

data_grp0_logspike_norm_tss = subset(class.df, conds_logspike_norm_tss == 0) # only metadata 0
data_grp1_logspike_norm_tss = subset(class.df, conds_logspike_norm_tss == 1) # only metadata 1

test_ttest_logspike_norm_tss = norm_tss_logspike
colnames_ttest_logspike_norm_tss = ifelse(conds_logspike_norm_tss == 0,"cond_zero","cond_one")
colnames_ttest_logspike_norm_tss
colnames(test_ttest_logspike_norm_tss) <- paste0(colnames_ttest_logspike_norm_tss, 1:ncol(test_ttest_logspike_norm_tss))
test_ttest_logspike_norm_tss
write.csv(test_ttest_logspike_norm_tss, file = "./first_dataset/ttest_output/test_ttest_logspike_norm_tss.csv", row.names = TRUE)

#!!!!!!!!!!!!   einfugen von "gene_name"
test_ttest_logspike_norm_tss <- read.csv("./first_dataset/ttest_output/test_ttest_logspike_norm_tss.csv")


dlong_logspike_norm_tss <- as.tbl(test_ttest_logspike_norm_tss) %>% 
  gather(key, value, -gene_name) %>% 
  mutate(group=ifelse(grepl("cond_zero",key), "control", "case"))


tt_pvalues_logspike_norm_tss_all = dlong_logspike_norm_tss %>% 
  group_by(gene_name) %>% 
  t_test(value~group)%>%
  adjust_pvalue(method = "BH")%>%
  add_significance()
tt_pvalues_logspike_norm_tss_all
write.csv(tt_pvalues_logspike_norm_tss_all, file = "./first_dataset//ttest_output/tt_pvalues_logspike_norm_tss_all.csv", row.names = FALSE)


tt_pvalues_logspike_norm_tss = dlong_logspike_norm_tss %>% 
  group_by(gene_name) %>% 
  summarise(p=t.test(value~group)$p.value)
tt_pvalues_logspike_norm_tss
write.csv(tt_pvalues_logspike_norm_tss, file = "./first_dataset/ttest_output/tt_pvalues_logspike_norm_tss.csv", row.names = FALSE)



############  RELAB #################
#build TAX file 
OTU_logspike_norm_relab = otu_table(norm_relab_logspike, taxa_are_rows = TRUE)
taxmat_logspike_norm_relab = matrix(sample(letters, 100, replace = TRUE), nrow = nrow(OTU_logspike_norm_relab), ncol = 7)
rownames(taxmat_logspike_norm_relab) <- rownames(OTU_logspike_norm_relab)
colnames(taxmat_logspike_norm_relab) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
TAX_logspike_norm_relab = tax_table(taxmat_logspike_norm_relab)
physeq_logspike_norm_relab = phyloseq(OTU_logspike_norm_relab, TAX_logspike_norm_relab)
write.table(OTU_logspike_norm_relab, file = "./first_dataset/OTU_spike.tsv", row.names=TRUE, sep="\t")

#Sample file with Metadata 
conds_logspike_norm_relab = c(1,	0,	1,	1,	1,	1,	0,	1,	1,	1,	0,	0,	1,	1,	0,	1,	0,	1,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	1,	1,	1,	1,	0,	0,	0,	1,	0,	1,	0,	0,	1,	1,	0,	1,	1,	1,	1,	0,	0,	0,	1,	1,	0,	0,	1,	0,	0,	0,	0,	1,	0,	0,	1,	0,	0,	0,	0,	0,	0,	1,	0,	0,	1,	1,	0,	1,	1,	0,	0,	1,	0,	0,	1,	0,	1,	1,	1,	0,	1,	1,	0,	1,	1,	1,	0,	0)
sample_data_logspike_norm_relab = sample_data(data.frame(Group = sample(conds_logspike_norm_relab, size=nsamples(physeq_logspike_norm_relab), 
                                                                        replace=TRUE),row.names=sample_names(physeq_logspike_norm_relab),stringsAsFactors=FALSE))
#write.csv(sample_data, file = "metadata.csv", row.names = TRUE)

#physeq class
physeq_logspike_norm_relab = phyloseq(OTU_logspike_norm_relab, TAX_logspike_norm_relab, sample_data_logspike_norm_relab)
physeq_logspike_norm_relab



smp_logspike_norm_relab = row.names=sample_names(physeq_logspike_norm_relab)
#o = c("OTU1","OTU2","OTU3","OTU4","OTU5","OTU6","OTU7","OTU8","OTU9","OTU10")
#o
class.df<- data.frame(conds_logspike_norm_relab, smp_logspike_norm_relab) 

data_grp0_logspike_norm_relab = subset(class.df, conds_logspike_norm_relab == 0) # only metadata 0
data_grp1_logspike_norm_relab = subset(class.df, conds_logspike_norm_relab == 1) # only metadata 1

test_ttest_logspike_norm_relab = norm_relab_logspike
colnames_ttest_logspike_norm_relab = ifelse(conds_logspike_norm_relab == 0,"cond_zero","cond_one")
colnames_ttest_logspike_norm_relab
colnames(test_ttest_logspike_norm_relab) <- paste0(colnames_ttest_logspike_norm_relab, 1:ncol(test_ttest_logspike_norm_relab))
test_ttest_logspike_norm_relab
write.csv(test_ttest_logspike_norm_relab, file = "./first_dataset/ttest_output/test_ttest_logspike_norm_relab.csv", row.names = TRUE)

#!!!!!!!!!!!!   einfugen von "gene_name"
test_logspike <- read.csv("./first_dataset/ttest_output/test_ttest_logspike_norm_relab.csv")


dlong_logspike <- as.tbl(test_logspike) %>% 
  gather(key, value, -gene_name) %>% 
  mutate(group=ifelse(grepl("cond_zero",key), "control", "case"))


tt_pvalues_logspike_all = dlong_logspike %>% 
  group_by(gene_name) %>% 
  t_test(value~group)%>%
  adjust_pvalue(method = "BH")%>%
  add_significance()
tt_pvalues_logspike_all
write.csv(tt_pvalues_logspike_all, file = "./first_dataset/ttest_output/tt_pvalues_logspike_norm_relab_all.csv", row.names = FALSE)


tt_pvalues_logspike = dlong_logspike %>% 
  group_by(gene_name) %>% 
  summarise(p=t.test(value~group)$p.value)
tt_pvalues_logspike
write.csv(tt_pvalues_logspike, file = "./first_dataset/ttest_output/tt_pvalues_logspike_norm_relab.csv", row.names = FALSE)


