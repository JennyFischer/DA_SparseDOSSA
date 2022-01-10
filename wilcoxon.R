library(microbiomeMarker)
library(microbiome)
library(compositions)
library(phyloseq)
library(ALDEx2)
library(ggpubr)
library(rstatix)

#whole otu file from output of sparseDOSSA
otu_t <- read.delim("~/Desktop/otu_table.tsv")

#OTU file for the null (i.e., not spiked-in) absolute abundances
#first seperate python script an manual deletion of first two columns
lognormal <- read.delim("~/Desktop/lognorm.tsv", header=FALSE, row.names=NULL)
rownames(lognormal) <- paste0("OTU", 1:nrow(lognormal))
colnames(lognormal) <- paste0("Sample", 1:ncol(lognormal))


#OTU file for spiked-in absolute abundancess
#first seperate python script an manual deletion of first two columns
logspike <- read.delim("~/Desktop/spike.tsv", header=FALSE, row.names=NULL)
rownames(logspike) <- paste0("OTU", 1:nrow(logspike))
colnames(logspike) <- paste0("Sample", 1:ncol(logspike))
logspike


################################################################################
#OTU file for the null (i.e., not spiked-in) absolute abundances
#build TAX file 
OTU_lognormal = otu_table(lognormal, taxa_are_rows = TRUE)
taxmat_lognormal = matrix(sample(letters, 100, replace = TRUE), nrow = nrow(OTU_lognormal), ncol = 7)
rownames(taxmat_lognormal) <- rownames(OTU_lognormal)
colnames(taxmat_lognormal) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
TAX_lognormal = tax_table(taxmat_lognormal)
physeq_lognormal = phyloseq(OTU_lognormal, TAX_lognormal)
                            
#Sample file with Metadata 
conds_lognormal = c(1,	0,	1,	1,	1,	1,	0,	1,	1,	1,	0,	0,	1,	1,	0,	1,	0,	1,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	1,	1,	1,	1,	0,	0,	0,	1,	0,	1,	0,	0,	1,	1,	0,	1,	1,	1,	1,	0,	0,	0,	1,	1,	0,	0,	1,	0,	0,	0,	0,	1,	0,	0,	1,	0,	0,	0,	0,	0,	0,	1,	0,	0,	1,	1,	0,	1,	1,	0,	0,	1,	0,	0,	1,	0,	1,	1,	1,	0,	1,	1,	0,	1,	1,	1,	0,	0)
sample_data_lognormal = sample_data(data.frame(Group = sample(conds_lognormal, size=nsamples(physeq_lognormal), 
                                                    replace=TRUE),row.names=sample_names(physeq_lognormal),stringsAsFactors=FALSE))
write.csv(sample_data, file = "metadata.csv", row.names = TRUE)

#physeq class
physeq_lognormal = phyloseq(OTU_lognormal, TAX_lognormal, sample_data_lognormal)
physeq_lognormal



smp_lognormal = row.names=sample_names(physeq_lognormal)
#o = c("OTU1","OTU2","OTU3","OTU4","OTU5","OTU6","OTU7","OTU8","OTU9","OTU10")
#o
class.df<- data.frame(conds_lognormal, smp_lognormal) 

data_grp0_lognormal = subset(class.df, conds_lognormal == 0) # only metadata 0
data_grp1_lognormal = subset(class.df, conds_lognormal == 1) # only metadata 1

test_wlc_lognormal = lognormal
colnames_wlc_lognormal = ifelse(conds_lognormal == 0,"cond_zero","cond_one")
colnames_wlc_lognormal
colnames(test_wlc_lognormal) <- paste0(colnames_wlc_lognormal, 1:ncol(test_wlc_lognormal))
test_wlc_lognormal
write.csv(test_wlc_lognormal, file = "/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/wilcoxon_output/test_wlc_lognormal.csv", row.names = TRUE)

#einfugen von "gene_name"
test_wlc_lognormal <- read.csv("~/Desktop/BA/Analysis_SparseDOSSA/outputs/wilcoxon_output/test_wlc_lognormal.csv")


dlong_lognormal <- as.tbl(test_wlc_lognormal) %>% 
  gather(key, value, -gene_name) %>% 
  mutate(group=ifelse(grepl("cond_zero",key), "control", "case"))
#write.csv(dlong, file = "dlong.csv", row.names = TRUE)
# plot the data
#dlong %>% 
 # ggplot(aes(x=group, y=value)) +
  #geom_boxplot()

#dlong %>% 
 # with(., wilcox.test(value ~ group))

dlong_lognormal %>% 
  ggplot(aes(x=gene_name, y=value, fill=group)) +
  geom_boxplot() +
  stat_compare_means(method= "wilcox.test")

wlx_pvalues_lognormal = dlong_lognormal %>% 
  group_by(gene_name) %>% 
  summarise(p=wilcox.test(value~group)$p.value)
wlx_pvalues_lognormal
#ordnen der OTU fehlt
#write.csv(wlx_pvalues, file = "wlx_pvalues.csv", row.names = TRUE)
write.csv(wlx_pvalues_lognormal, file = "~/Desktop/BA/Analysis_SparseDOSSA/outputs/wilcoxon_output/wlx_pvalues_lognormal.csv", row.names = FALSE)

################################################################################
#OTU file for spiked-in absolute abundancess
#build TAX file 
OTU_logspike = otu_table(logspike, taxa_are_rows = TRUE)
taxmat_logspike = matrix(sample(letters, 100, replace = TRUE), nrow = nrow(OTU_logspike), ncol = 7)
rownames(taxmat_logspike) <- rownames(OTU_logspike)
colnames(taxmat_logspike) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
TAX_logspike = tax_table(taxmat_logspike)
physeq_logspike = phyloseq(OTU_logspike, TAX_logspike)

#Sample file with Metadata 
conds_logspike = c(1,	0,	1,	1,	1,	1,	0,	1,	1,	1,	0,	0,	1,	1,	0,	1,	0,	1,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	1,	1,	1,	1,	0,	0,	0,	1,	0,	1,	0,	0,	1,	1,	0,	1,	1,	1,	1,	0,	0,	0,	1,	1,	0,	0,	1,	0,	0,	0,	0,	1,	0,	0,	1,	0,	0,	0,	0,	0,	0,	1,	0,	0,	1,	1,	0,	1,	1,	0,	0,	1,	0,	0,	1,	0,	1,	1,	1,	0,	1,	1,	0,	1,	1,	1,	0,	0)
sample_data_logspike = sample_data(data.frame(Group = sample(conds_logspike, size=nsamples(physeq_logspike), 
                                                              replace=TRUE),row.names=sample_names(physeq_logspike),stringsAsFactors=FALSE))
#write.csv(sample_data, file = "metadata.csv", row.names = TRUE)

#physeq class
physeq_logspike = phyloseq(OTU_logspike, TAX_logspike, sample_data_logspike)
physeq_logspike



smp_logspike = row.names=sample_names(physeq_logspike)
#o = c("OTU1","OTU2","OTU3","OTU4","OTU5","OTU6","OTU7","OTU8","OTU9","OTU10")
#o
class.df<- data.frame(conds_logspike, smp_logspike) 

data_grp0_logspike = subset(class.df, conds_logspike == 0) # only metadata 0
data_grp1_logspike = subset(class.df, conds_logspike == 1) # only metadata 1

test_wlc_logspike = logspike
colnames_wlc_logspike = ifelse(conds_logspike == 0,"cond_zero","cond_one")
colnames_wlc_logspike
colnames(test_wlc_logspike) <- paste0(colnames_wlc_logspike, 1:ncol(test_wlc_logspike))
test_wlc_logspike
write.csv(test_wlc_logspike, file = "~/Desktop/BA/Analysis_SparseDOSSA/outputs/wilcoxon_output/test_wlc_logspike.csv", row.names = TRUE)

#einfugen von "gene_name"
test_wlc_logspike <- read.csv("~/Desktop/BA/Analysis_SparseDOSSA/outputs/wilcoxon_output/test_wlc_logspike.csv")


dlong_logspike <- as.tbl(test_wlc_logspike) %>% 
  gather(key, value, -gene_name) %>% 
  mutate(group=ifelse(grepl("cond_zero",key), "control", "case"))
#write.csv(dlong, file = "dlong.csv", row.names = TRUE)
# plot the data
#dlong %>% 
# ggplot(aes(x=group, y=value)) +
#geom_boxplot()

#dlong %>% 
# with(., wilcox.test(value ~ group))

dlong_logspike %>% 
  ggplot(aes(x=gene_name, y=value, fill=group)) +
  geom_boxplot() +
  stat_compare_means(method= "wilcox.test")

wlx_pvalues_logspike = dlong_logspike %>% 
  group_by(gene_name) %>% 
  summarise(p=wilcox.test(value~group)$p.value)
wlx_pvalues_logspike
#ordnen der OTU fehlt
#write.csv(wlx_pvalues, file = "wlx_pvalues.csv", row.names = TRUE)
write.csv(wlx_pvalues_logspike, file = "~/Desktop/BA/Analysis_SparseDOSSA/outputs/wilcoxon_output/wlx_pvalues_logspike.csv", row.names = FALSE)

