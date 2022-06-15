lognormal <- read.delim("./first_dataset/OTU_lognormal.tsv", header=TRUE, row.names=NULL) #BA/Analysis_SparseDOSSA/
logspike <- read.delim("./first_dataset/OTU_logspike.tsv", header=TRUE, row.names=NULL) #BA/Analysis_SparseDOSSA/

df_norm<-data.frame(lognormal)
pdf( "./first_dataset/lognormal_OTU_distribution_exact.pdf")
norm_dist = barplot(as.matrix(df_norm[2:5]),col = rainbow(100),ylab="Counts",
        main="OTU's Distribution of relative abundance dataset") #,xlab="SampleID"
dev.off()

logspike
df_spike<-data.frame(logspike)
pdf( "./first_dataset/logspike_OTU_distribution_exact.pdf")
spike_dist = barplot(as.matrix(df_spike[2:5]),col = rainbow(100),ylab="Counts",
        main="OTU's Distribution of absolut abundance dataset")
dev.off()


