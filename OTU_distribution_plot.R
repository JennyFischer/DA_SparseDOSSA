lognormal <- read.delim("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/OTU_lognormal.tsv", header=TRUE, row.names=NULL)
logspike <- read.delim("/Users/jf/Desktop/BA/Analysis_SparseDOSSA/OTU_spike.tsv", header=TRUE, row.names=NULL)

df_norm<-data.frame(lognormal)
pdf( "/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/plots/lognormal_OTU_distribution.pdf")
norm_dist = barplot(as.matrix(df_norm),col = rainbow(100),ylab="Counts",
        main="OTU's Distribution of relative abundance dataset") #,xlab="SampleID"
dev.off()


df_spike<-data.frame(logspike)
pdf( "/Users/jf/Desktop/BA/Analysis_SparseDOSSA/outputs/plots/logspike_OTU_distribution.pdf")
spike_dist = barplot(as.matrix(df_spike),col = rainbow(100),ylab="Counts",
        main="OTU's Distribution of relative spike dataset")
dev.off()
