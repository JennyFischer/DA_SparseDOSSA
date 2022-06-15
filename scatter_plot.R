library(ggplot2)
library(dplyr)

#lognormal <- read.delim("~/Desktop/OTU_lognormal.tsv", header=TRUE, row.names=NULL)
lognormal <- read.delim("./first_dataset/combined_data_null.tsv", header=TRUE, row.names=NULL)

lognormal <- as.data.frame(t(lognormal))
lognormal =lognormal[-c(0,1,2,3),]
write.table(lognormal, file = "./first_dataset/testlogn.tsv", row.names=TRUE, sep="\t")
lognormal= read.table(file = './first_dataset/testlogn.tsv', sep = '\t', header = TRUE)
lognormal
#spike
logspike <- read.delim("./first_dataset/combined_data_spike.tsv", header=TRUE, row.names=NULL)
logspike  <- as.data.frame(t(logspike))
logspike  =logspike [-c(0,1,2,3),]
write.table(logspike, file = "./first_dataset/testlogs.tsv", row.names=TRUE, sep="\t")
logspike= read.table(file = './first_dataset/testlogs.tsv', sep = '\t', header = TRUE)


my_list1 <- as.vector(lognormal$Test) 
my_list1

l = my_list1[1:100]
l 

my_list2 <- as.vector(lognormal$V2) 
my_list2
p = my_list2[1:100]
p

my_list3 <- as.vector(logspike$V2) 
my_list3
d = my_list3[1:100]
d

class(t)
a.country <- rep(c("Relative Abundance"), 100)
a.value   <- as.numeric(p)#c( 274,   0, 479,   0 )) 
#a.value   <- as.numeric(c(1, 4, 3, 5))
a.year    <-l# c("OTU1", "OTU2", "OTU3", "OTU4")
dtframe.a <- data.frame(a.country, a.year, a.value)

#Create dataframe B
b.country <- rep(c("Absolute Abundance"), 100)
b.value   <- as.numeric(d)
b.year    <- l #c("OTU1", "OTU2", "OTU3", "OTU4")
dtframe.b <- data.frame(b.country, b.year, b.value)
dtframe.b 
# Use ggplot2 to plot data from 2 dataframes
require(ggplot2)


OTUs =a.year
Counts=a.value
Abundances= a.country
ggplot() +
  geom_point(data=dtframe.a, aes(OTUs, Counts, color= Abundances)) +
  geom_point(data=dtframe.b, aes(b.year, b.value, color= b.country))+
  theme(axis.text.x=element_blank())



library(reshape2)
melted_logn <- melt(lognormal, id = "Test")
melted_logn
melted_logs <- melt(logspike, id = "Test")
melted_logs


my_list1 <- as.vector(lognormal$Test) 
my_list1

l = my_list1[1:100]
l 
l = rep(c(l), 10)
l
my_list2 <- as.vector(melted_logn$value) 
my_list2
p = my_list2[1:1000]
p

my_list3 <- as.vector(melted_logs$value) 
my_list3
d = my_list3[1:1000]
d

class(t)
a.country <- rep(c("Relative Abundance"), 1000)
a.value   <- as.numeric(p)#c( 274,   0, 479,   0 )) 
#a.value   <- as.numeric(c(1, 4, 3, 5))
a.year    <-l# c("OTU1", "OTU2", "OTU3", "OTU4")
dtframe.a <- data.frame(a.country, a.year, a.value)

#Create dataframe B
b.country <- rep(c("Absolute Abundance"), 1000)
b.value   <- as.numeric(d)
b.year    <- l #c("OTU1", "OTU2", "OTU3", "OTU4")
dtframe.b <- data.frame(b.country, b.year, b.value)
dtframe.b 
# Use ggplot2 to plot data from 2 dataframes
require(ggplot2)


OTUs =a.year
Counts=a.value
Abundances= a.country
ggplot() +
  geom_point(data=dtframe.a, aes(OTUs, Counts, color= Abundances)) +
  geom_point(data=dtframe.b, aes(b.year, b.value, color= b.country))+
  theme(axis.text.x=element_blank())




o =rowSums(lognormal)
o
length(o)
length(which(o==0))
r =rowSums(logspike)
length(which(r==0))
r
length(r)
