
# using http://www.statmethods.net/advstats/cluster.html to run hclust on tpm data
# Philippa Borrill
# 5.6.2017

library("dendextend")
library("pvclust")

setwd("Y:\\expression_browser\\WGA\\hclust\\")
#setwd("/nbi/Research-Groups/NBI/Cristobal-Uauy/expression_browser/WGA/hclust")

tpm_data <- read.delim(file="Development_mean_tpm_per_chr_group.txt", sep= "\t", header=T)

tpm_data[1:4,1:4]

setwd("Y:\\expression_browser\\WGA\\hclust\\for_manuscript")
#setwd("/nbi/Research-Groups/NBI/Cristobal-Uauy/expression_browser/WGA/hclust/for_manuscript_cluster")

# used log2 data
tpm_data_plus1 <- tpm_data+1
tpm_data_plus1[1:4,1:4]

rownames(tpm_data_plus1) <- tpm_data_plus1[,1]
tpm_data_plus1 <- tpm_data_plus1[,-1]
tpm_data_plus1[1:4,1:4]

mydata <- log2(tpm_data_plus1)
mydata[1:4,1:4]
dim(mydata)

# transpose data
mydata <- t(mydata)
mydata[1:4,1:4]

# Average Hierarchical Clustering (clusters rows)
d <- dist(mydata, method = "euclidean") # distance matrix
fit <- hclust(d, method="average")

# now plot with colours
# now plot dendrogram with colours according to homoeologue

rownames(mydata)
genomes <- rep(c("A", "B", "D"),each= 70 )
genomes
# convert sampleTree to dendrogram
dend <- as.dendrogram(fit)
head(genomes)
head(genomes[(order.dendrogram(dend))])
colorCodes <- c(A= "#579D1C", B = "#4B1F6F", D = "#FF950E")
labels_colors(dend) <- colorCodes[(genomes[(order.dendrogram(dend))])]
labels_colors(dend)


# now plot dendrogram with colours according to homoeologue with branch lengths and high level tissue
tissue <- read.csv("high_level_tissue.csv")
head(tissue)
unique(tissue$High_level_tissue)
unique(tissue$Intermediate)

# now make column with right colours
tissue$colour <- tissue$High_level_tissue
tissue$colour <- gsub("grain","#dfc27d",tissue$colour)
tissue$colour <- gsub("leaves","#7fbc41",tissue$colour)
tissue$colour <- gsub("roots","#a6611a",tissue$colour)
tissue$colour <- gsub("spike","#01665e",tissue$colour)
unique(tissue$colour)



## now try bootstrapping using 1,000 gene sub-samples
mydata[1:4,1:4]
 # select 1000 random columns
sample1 <- mydata[,sample(1:ncol(mydata),1000,replace=FALSE)]
sample1[1:4,1:4]
dim(sample1)

set.seed(13134)
result <- pvclust(t(sample1), method.dist ="euclidean", method.hclust = "average", nboot=10)
dend_res <- as.dendrogram(result)
colorCodes <- c(A= "#579D1C", B = "#4B1F6F", D = "#FF950E")
labels_colors(dend_res) <- colorCodes[(genomes[(order.dendrogram(dend_res))])]
hung_dend <- hang.dendrogram(dend_res, hang=0.01)


# now calculate p-values

# run with 1,000 bootstraps
set.seed(13134)
result <- pvclust(t(mydata), method.dist ="euclidean", method.hclust = "average", nboot=1000)

saveRDS(result, file="pvclust_result_1000_bootstraps.RDS") # save result

head(result)

dend_res <- as.dendrogram(result)
colorCodes <- c(A= "#579D1C", B = "#4B1F6F", D = "#FF950E")
labels_colors(dend_res) <- colorCodes[(genomes[(order.dendrogram(dend_res))])]

pdf(file="hclust_Development_mean_tpm_per_chr_group_coloured_by_homoeologue_with_lengths_pvalues_high_level_tissue_1000bootstraps.pdf", width=30, height=15)
par(cex = 0.6);
par(mar = c(20,4,2,0))
dend_res %>% hang.dendrogram(0.01) %>%
  plot(main= "Cluster dendrogram with AU/BP values (%)") 
colored_bars(colors=tissue$colour, dend=hung_dend, rowLabels = "High level tissue" )
result %>% text
dev.off()


