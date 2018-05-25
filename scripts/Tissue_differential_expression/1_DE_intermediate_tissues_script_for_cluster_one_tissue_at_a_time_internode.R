# Aim is to calculate the number of DE genes between intermediate tissues in Bayer all vs all 
# 18.7.2017
# Philippa Borrill


#### prep data ######
# setwd
setwd("/nbi/group-data/NBI/Cristobal-Uauy/expression_browser/WGA/Bayer_analysis/DE")

library("matrixStats")
library("DESeq2")

# read in counts
counts <- read.csv(file="/nbi/group-data/NBI/Cristobal-Uauy/expression_browser/WGA/WGCNA/expressionValuesPerGene/Development_count.tsv", header=T, sep = "\t")
head(counts)
dim(counts)

# read in metadata
metadata <- read.csv(file="/nbi/group-data/NBI/Cristobal-Uauy/expression_browser/WGA/WGCNA/HC_only/WGCNA_metadata_cluster_with_intermed_coarse.csv", header=T)
head(metadata)
metadata[1:4,1:4]
colnames(metadata)
metadata <- metadata[metadata$Provider == "M. Davey BCS",] 
dim(metadata)

# now need to make metadata into the same order as the counts dataframe
df_counts_order <- as.data.frame(colnames(counts))
head(df_counts_order)

ordered_metadata <- merge(df_counts_order, metadata[,c(4,22)], by.x="colnames(counts)", by.y="Sample.IDs")
head(ordered_metadata)
colnames(ordered_metadata) <- c("Sample.IDs","Intermediate")

# check they are in the same order
ordered_metadata$Sample.IDs == colnames(counts)

###### now load into DESeq2 #####

# get rid of any genes which don't pass tpm filter

tpms <- read.csv(file="/nbi/group-data/NBI/Cristobal-Uauy/expression_browser/WGA/WGCNA/expressionValuesPerGene/Development_tpm.tsv", header=T, sep = "\t")
dim(tpms)
dim(counts)
head(colnames(tpms))
head(colnames(counts))
head(rownames(tpms))
head(rownames(counts))

# Filter count data to only keep genes expressed >0.5 tpm in at least 1 % of samples
tpm_threshold <- 0.5 

# set what is 1% of samples
perc_samples <- nrow(ordered_metadata)/100


#checked the filtering works as expected 15-05-2017
count_filt <- counts[rowCounts(as.matrix(tpms>tpm_threshold))>perc_samples,] # select only rows which have expr >0.5 tpm in > 1 % of samples (2.09)
nrow(counts)
nrow(count_filt)

# get rid of LC genes
count_filt[1:4,1:4]
dim(count_filt)
length(grepl("LC",rownames(count_filt)))
length(!grepl("LC",rownames(count_filt)))
countsHC <- count_filt[(!grepl("LC",rownames(count_filt))),]

countsHC[1:4,1:4]
dim(countsHC)

head(sort(rowSums(countsHC)))

# round numbers
bckCountTable <- round(countsHC)
head(bckCountTable)

# get rid of any rows which have a sum of 0
bckCountTable <- bckCountTable[rowSums(bckCountTable) != 0,]
dim(bckCountTable)

# check metadata and columns are in the same order
ordered_metadata$Sample.IDs == colnames(bckCountTable)

samples <- data.frame(row.names=colnames(bckCountTable), condition=as.factor(as.character(ordered_metadata$Intermediate)))

samples

bckCDS <- DESeqDataSetFromMatrix(countData = bckCountTable, colData=samples, design=~condition)
bckCDS

# check conditions are ok
colData(bckCDS)$condition

#make sure the seedling aerial tissues is used as the reference condition :
colData(bckCDS)$condition <- relevel(colData(bckCDS)$condition, "seedling aerial tissues")

# run DeSeq2
bckCDS_1 <- DESeq(bckCDS)


## now want to differential expression for all tissues vs all
unique(samples$condition) # these are the tissues I have


#set tissue to compare to
j <- as.character(unique(samples$condition)[6])
j

# for all conditions except the already selected condition
for (i in unique(samples$condition)[unique(samples$condition)!= j]) { 
  
      bck_res <- results(bckCDS_1,contrast=c("condition",i,j))
  
  # sort results on padj
  ordered_res <- bck_res[order(bck_res$padj),]
  head(ordered_res)
  
  ordered_res_na.rm <- na.omit(ordered_res)
  
  ordered_res_na.rm_padj0.05 <- ordered_res_na.rm[ordered_res_na.rm$padj<0.05,] 
  
  # output ordered_res to csv
  write.csv(ordered_res,file=paste(i,"_vs_",j,"_results_padj_0.05.csv"))

}

