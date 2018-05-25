# get expression data for nullitetras for 1:1:1 triads
# Philippa Borrill
# 28.9.2017



##### Load the data #####

library("matrixStats")

setwd("Y:\\expression_browser\\WGA\\nullitetras\\")

loadValuesFromExperiment<-function(metadata, folder, unit="tpm", values=c("Development")){
  metadata$Sample.IDs <- gsub("-",".",metadata$Sample.IDs)
  #  print(metadata$Sample.IDs)
  v<-values[1]
  v<-gsub(" ","_",v)
  path<-paste0(folder,"\\",v,"_",unit,".tsv")
  ret<-read.table(path, row.names = 1, header= TRUE)
  #  print(v)
  for(i in 2:length(values)){
    v<-values[i]
    v<-gsub(" ","_",v)
    #   print(v)
    path<-paste0(folder,"\\",v,"_",unit,".tsv")
    tmp<-read.table(path, row.names = 1, header= TRUE)
    ret<-cbind(ret,tmp)
  }
  #  print(colnames(ret))
  md<-metadata[metadata$Sample.IDs%in%colnames(ret),]
  ret<-ret[,as.character(md$Sample.IDs),]
  list(ret,md)
}

folder<-"Y:\\expression_browser\\collaborators\\kallisto\\1145_samples\\RefSeq_1.0\\ByGene"

metadata_file<- "Y:\\expression_browser\\WGA\\nullitetras\\Metadata_Oct2_with_missing_samples.tab"

metadata <- read.csv(metadata_file, row.names = 1, sep="\t")
metadata[1:10,1:4]
dim(metadata)
# filter out any samples we don't want to include

metadata <-  metadata[metadata[,"study_title"] == "SRP028357",]

metadata[1:10,1:4]
print("nrow of metadata")
print(nrow(metadata))

#tpms  <-loadValuesFromExperiment(metadata, folder, unit="tpm",  values=unique(metadata$study_title)) # for some reason get the "NA" problem

tpms  <- read.csv(file="Y:\\expression_browser\\collaborators\\kallisto\\1145_samples\\RefSeq_1.0\\ByGene\\SRP028357_tpm.tsv", header=T, sep = "\t")

tpms[1:4,1:4]
dim(tpms)
head(colnames(tpms))
head(rownames(tpms))

print("dimension of tpms")
print(dim(tpms))

# get rid of LC names
tpmsHC <- tpms[(!grepl("LC",tpms$gene)),]
print(dim(tpmsHC))

# now want to add in metadata
tpmsHC[1:4,1:4]

metadata_useful <- (metadata)[,c(3,7,10)]
#metadata_useful
head(metadata_useful)
dim(metadata_useful)

metadata_useful[metadata_useful$Variety =="Chinese Spring, N1BT1A",]

# convert tpmsHC to long format
library("reshape2")
melted_tpmsHC <- melt(tpmsHC, id.vars=c("gene"))
head(melted_tpmsHC)

# add in metadata
head(metadata_useful)
merged_melted_tpmsHC <- merge(melted_tpmsHC, metadata_useful, by.x="variable", by.y="Sample.IDs")
head(merged_melted_tpmsHC)
colnames(merged_melted_tpmsHC)[1:3] <- c("sample","gene","tpm")
head(merged_melted_tpmsHC)
tail(merged_melted_tpmsHC)
merged_melted_tpmsHC$tpm <- as.numeric(merged_melted_tpmsHC$tpm) # convert to numeric


# now calculate average per gene per variety per tissue
library("plyr")
mean_tpmsHC <- ddply(merged_melted_tpmsHC,
                               c("Variety","gene","High.level.tissue"),
                               summarise,
                               mean= mean(tpm))
head(mean_tpmsHC)

# read in homoeologue info table
hom.df <- read.csv(file="Y:\\expression_browser\\WGA\\data_tables\\wheat.homeolog_groups.release.nonTE.csv", header=T)
head(hom.df)

hom.df_1_1_1 <- hom.df[hom.df$cardinality_abs == "1:1:1" & hom.df$HC.LC =="HC-only",]
head(hom.df_1_1_1)
dim(hom.df_1_1_1)

# need to select triads with expression > 0.5 tpm in Chinese Spring
tpm_threshold <- 0.5 
head(mean_tpmsHC)

# first do this for leaf:
mean_tpmsHC_CS_leaf <- mean_tpmsHC[mean_tpmsHC$Variety == "Chinese Spring" & mean_tpmsHC$High.level.tissue == "leaves/shoots",]
head(mean_tpmsHC_CS_leaf)

leaf_CS_triad_expr <- merge(hom.df_1_1_1[,7:9],mean_tpmsHC_CS_leaf, by.x = "A", by.y="gene" )
head(leaf_CS_triad_expr)
leaf_CS_triad_expr <- merge(leaf_CS_triad_expr,mean_tpmsHC_CS_leaf, by.x = "B", by.y="gene", suffixes = c(".A",".B") )
head(leaf_CS_triad_expr)
leaf_CS_triad_expr <- merge(leaf_CS_triad_expr,mean_tpmsHC_CS_leaf, by.x = "D", by.y="gene")
head(leaf_CS_triad_expr)
leaf_CS_triad_expr <- leaf_CS_triad_expr[,c(3,2,1,6,9,12)]
colnames(leaf_CS_triad_expr)[6] <- "mean.D"
head(leaf_CS_triad_expr)
leaf_CS_triad_expr$sum <- (leaf_CS_triad_expr$mean.A +leaf_CS_triad_expr$mean.B + leaf_CS_triad_expr$mean.D )
head(leaf_CS_triad_expr)

leaf_CS_triad_expr_over0.5 <- leaf_CS_triad_expr[leaf_CS_triad_expr$sum > tpm_threshold,]
dim(leaf_CS_triad_expr)
dim(leaf_CS_triad_expr_over0.5)

# do some triads not have a 1A, 1B and 1D?
dim(leaf_CS_triad_expr_over0.5[grepl("CS1A",leaf_CS_triad_expr_over0.5$A)& 
                                 grepl("CS1B",leaf_CS_triad_expr_over0.5$B) & 
                                 grepl("CS1D",leaf_CS_triad_expr_over0.5$D),])

# yes some triads did not have a 1A, 1B and 1D therefore keep only triads with 1A, 1B and 1D

leaf_CS_triad_expr_over0.5 <- (leaf_CS_triad_expr_over0.5[grepl("CS1A",leaf_CS_triad_expr_over0.5$A)& 
                                                            grepl("CS1B",leaf_CS_triad_expr_over0.5$B) & 
                                                            grepl("CS1D",leaf_CS_triad_expr_over0.5$D),])
head(leaf_CS_triad_expr_over0.5)
dim(leaf_CS_triad_expr_over0.5)

leaf_triads <- c(as.character(leaf_CS_triad_expr_over0.5$A), as.character(leaf_CS_triad_expr_over0.5$B), as.character(leaf_CS_triad_expr_over0.5$D))
head(leaf_triads)
length(leaf_triads)
tail(leaf_triads)

# second do this for roots:
mean_tpmsHC_CS_root <- mean_tpmsHC[mean_tpmsHC$Variety == "Chinese Spring" & mean_tpmsHC$High.level.tissue == "roots",]
head(mean_tpmsHC_CS_root)

root_CS_triad_expr <- merge(hom.df_1_1_1[,7:9],mean_tpmsHC_CS_root, by.x = "A", by.y="gene" )
head(root_CS_triad_expr)
root_CS_triad_expr <- merge(root_CS_triad_expr,mean_tpmsHC_CS_root, by.x = "B", by.y="gene", suffixes = c(".A",".B") )
head(root_CS_triad_expr)
root_CS_triad_expr <- merge(root_CS_triad_expr,mean_tpmsHC_CS_root, by.x = "D", by.y="gene")
head(root_CS_triad_expr)
root_CS_triad_expr <- root_CS_triad_expr[,c(3,2,1,6,9,12)]
colnames(root_CS_triad_expr)[6] <- "mean.D"
head(root_CS_triad_expr)
root_CS_triad_expr$sum <- (root_CS_triad_expr$mean.A +root_CS_triad_expr$mean.B + root_CS_triad_expr$mean.D )
head(root_CS_triad_expr)

root_CS_triad_expr_over0.5 <- root_CS_triad_expr[root_CS_triad_expr$sum > tpm_threshold,]
dim(root_CS_triad_expr)
head(root_CS_triad_expr_over0.5)
dim(root_CS_triad_expr_over0.5)

# do some triads not have a 1A, 1B and 1D?
dim(root_CS_triad_expr_over0.5[grepl("CS1A",root_CS_triad_expr_over0.5$A)& 
                                 grepl("CS1B",root_CS_triad_expr_over0.5$B) & 
                                 grepl("CS1D",root_CS_triad_expr_over0.5$D),])

# yes some triads did not have a 1A, 1B and 1D therefore keep only triads with 1A, 1B and 1D

root_CS_triad_expr_over0.5 <- (root_CS_triad_expr_over0.5[grepl("CS1A",root_CS_triad_expr_over0.5$A)& 
                              grepl("CS1B",root_CS_triad_expr_over0.5$B) & 
                              grepl("CS1D",root_CS_triad_expr_over0.5$D),])
head(root_CS_triad_expr_over0.5)
dim(root_CS_triad_expr_over0.5)

root_triads <- c(as.character(root_CS_triad_expr_over0.5$A), as.character(root_CS_triad_expr_over0.5$B), as.character(root_CS_triad_expr_over0.5$D))
head(root_triads)
length(root_triads)
tail(root_triads)


# Filter  data to only keep genes which are in 1:1:1 triads with triad sum expression >0.5 tpm 
head(mean_tpmsHC)

genes_over_0.5WT_leaves <- mean_tpmsHC[mean_tpmsHC$Variety =="Chinese Spring"  & 
                                    mean_tpmsHC$High.level.tissue == "leaves/shoots" &
                                    mean_tpmsHC$gene %in% leaf_triads,]
head(genes_over_0.5WT_leaves)

genes_over_0.5WT_roots <- mean_tpmsHC[mean_tpmsHC$Variety =="Chinese Spring" &
                                        mean_tpmsHC$High.level.tissue == "roots" &
                                        mean_tpmsHC$gene %in% root_triads,]
head(genes_over_0.5WT_roots)

length(unique(mean_tpmsHC$gene))
length(unique(genes_over_0.5WT_leaves$gene))
length(unique(genes_over_0.5WT_roots$gene))

# now separate out only chromosome1 or 5
# chr1
genes_over_0.5WT_leaves_chr1A <- genes_over_0.5WT_leaves[grepl("CS1A",genes_over_0.5WT_leaves$gene),]
head(genes_over_0.5WT_leaves_chr1A)
tail(genes_over_0.5WT_leaves_chr1A)

genes_over_0.5WT_leaves_chr1B <- genes_over_0.5WT_leaves[grepl("CS1B",genes_over_0.5WT_leaves$gene),]

genes_over_0.5WT_leaves_chr1D <- genes_over_0.5WT_leaves[grepl("CS1D",genes_over_0.5WT_leaves$gene),]

genes_over_0.5WT_roots_chr1A <- genes_over_0.5WT_roots[grepl("CS1A",genes_over_0.5WT_roots$gene),]

genes_over_0.5WT_roots_chr1B <- genes_over_0.5WT_roots[grepl("CS1B",genes_over_0.5WT_roots$gene),]

genes_over_0.5WT_roots_chr1D <- genes_over_0.5WT_roots[grepl("CS1D",genes_over_0.5WT_roots$gene),]


# make a dataframe for each heatmap to plot using only genes expressed in correct tissue on correct chr in WT

levels(droplevels(mean_tpmsHC$Variety))

# chr1A leaves
long_data_chr1A_leaves <- mean_tpmsHC[!grepl("N5",mean_tpmsHC$Variety) & mean_tpmsHC$High.level.tissue == "leaves/shoots" &
                                        mean_tpmsHC$gene %in% genes_over_0.5WT_leaves_chr1A$gene,]
head(long_data_chr1A_leaves)
tail(long_data_chr1A_leaves)
levels(droplevels(long_data_chr1A_leaves$Variety))

library(tidyr)
wide_data_chr1A_leaves <- spread(long_data_chr1A_leaves[,c(1,2,4)], key = Variety, value = mean)
head(wide_data_chr1A_leaves)

rownames(wide_data_chr1A_leaves) <- wide_data_chr1A_leaves$gene
wide_data_chr1A_leaves <- wide_data_chr1A_leaves[,-1]
matrix_wide_data_chr1A_leaves <- as.matrix(wide_data_chr1A_leaves)
head(matrix_wide_data_chr1A_leaves)

as.numeric(colMeans(matrix_wide_data_chr1A_leaves))[1]
## now set output directory
setwd("Y:\\expression_browser\\WGA\\nullitetras\\1_1_1_triads")


## now do for all chr combinations ###
## for leaves and roots
library(tidyr)
library("NMF")
library(gplots)
nmf.options(grid.patch=TRUE)

chrs <- c("1A","1B","1D")
tissues <- c("leaves","roots")

nullitetra_summary <- data.frame(chromosome=character(), tissue=character(), Chinese_Spring=numeric(),
                                 NulliATetraB=numeric(), NulliATetraD=numeric(), NulliBTetraA=numeric(),
                                 NulliBTetraD=numeric(), NulliDTetraA=numeric(), NulliDTetraB=numeric(), stringsAsFactors=F)
nullitetra_summary

for (chr in chrs) { 
  for (tissue in tissues) {

   chrom_inc <- substr(chr,1,1)
   
   
long_data <- mean_tpmsHC[(grepl(paste0("N",chrom_inc),mean_tpmsHC$Variety) | mean_tpmsHC$Variety == "Chinese Spring") & grepl(tissue, mean_tpmsHC$High.level.tissue) &
                           mean_tpmsHC$gene %in% get(paste0("genes_over_0.5WT_",tissue,"_chr",chr))$gene,]

wide_data <- spread(long_data[,c(1,2,4)], key = Variety, value = mean)
head(wide_data)

rownames(wide_data) <- wide_data$gene
wide_data <- wide_data[,-1]
matrix_wide_data <- as.matrix(wide_data)


to_add <- data.frame(chromosome= as.character(chr), tissue = as.character(tissue), Chinese_Spring= as.numeric(colMeans(matrix_wide_data))[1],
            NulliATetraB=as.numeric(colMeans(matrix_wide_data))[2], NulliATetraD = as.numeric(colMeans(matrix_wide_data))[3],
           NulliBTetraA=as.numeric(colMeans(matrix_wide_data))[4],NulliBTetraD=as.numeric(colMeans(matrix_wide_data))[5], 
           NulliDTetraA=as.numeric(colMeans(matrix_wide_data))[6], NulliDTetraB=as.numeric(colMeans(matrix_wide_data))[7],
           stringsAsFactors=F)

nullitetra_summary <- rbind(nullitetra_summary, (to_add))
  }
}

nullitetra_summary

write.csv(file="1_1_1_triads_nullitetra_mean_tpm.csv",nullitetra_summary)


### 20.2.2018 ####
# now we want to know the distribution across all triads
# already have the CS tpm for the triads:
head(leaf_CS_triad_expr_over0.5)
dim(leaf_CS_triad_expr_over0.5)

# need to add columns for expression in tpm of A, B and D gene in nulliA, nulliB and nulliD respectively
# get this for each chr and tissue at a time:

# chr1A leaves
long_data_chr1A_leaves <- mean_tpmsHC[grepl("N1A",mean_tpmsHC$Variety) & mean_tpmsHC$High.level.tissue == "leaves/shoots" &
                                        mean_tpmsHC$gene %in% genes_over_0.5WT_leaves_chr1A$gene,]
head(long_data_chr1A_leaves)
tail(long_data_chr1A_leaves)
levels(droplevels(long_data_chr1A_leaves$Variety))

# calculate mean expr per gene in N1A 
leaf_N1A <- ddply(long_data_chr1A_leaves,
                     c("gene"),
                     summarise,
                     mean= mean(mean))
head(leaf_N1A)
dim(leaf_N1A)
# chr1B leaves
long_data_chr1B_leaves <- mean_tpmsHC[grepl("N1B",mean_tpmsHC$Variety) & mean_tpmsHC$High.level.tissue == "leaves/shoots" &
                                        mean_tpmsHC$gene %in% genes_over_0.5WT_leaves_chr1B$gene,]
head(long_data_chr1B_leaves)
tail(long_data_chr1B_leaves)
levels(droplevels(long_data_chr1B_leaves$Variety))

# calculate mean expr per gene in N1B 
leaf_N1B <- ddply(long_data_chr1B_leaves,
                  c("gene"),
                  summarise,
                  mean= mean(mean))
head(leaf_N1B)

# chr1D leaves
long_data_chr1D_leaves <- mean_tpmsHC[grepl("N1D",mean_tpmsHC$Variety) & mean_tpmsHC$High.level.tissue == "leaves/shoots" &
                                        mean_tpmsHC$gene %in% genes_over_0.5WT_leaves_chr1D$gene,]
head(long_data_chr1D_leaves)
tail(long_data_chr1D_leaves)
levels(droplevels(long_data_chr1D_leaves$Variety))

# calculate mean expr per gene in N1D 
leaf_N1D <- ddply(long_data_chr1D_leaves,
                  c("gene"),
                  summarise,
                  mean= mean(mean))
head(leaf_N1D)

# chr1A roots
long_data_chr1A_roots <- mean_tpmsHC[grepl("N1A",mean_tpmsHC$Variety) & mean_tpmsHC$High.level.tissue == "roots" &
                                        mean_tpmsHC$gene %in% genes_over_0.5WT_roots_chr1A$gene,]
head(long_data_chr1A_roots)
tail(long_data_chr1A_roots)
levels(droplevels(long_data_chr1A_roots$Variety))

# calculate mean expr per gene in N1A 
root_N1A <- ddply(long_data_chr1A_roots,
                  c("gene"),
                  summarise,
                  mean= mean(mean))
head(root_N1A)
dim(root_N1A)
# chr1B roots
long_data_chr1B_roots <- mean_tpmsHC[grepl("N1B",mean_tpmsHC$Variety) & mean_tpmsHC$High.level.tissue == "roots" &
                                       mean_tpmsHC$gene %in% genes_over_0.5WT_roots_chr1B$gene,]
head(long_data_chr1B_roots)
tail(long_data_chr1B_roots)
levels(droplevels(long_data_chr1B_roots$Variety))

# calculate mean expr per gene in N1B 
root_N1B <- ddply(long_data_chr1B_roots,
                  c("gene"),
                  summarise,
                  mean= mean(mean))
head(root_N1B)


# chr1D roots
long_data_chr1D_roots <- mean_tpmsHC[grepl("N1D",mean_tpmsHC$Variety) & mean_tpmsHC$High.level.tissue == "roots" &
                                       mean_tpmsHC$gene %in% genes_over_0.5WT_roots_chr1D$gene,]
head(long_data_chr1D_roots)
tail(long_data_chr1D_roots)
levels(droplevels(long_data_chr1D_roots$Variety))

# calculate mean expr per gene in N1D 
root_N1D <- ddply(long_data_chr1D_roots,
                  c("gene"),
                  summarise,
                  mean= mean(mean))
head(root_N1D)


# now merge together into 1 data set for leaf and 1 for root:
## leaf ####
head(leaf_CS_triad_expr_over0.5)
head(leaf_N1A)
head(leaf_N1B)
head(leaf_N1D)

leaf_triad_expr_in_nullis <- merge(leaf_CS_triad_expr_over0.5,leaf_N1A, by.x = "A", by.y="gene")
head(leaf_triad_expr_in_nullis)
leaf_triad_expr_in_nullis <- merge(leaf_triad_expr_in_nullis,leaf_N1B, by.x = "B", by.y="gene", suffixes = c(".nulliA",".nulliB"))
head(leaf_triad_expr_in_nullis)
leaf_triad_expr_in_nullis <- merge(leaf_triad_expr_in_nullis,leaf_N1D, by.x = "D", by.y="gene")
head(leaf_triad_expr_in_nullis)
dim(leaf_triad_expr_in_nullis)
colnames(leaf_triad_expr_in_nullis)[10] <- "mean.nulliD"
leaf_triad_expr_in_nullis <- leaf_triad_expr_in_nullis[,c(3,2,1,4,5,6,8,9,10)]
head(leaf_triad_expr_in_nullis)
colnames(leaf_triad_expr_in_nullis) <- c("A","B","D","CS_tpmA","CS_tpmB","CS_tpmD","tpm_nulliA","tpm_nulliB","tpm_nulliD")
head(leaf_triad_expr_in_nullis)

# add column with % mismapping in nulliA, nulliB, nulliD and average
leaf_triad_expr_in_nullis$perc_nulliA_mismap <- leaf_triad_expr_in_nullis$tpm_nulliA/leaf_triad_expr_in_nullis$CS_tpmA *100
leaf_triad_expr_in_nullis$perc_nulliB_mismap <- leaf_triad_expr_in_nullis$tpm_nulliB/leaf_triad_expr_in_nullis$CS_tpmB *100
leaf_triad_expr_in_nullis$perc_nulliD_mismap <- leaf_triad_expr_in_nullis$tpm_nulliD/leaf_triad_expr_in_nullis$CS_tpmD *100
leaf_triad_expr_in_nullis$av_perc_nulli_mismap <- (leaf_triad_expr_in_nullis$perc_nulliA_mismap +
                                                     leaf_triad_expr_in_nullis$perc_nulliB_mismap + 
                                                     leaf_triad_expr_in_nullis$perc_nulliD_mismap)/3
head(leaf_triad_expr_in_nullis)

# need to remove any triads where A, or B or D is 0 in CS (otherwise I get an infinite fold change in the nulli line)
leaf_triad_expr_in_nullis_no_zeroCS <- leaf_triad_expr_in_nullis[leaf_triad_expr_in_nullis$CS_tpmA != 0 &
                                                                   leaf_triad_expr_in_nullis$CS_tpmB != 0 &
                                                                   leaf_triad_expr_in_nullis$CS_tpmD != 0,]
# want to add the triad ID Ricardo uses:
head(hom.df_1_1_1)
head(leaf_triad_expr_in_nullis_no_zeroCS)
leaf_triad_expr_in_nullis_no_zeroCS <- merge(leaf_triad_expr_in_nullis_no_zeroCS, hom.df_1_1_1[,c(1,7)], by.x = "A", by.y="A")
head(leaf_triad_expr_in_nullis_no_zeroCS)
dim(leaf_triad_expr_in_nullis_no_zeroCS)

# save the data 
write.csv(file="leaf_triad_expr_in_nullis_no_zeroCS.csv", leaf_triad_expr_in_nullis_no_zeroCS, row.names = F)


## root ####
head(root_CS_triad_expr_over0.5)
head(root_N1A)
head(root_N1B)
head(root_N1D)

root_triad_expr_in_nullis <- merge(root_CS_triad_expr_over0.5,root_N1A, by.x = "A", by.y="gene")
head(root_triad_expr_in_nullis)
root_triad_expr_in_nullis <- merge(root_triad_expr_in_nullis,root_N1B, by.x = "B", by.y="gene", suffixes = c(".nulliA",".nulliB"))
head(root_triad_expr_in_nullis)
root_triad_expr_in_nullis <- merge(root_triad_expr_in_nullis,root_N1D, by.x = "D", by.y="gene")
head(root_triad_expr_in_nullis)
dim(root_triad_expr_in_nullis)
colnames(root_triad_expr_in_nullis)[10] <- "mean.nulliD"
root_triad_expr_in_nullis <- root_triad_expr_in_nullis[,c(3,2,1,4,5,6,8,9,10)]
head(root_triad_expr_in_nullis)
colnames(root_triad_expr_in_nullis) <- c("A","B","D","CS_tpmA","CS_tpmB","CS_tpmD","tpm_nulliA","tpm_nulliB","tpm_nulliD")
head(root_triad_expr_in_nullis)

# add column with % mismapping in nulliA, nulliB, nulliD and average
root_triad_expr_in_nullis$perc_nulliA_mismap <- root_triad_expr_in_nullis$tpm_nulliA/root_triad_expr_in_nullis$CS_tpmA *100
root_triad_expr_in_nullis$perc_nulliB_mismap <- root_triad_expr_in_nullis$tpm_nulliB/root_triad_expr_in_nullis$CS_tpmB *100
root_triad_expr_in_nullis$perc_nulliD_mismap <- root_triad_expr_in_nullis$tpm_nulliD/root_triad_expr_in_nullis$CS_tpmD *100
root_triad_expr_in_nullis$av_perc_nulli_mismap <- (root_triad_expr_in_nullis$perc_nulliA_mismap +
                                                     root_triad_expr_in_nullis$perc_nulliB_mismap + 
                                                     root_triad_expr_in_nullis$perc_nulliD_mismap)/3
head(root_triad_expr_in_nullis)

# need to remove any triads where A, or B or D is 0 in CS (otherwise I get an infinite fold change in the nulli line)
root_triad_expr_in_nullis_no_zeroCS <- root_triad_expr_in_nullis[root_triad_expr_in_nullis$CS_tpmA != 0 &
                                                                   root_triad_expr_in_nullis$CS_tpmB != 0 &
                                                                   root_triad_expr_in_nullis$CS_tpmD != 0,]
# want to add the triad ID Ricardo uses:
head(hom.df_1_1_1)
head(root_triad_expr_in_nullis_no_zeroCS)
root_triad_expr_in_nullis_no_zeroCS <- merge(root_triad_expr_in_nullis_no_zeroCS, hom.df_1_1_1[,c(1,7)], by.x = "A", by.y="A")
head(root_triad_expr_in_nullis_no_zeroCS)
dim(root_triad_expr_in_nullis_no_zeroCS)

# save the data 
write.csv(file="root_triad_expr_in_nullis_no_zeroCS.csv", root_triad_expr_in_nullis_no_zeroCS, row.names = F)


