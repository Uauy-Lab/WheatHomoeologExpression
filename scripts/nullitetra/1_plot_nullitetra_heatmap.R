# plot nullitetra tpm for WGA
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
metadata_useful
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


# Filter  data to only keep genes expressed >0.5 tpm in WT
tpm_threshold <- 0.5 

genes_over_0.5WT_leaves <- mean_tpmsHC[mean_tpmsHC$Variety =="Chinese Spring" &
                                  mean_tpmsHC$mean > tpm_threshold & mean_tpmsHC$High.level.tissue == "leaves/shoots",]
genes_over_0.5WT_roots <- mean_tpmsHC[mean_tpmsHC$Variety =="Chinese Spring" &
                                         mean_tpmsHC$mean > tpm_threshold & mean_tpmsHC$High.level.tissue == "roots",]


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

#chr5
genes_over_0.5WT_leaves_chr5A <- genes_over_0.5WT_leaves[grepl("CS5A",genes_over_0.5WT_leaves$gene),]

genes_over_0.5WT_leaves_chr5B <- genes_over_0.5WT_leaves[grepl("CS5B",genes_over_0.5WT_leaves$gene),]

genes_over_0.5WT_leaves_chr5D <- genes_over_0.5WT_leaves[grepl("CS5D",genes_over_0.5WT_leaves$gene),]

genes_over_0.5WT_roots_chr5A <- genes_over_0.5WT_roots[grepl("CS5A",genes_over_0.5WT_roots$gene),]

genes_over_0.5WT_roots_chr5B <- genes_over_0.5WT_roots[grepl("CS5B",genes_over_0.5WT_roots$gene),]

genes_over_0.5WT_roots_chr5D <- genes_over_0.5WT_roots[grepl("CS5D",genes_over_0.5WT_roots$gene),]

# make a dataframe for each heatmap to plot using only genes expressed in correct tissue on correct chr in WT

levels(droplevels(mean_tpmsHC$Variety))

#tail(mean_tpmsHC[mean_tpmsHC$gene %in% genes_over_0.5WT_leaves_chr1A$gene & mean_tpmsHC$Variety =="Chinese Spring, N1BT1A",])

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


library("NMF")
library(gplots)
nmf.options(grid.patch=TRUE)
pdf(file="chr1_leaves_test.pdf", width=6, height=6)
aheatmap(matrix_wide_data_chr1A_leaves, scale="row", Rowv = NA, Colv = NA,labRow = NA,color=colorpanel(20,'yellow','red'))
dev.off()

## now do for all chr combinations ###
## for leaves and roots
library(tidyr)
library("NMF")
library(gplots)
nmf.options(grid.patch=TRUE)

chrs <- c("1A","1B","1D","5A","5B","5D")
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

pdf(file=paste0("chr",chr,"_",tissue,".pdf"), width=6, height=6)
aheatmap(matrix_wide_data, scale="row", Rowv = NA, Colv = NA,labRow = NA,color=colorpanel(20,'yellow','red'))
dev.off()

to_add <- data.frame(chromosome= as.character(chr), tissue = as.character(tissue), Chinese_Spring= as.numeric(colMeans(matrix_wide_data))[1],
            NulliATetraB=as.numeric(colMeans(matrix_wide_data))[2], NulliATetraD = as.numeric(colMeans(matrix_wide_data))[3],
           NulliBTetraA=as.numeric(colMeans(matrix_wide_data))[4],NulliBTetraD=as.numeric(colMeans(matrix_wide_data))[5], 
           NulliDTetraA=as.numeric(colMeans(matrix_wide_data))[6], NulliDTetraB=as.numeric(colMeans(matrix_wide_data))[7],
           stringsAsFactors=F)

nullitetra_summary <- rbind(nullitetra_summary, (to_add))
  }
}

nullitetra_summary

write.csv(file="nullitetra_mean_tpm.csv",nullitetra_summary)
