# Want to summarise number of DE genes per tissue
# Philippa Borrill
# 24.7.2017

setwd("Y:\\expression_browser\\WGA\\Bayer_analysis\\DE\\")

tpms <- read.csv(file="Y:\\expression_browser\\WGA\\WGCNA\\expressionValuesPerGene\\Development_tpm.tsv", header=T, sep = "\t")
dim(tpms)

# read in metadata
metadata <- read.csv(file="Y:\\expression_browser\\WGA\\WGCNA\\HC_only\\WGCNA_metadata_cluster_with_intermed_coarse.csv", header=T)
head(metadata)
metadata[1:4,1:4]
colnames(metadata)
metadata <- metadata[metadata$Provider == "M. Davey BCS",] 
dim(metadata)

# now need to make metadata into the same order as the counts dataframe
df_counts_order <- as.data.frame(colnames(tpms))
head(df_counts_order)

ordered_metadata <- merge(df_counts_order, metadata[,c(4,22)], by.x="colnames(tpms)", by.y="Sample.IDs")
head(ordered_metadata)
colnames(ordered_metadata) <- c("Sample.IDs","Intermediate")


# check metadata and columns are in the same order
ordered_metadata$Sample.IDs == colnames(tpms)


#samples <- data.frame(row.names=c("RNAi_rep1","RNAi_rep2","RNAi_rep3","RNAi_rep4","WT_rep1","WT_rep2", "WT_rep3"), condition=as.factor(c(rep("RNAi",4),rep("WT",3))))
samples <- data.frame(row.names=colnames(tpms), condition=as.factor(as.character(ordered_metadata$Intermediate)))

samples

#i <- "anther"
#j <- "awns"

fold_changes <- data.frame(tissue1=character(),vs_tissue2=character(), 
                           FC_over_2UP_padj_0.05=numeric(), FC_over_2DOWN_padj_0.05=numeric(),
                           FC_over_2UP_padj_0.01=numeric(), FC_over_2DOWN_padj_0.01=numeric(),
                           FC_over_2UP_padj_0.001=numeric(), FC_over_2DOWN_padj_0.001=numeric(),
                           FC_over_1.5UP_padj_0.05=numeric(), FC_over_1.5DOWN_padj_0.05=numeric(),
                           FC_over_1.5UP_padj_0.01=numeric(), FC_over_1.5DOWN_padj_0.01=numeric(),
                           FC_over_1.5UP_padj_0.001=numeric(), FC_over_1.5DOWN_padj_0.001=numeric(),
                           stringsAsFactors=F)


# for all conditions (22 conditions)
for (j in unique(samples$condition)) {
  
  # for all conditions except the already selected condition
  for (i in unique(samples$condition)[unique(samples$condition)!= j]) { 
    
    bck_res <- read.csv(file=paste(i,"_vs_",j,"_results_padj_0.05.csv"))
    
    # sort results on padj
    ordered_res <- bck_res[order(bck_res$padj),]
    head(ordered_res)
    
    ordered_res_na.rm <- na.omit(ordered_res)
    
    ordered_res_na.rm_padj0.05 <- ordered_res_na.rm[ordered_res_na.rm$padj<0.05,] 
    
    # calculate number of genes over certain threshold
    FC_over_2UP_padj_0.05 <- nrow(ordered_res_na.rm_padj0.05[(ordered_res_na.rm_padj0.05$padj<0.05) & (ordered_res_na.rm_padj0.05$log2FoldChange > log2(2)),])
    FC_over_2DOWN_padj_0.05 <- nrow(ordered_res_na.rm_padj0.05[(ordered_res_na.rm_padj0.05$padj<0.05) & (ordered_res_na.rm_padj0.05$log2FoldChange < log2(0.5)),])
    FC_over_2UP_padj_0.01 <- nrow(ordered_res_na.rm_padj0.05[(ordered_res_na.rm_padj0.05$padj<0.01) & (ordered_res_na.rm_padj0.05$log2FoldChange > log2(2)),])
    FC_over_2DOWN_padj_0.01 <- nrow(ordered_res_na.rm_padj0.05[(ordered_res_na.rm_padj0.05$padj<0.01) & (ordered_res_na.rm_padj0.05$log2FoldChange < log2(0.5)),])
    FC_over_2UP_padj_0.001 <- nrow(ordered_res_na.rm_padj0.05[(ordered_res_na.rm_padj0.05$padj<0.001) & (ordered_res_na.rm_padj0.05$log2FoldChange > log2(2)),])
    FC_over_2DOWN_padj_0.001 <- nrow(ordered_res_na.rm_padj0.05[(ordered_res_na.rm_padj0.05$padj<0.001) & (ordered_res_na.rm_padj0.05$log2FoldChange < log2(0.5)),])
    FC_over_1.5UP_padj_0.05 <- nrow(ordered_res_na.rm_padj0.05[(ordered_res_na.rm_padj0.05$padj<0.05) & (ordered_res_na.rm_padj0.05$log2FoldChange > log2(1.5)),])
    FC_over_1.5DOWN_padj_0.05 <- nrow(ordered_res_na.rm_padj0.05[(ordered_res_na.rm_padj0.05$padj<0.05) & (ordered_res_na.rm_padj0.05$log2FoldChange < log2(2/3)),])
    FC_over_1.5UP_padj_0.01 <- nrow(ordered_res_na.rm_padj0.05[(ordered_res_na.rm_padj0.05$padj<0.01) & (ordered_res_na.rm_padj0.05$log2FoldChange > log2(1.5)),])
    FC_over_1.5DOWN_padj_0.01 <- nrow(ordered_res_na.rm_padj0.05[(ordered_res_na.rm_padj0.05$padj<0.01) & (ordered_res_na.rm_padj0.05$log2FoldChange < log2(2/3)),])
    FC_over_1.5UP_padj_0.001 <- nrow(ordered_res_na.rm_padj0.05[(ordered_res_na.rm_padj0.05$padj<0.001) & (ordered_res_na.rm_padj0.05$log2FoldChange > log2(1.5)),])
    FC_over_1.5DOWN_padj_0.001 <- nrow(ordered_res_na.rm_padj0.05[(ordered_res_na.rm_padj0.05$padj<0.001) & (ordered_res_na.rm_padj0.05$log2FoldChange < log2(2/3)),])
    
    identifier <- paste0(i,"_vs_",j)
    
    fold_changes[nrow(fold_changes)+1,] <- c(as.character(i),as.character(j),
                                             FC_over_2UP_padj_0.05,FC_over_2DOWN_padj_0.05,
                                             FC_over_2UP_padj_0.01,FC_over_2DOWN_padj_0.01,
                                             FC_over_2UP_padj_0.001,FC_over_2DOWN_padj_0.001,
                                             FC_over_1.5UP_padj_0.05,FC_over_1.5DOWN_padj_0.05,
                                             FC_over_1.5UP_padj_0.01,FC_over_1.5DOWN_padj_0.01,
                                             FC_over_1.5UP_padj_0.001,FC_over_1.5DOWN_padj_0.001)
    
  }
}
head(fold_changes)

write.csv(file="number_of_genes_up_down_reg_various_padj.csv", fold_changes)
