
## Carry out analysis of TE family enrichment in gene promoters based on homoeolog expression bias variation.

library(tidyr)

retro_table <- read.csv("Y:/Sophie/transcriptome_paper/retro/retro_promoter_table_ATGstart_withTE.txt", sep="\t", header=FALSE)
colnames(retro_table) <- c("gene.id","distance","TE")
##split the TE by decimal point
retro_table <- separate(data=retro_table, col=TE, into=c("TE","x"), sep="\\.")

##do by movement categories first
setwd("Y:/expression_browser/dataset_reports/20170911_all_reports/02.movement/HC_CS_no_stress/")
top_10_genes <- read.csv("HC_CS_no_stress_movement_top_10pc.txt")
low_10_genes <- read.csv("HC_CS_no_stress_movement_low_10pc.txt")
mid_80_genes <- read.csv("HC_CS_no_stress_movement_middle_80pc.txt")
subset_genes <- rbind(top_10_genes, low_10_genes, mid_80_genes)
#get number of genes that have each TE at 5Kb level
all_genes <- aggregate(gene.id ~ TE, data=retro_table, FUN=length)

#subset into top, mid, and 80
retro_table_top10 <- subset(retro_table, gene.id %in% top_10_genes$x)
retro_table_low10 <- subset(retro_table, gene.id %in% low_10_genes$x)
retro_table_mid80 <- subset(retro_table, gene.id %in% mid_80_genes$x)
retro_table_subset <- subset(retro_table, gene.id %in% subset_genes$x)

all_genes_top10 <- aggregate(gene.id ~ TE, data=retro_table_top10, FUN=length)
all_genes_low10 <- aggregate(gene.id ~ TE, data=retro_table_low10, FUN=length)
all_genes_mid80 <- aggregate(gene.id ~ TE, data=retro_table_mid80, FUN=length)
all_genes_subset <- aggregate(gene.id ~ TE, data=retro_table_subset, FUN=length)

##collapse rows that have a decimal point

##form into table
final_table <- as.data.frame(unique(retro_table$TE))
colnames(final_table) <- c("TE")
final_table <- merge(final_table, all_genes, by="TE")
final_table <- merge(final_table, all_genes_subset, by="TE")
final_table <- merge(final_table, all_genes_top10, by="TE")
final_table <- merge(final_table, all_genes_mid80, by="TE")
final_table <- merge(final_table, all_genes_low10, by="TE")
colnames(final_table) <- c("TE","NumAllGenesWithTE","NumSubsetGenesWithTE","NumTop10WithTE","NumMid80WithTE","NumLow10WithTE")

##add the percentages in now
final_table$PercentAllGenes <- final_table$NumAllGenesWithTE/length(unique(retro_table$gene.id))
final_table$PercentSubsetGenes <- final_table$NumSubsetGenesWithTE/length(unique(retro_table_subset$gene.id))
final_table$PercentTop10 <- final_table$NumTop10WithTE/length(unique(retro_table_top10$gene.id))
final_table$PercentLow10 <- final_table$NumLow10WithTE/length(unique(retro_table_low10$gene.id))
final_table$PercentMid80 <- final_table$NumMid80WithTE/length(unique(retro_table_mid80$gene.id))

final_table$FractionTop10 <- final_table$NumTop10WithTE/final_table$NumSubsetGenesWithTE
final_table$FractionMid80 <- final_table$NumMid80WithTE/final_table$NumSubsetGenesWithTE
final_table$FractionLow10 <- final_table$NumLow10WithTE/final_table$NumSubsetGenesWithTE

final_table$Log2Top10 <- log2(final_table$PercentTop10/final_table$PercentSubsetGenes)
final_table$Log2Low10 <- log2(final_table$PercentLow10/final_table$PercentSubsetGenes)
final_table$Log2Mid80 <- log2(final_table$PercentMid80/final_table$PercentSubsetGenes)

##add chi-square test
expected = c(0.10,0.80,0.10)
chi <- cbind(final_table, t(apply(final_table[4:6], 1, function(x) {
  ch <- chisq.test(x, p=expected)
  c(unname(ch$statistic), ch$p.value)})))
colnames(chi)[19] <-"p.value"
##add adjusted p-value
chi$padj <- p.adjust(chi$p.value, method="BH")

final_table$chisquare <- chi[,20]

setwd("Y://Sophie/transcriptome_paper/Analysis_of_Promoters/Retro/")
write.table(final_table, file="10-80-10_Movement_Categories_By_TE_Family.txt", row.names = FALSE, sep = "\t")