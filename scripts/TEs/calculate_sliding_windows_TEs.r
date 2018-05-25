setwd("W:/Jemima/companion_paper/TEs/")
library(stringr)
library(ggplot2)
library(zoo)
#to speed things up we will only calculate the sliding windows for the genes in 1-1-1 triads 
all_111_triads_stacked <- "W:/Jemima/companion_paper/triads_1-1-1_HC_only_stacked.tsv"
triads_stacked <- read.table(all_111_triads_stacked, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gene_list <- triads_stacked$gene

##the TE info is in files by chromosome, but we can start with just 1A in the first instance

TE_data_path <- "W:/Jemima/companion_paper/TEs/byChr/"
TE_files <- list.files(TE_data_path)
all_TE_data <- data.frame(seqId=character(), geneId=character(), geneStart=numeric(), geneStop=numeric(), geneStrand=numeric(), teId=character(), teCompo=character(), teStatus=character(), tePosition=numeric(), distanceToTe=numeric())

for (file in TE_files){
	data_temp <- read.table(paste(TE_data_path, file, sep = ""), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
	all_TE_data <- rbind(all_TE_data, data_temp)
}

#TE_data <- read.table(TE_data_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

TE_111_NAincl <- subset(all_TE_data, geneId %in% triads_stacked$gene)

#now remove all the genes that don't have any NA values
TE_111 <- TE_111_NAincl[!is.na(TE_111_NAincl$teId),]
#remove inconsistent beginning of string in teId
TE_111$teId <- gsub("seq_", "", TE_111$teId)

TE_111$chunk_start <- str_split_fixed(TE_111$teId, "_", 10)[,6]
#TE_111$chunk_start[TE_111$chunk_start==""] <- "NA"
TE_111$chunk_start <- as.numeric(TE_111$chunk_start)

TE_111$chunk_end <- str_split_fixed(TE_111$teId, "_", 10)[,7]
#TE_111$chunk_end[TE_111$chunk_end==""] <- "NA"
TE_111$chunk_end <- as.numeric(TE_111$chunk_end)

TE_111$TE_size <- TE_111$chunk_end - TE_111$chunk_start

#work out TE start and end and relative positions - plus strand
TE_111_plus <- TE_111[TE_111$geneStrand == "+",]
TE_111_plus$TE_start <- TE_111_plus$geneStart - TE_111_plus$distanceToTe - TE_111_plus$TE_size
TE_111_plus$TE_end <- TE_111_plus$geneStart - TE_111_plus$distanceToTe
TE_111_plus$TE_start_rel <- 0 - TE_111_plus$distanceToTe - TE_111_plus$TE_size
TE_111_plus$TE_end_rel <- 0 - TE_111_plus$distanceToTe

#work out TE start and end and relative positions - minus strand
TE_111_minus <- TE_111[TE_111$geneStrand == "-",]
TE_111_minus$TE_start <- TE_111_minus$geneStart + TE_111_minus$distanceToTe
TE_111_minus$TE_end <- TE_111_minus$geneStart + TE_111_minus$distanceToTe + TE_111_minus$TE_size
TE_111_minus$TE_start_rel <- 0 - TE_111_minus$distanceToTe
TE_111_minus$TE_end_rel <- 0 - TE_111_minus$distanceToTe - TE_111_minus$TE_size

#recombine all genes
TE_111_all <- rbind(TE_111_plus, TE_111_minus)
#set max rel start and end to -5000
TE_111_all[which(TE_111_all$TE_start_rel < -5000), "TE_start_rel"] <- -5000
TE_111_all[which(TE_111_all$TE_end_rel < -5000), "TE_end_rel"] <- -5000

write.table(TE_111_all, file = "W:/Jemima/companion_paper/TEs/all_TE_position_size_relative_positions_ATG_start.tsv", sep = "\t", col.names = TRUE, row.name = FALSE, quote = FALSE)

TE_111_all_simple <- TE_111_all[,c("geneId", "TE_start_rel", "TE_end_rel")]
#set up windows
window_size <- 100
step_size <- 10
window_start <- seq(-5000, (0-window_size), by = step_size)
window_end <- window_start + window_size -1
windows <- data.frame(window_start, window_end)
windows$midpoint <- (windows$window_start + windows$window_end - 1)/2
positions <- data.frame(position=seq(-5000, 0), TE = rep(0, length = 5001))

midpoint <- windows$midpoint
len <- length(midpoint)

#specify the size of the data frame to speed up the process
#must be 3 columns and the number of genes*the number of midpoints
rolling_counts <- data.frame(matrix(NA, nrow = length(gene_list)*len, ncol = 3))
colnames(rolling_counts) <- c("midpoint", "roll_count" , "gene")
rolling_counts$midpoint <- rep(midpoint, length(gene_list))

for(i in seq(1, length(gene_list))){
	gene <- gene_list[i]
	print(paste(gene, ": ",i, " of ", length(gene_list), sep = ""))
	gene_TEs <- TE_111_all_simple[TE_111_all_simple$geneId == gene,]
	positions_gene <- positions
	for (j in seq(1, nrow(gene_TEs))){
		TE_start <- gene_TEs[j,"TE_start_rel"]
		TE_end <- gene_TEs[j,"TE_end_rel"]
		positions_gene[which((positions_gene$position >= TE_start) & (positions_gene$position <= TE_end)), "TE"] <- 1
		positions_gene[which((positions_gene$position <= TE_start) & (positions_gene$position >= TE_end)), "TE"] <- 1
		}
	roll_averages <- rollapply(positions_gene$TE, width = window_size, by = step_size, FUN = sum)
	y <- len*i
	x <- y-len+1
	rolling_counts[c(x:y),"gene"] <- gene
	rolling_counts[c(x:y),"roll_count"] <- roll_averages
	}
	
write.table(rolling_counts, file = "TE_5kb_upstream_rolling_window_100window_10step.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

