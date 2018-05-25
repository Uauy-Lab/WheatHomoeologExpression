##
##
## Analysis of synonymous and non-synonymous mutation rate between homoeologs
##
##



##usage kaks(x, verbose = FALSE, debug = FALSE, forceUpperCase = TRUE, rmgap = TRUE)
##
## input alignment file as x; need to be aligned at amino-acid level, then converted back to the relevant codons. (consider reverse.align)
## would be good to align all three homoeologs at one, and then kaks can parse them all for each triad


library(seqinr)

## function to get rev_align file
## Note: part of this function is to remove any codons that are "gaps" in the alignment; this is required for the Ka/Ks package

rev_align <- function(files_row){
  file_out <- paste0(files_row[2],"revalign.fa")
  reverse.align(nucl.file = files_row[2], protaln.file = files_row[1],
                input.format = 'fasta', out.file = file_out)
}

#function to get Ka/Ks values for the homoeolog comparisons
kaks_func <- function(files_row){
  align_file <- read.alignment(file=files_row[3], format="fasta")
  kaks_output <- kaks(align_file)
  ##extract ka and ks values for each sample and divide to get ka/ks
  ka <- as.matrix(kaks_output[["ka"]])
  ks <- as.matrix(kaks_output[["ks"]])
  ##get the ka/ks value for each of the three comparisons
  kaksAB <- ka[1,2]/ks[1,2]
  kaksAD <- ka[1,3]/ks[1,3]
  kaksBD <- ka[2,3]/ks[2,3]
  ##write line to file
  #triad number
  triad <- strsplit(as.character(files_row[3]),"\\.")[[1]][1]
  # sequence names
  name_A <- align_file$nam[1]
  name_B <- align_file$nam[2]
  name_D <- align_file$nam[3]
  line <- paste(triad,kaksAB,kaksAD,kaksBD,name_A,name_B,name_D,sep="\t")
  write(line, file="Y:/Sophie/transcriptome_paper/Analysis_of_Promoters/kaks_output_non_syntenic.txt",append=TRUE)
}

##make a list of all folders in file
##Here change the working directory to that containing the un-zipped "NonSyntenicIdentities" folder 

setwd("Y:/Sophie/transcriptome_paper/non_syntenic_identities/NonSyntenicIdentities/alignments/")
folders <- list.files("./")
for (i in folders) {
  dir = paste0("Y:/Sophie/transcriptome_paper/non_syntenic_identities/NonSyntenicIdentities/alignments/",i)
  setwd(dir)
  pep_files <- as.data.frame(list.files("./",pattern="*.pep.fa"))
  cds_files <- list.files("./",pattern="*.cds.fa")
  files_table <- cbind.data.frame(pep_files,cds_files)
  colnames(files_table) <- c("pep","cds")
  ##function to get the reverse aligned file for each row in data frame
  apply(files_table, 1, rev_align)
  ##now get list of all rev_aligned files
  files_table <- cbind.data.frame(files_table, list.files("./",pattern="*revalign.fa"))
  colnames(files_table) <- c("pep","cds","revalign")
  ##now carry out the kaks function
  apply(files_table, 1, kaks_func)
}

##
## Now compare to the top 10, low 10, and middle 80 of the triads
##

setwd("Y:/Sophie/transcriptome_paper/Analysis_of_Promoters/")
kaks <- read.csv("kaks_output_non_syntenic.txt",sep="\t")

##replace all "Inf" with 10 (highest number available; indicates that Ks was zero); replace NaN and - numbers with NA
library(schoolmath)
kaks[is.na(kaks)] <- NA
is.na(kaks) <- sapply(kaks, is.infinite)
kaks[is.na(kaks)] <- 10
kaks[kaks<0] <- NA

##obtain "average kaks" for the three pairwise comparisons as well
kaks$average <- rowMeans(kaks[,2:4])
write.table(kaks, "kaks_output_nonsyntenic_withaverage.txt",sep="\t", row.names = F)
kaks_nonsyntenic <- read.csv("kaks_output_nonsyntenic_withaverage.txt",sep="\t")

kaks <- read.csv("kaks_output_withaverage.txt",sep="\t")

##here set directory to the relevant folder containing the correct gene subset files as named below
setwd("Y:/expression_browser/dataset_reports/20170911_all_reports/02.movement/HC_CS_no_stress/")
top_10_genes <- read.csv("HC_CS_no_stress_movement_top_10pc.txt")
low_10_genes <- read.csv("HC_CS_no_stress_movement_low_10pc.txt")
mid_80_genes <- read.csv("HC_CS_no_stress_movement_middle_80pc.txt")

setwd("Y:/Sophie/transcriptome_paper/non_syntenic_identities/NonSyntenicIdentities/")
nonsyn_genes <- read.csv("NonSyn_A_homs_HCCS_nostress.txt")

##filter kaks for triads in top 10, low 10, mid 80
kaks_top10 <- kaks[kaks$A %in% top_10_genes$x,c("Triad","kaksAB","kaksAD","kaksBD","average")]
kaks_low10 <- kaks[kaks$A %in% low_10_genes$x,c("Triad","kaksAB","kaksAD","kaksBD","average")]
kaks_mid80 <- kaks[kaks$A %in% mid_80_genes$x,c("Triad","kaksAB","kaksAD","kaksBD","average")]

kaks_nonsyntenic <- kaks_nonsyntenic[kaks_nonsyntenic$A %in% nonsyn_genes$x, c("Triad","kaksAB","kaksAD","kaksBD","average")]

##get statistics on kaks_nonsyntenic
mean(kaks_nonsyntenic$average, na.rm=TRUE) ##0.3885003
sd(kaks_nonsyntenic$average, na.rm=TRUE) ##0.4187972
sd(kaks_nonsyntenic$average, na.rm=TRUE)/sqrt(length(kaks_nonsyntenic$average[!is.na(kaks_nonsyntenic$average)])) ##0.01456301
summary(kaks_nonsyntenic$average)

##Mann-Whitney
## compare all three kaks measures, and the average kaks
## all have the same default parameters as the first example (shown explicitly)
## note: due to the way in which numbers are "printed" in R, all p-values are capped at 
## <2.2e-16; when the output of the wilcox test is saved and the p-value manually accessed
## you can see the variation in the p-values, as expected. However, for this case, having 
## such small distinctions is immaterial, and sticking with the 2.2e-16 values should be fine.

wilcox.test(kaks_low10$average, kaks_nonsyntenic$average, alternative="two.sided", paired=FALSE, exact=FALSE)
wilcox.test(kaks_low10$kaksAB, kaks_nonsyntenic$kaksAB)
wilcox.test(kaks_low10$kaksAD, kaks_nonsyntenic$kaksAD)
wilcox.test(kaks_low10$kaksBD, kaks_nonsyntenic$kaksBD)
wilcox.test(kaks_mid80$average, kaks_nonsyntenic$average)
wilcox.test(kaks_mid80$kaksAB, kaks_nonsyntenic$kaksAB)
wilcox.test(kaks_mid80$kaksAD, kaks_nonsyntenic$kaksAD)
wilcox.test(kaks_mid80$kaksBD, kaks_nonsyntenic$kaksBD)
wilcox.test(kaks_top10$average, kaks_nonsyntenic$average)
wilcox.test(kaks_top10$kaksAB, kaks_nonsyntenic$kaksAB)
wilcox.test(kaks_top10$kaksAD, kaks_nonsyntenic$kaksAD)
wilcox.test(kaks_top10$kaksBD, kaks_nonsyntenic$kaksBD)

##then adjust p-value for multiple comparisons
library(tidyr)
melt_kaks_low10 <- gather(kaks_low10, variable, value, -Triad)
melt_kaks_top10 <- gather(kaks_top10, variable, value, -Triad)
melt_kaks_mid80 <- gather(kaks_mid80, variable, value, -Triad)
melt_kaks_nonsyntenic <- gather(kaks_nonsyntenic, variable, value, -Triad)

melt_kaks_low10$variable <- paste(melt_kaks_low10$variable,"_low10",sep="")
melt_kaks_top10$variable <- paste(melt_kaks_top10$variable,"_top10",sep="")
melt_kaks_mid80$variable <- paste(melt_kaks_mid80$variable,"_mid80",sep="")
melt_kaks_nonsyntenic$variable <- paste(melt_kaks_nonsyntenic$variable,"_NonSyn",sep="")


melt_kaks<- rbind(melt_kaks_low10,melt_kaks_mid80, melt_kaks_top10, melt_kaks_nonsyntenic)
melt_kaks$variable <- as.factor(melt_kaks$variable)

library(ggplot2)
p1 <- ggplot(melt_kaks[melt_kaks$variable %in% c("average_top10","average_low10","average_mid80", "average_NonSyn"),], aes(x=value, fill=variable))+
  geom_density(alpha=0.4) +
  theme(axis.text.x = element_text(angle=90)) +
  scale_x_continuous(limits = c(0, 2)) +
  ggtitle("Average Ka/Ks")

p2 <- ggplot(melt_kaks[melt_kaks$variable %in% c("kaksAB_top10","kaksAB_low10","kaksAB_mid80", "kaksAB_NonSyn"),], aes(x=value, fill=variable))+
  geom_density(alpha=0.4) +
  theme(axis.text.x = element_text(angle=90)) +
  scale_x_continuous(limits = c(0, 2)) +
  ggtitle("A by B pairwise Ka/Ks")

p3 <- ggplot(melt_kaks[melt_kaks$variable %in% c("kaksAD_top10","kaksAD_low10","kaksAD_mid80", "kaksAD_NonSyn"),], aes(x=value, fill=variable))+
  geom_density(alpha=0.4) +
  theme(axis.text.x = element_text(angle=90)) +
  scale_x_continuous(limits = c(0, 2)) +
  ggtitle("A by D pairwise Ka/Ks")

p4 <- ggplot(melt_kaks[melt_kaks$variable %in% c("kaksBD_top10","kaksBD_low10","kaksBD_mid80", "kaksBD_NonSyn"),], aes(x=value, fill=variable))+
  geom_density(alpha=0.4) +
  theme(axis.text.x = element_text(angle=90)) +
  scale_x_continuous(limits = c(0, 2)) +
  ggtitle("B by D pairwise Ka/Ks")

library(ggpubr)
ggarrange(p1,p2,p3,p4, ncol=2,nrow=2)

ggplot(melt_kaks, aes(x=variable, y=value, fill=variable))+
  geom_boxplot() +
  theme(axis.text.x = element_text(angle=90)) +
  scale_y_continuous(limits = c(0, 2))
