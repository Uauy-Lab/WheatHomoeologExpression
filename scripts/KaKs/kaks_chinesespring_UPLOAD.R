##
##
## Analysis of syn and non-syn mutation rate between homoeologs
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
  kaAB <- ka[1,2]
  kaAD <- ka[1,3]
  kaBD <- ka[2,3]
  ksAB <- ks[1,2]
  ksAD <- ks[1,3]
  ksBD <- ks[2,3]
  ##write line to file
  #triad number
  triad <- strsplit(as.character(files_row[3]),"\\.")[[1]][1]
  # sequence names
  name_A <- align_file$nam[1]
  name_B <- align_file$nam[2]
  name_D <- align_file$nam[3]
  line <- paste(triad,kaksAB,kaksAD,kaksBD,kaAB,kaAD,kaBD,ksAB,ksAD,ksBD,name_A,name_B,name_D,sep="\t")
  #write(line, file="Y:/Sophie/transcriptome_paper/Analysis_of_Promoters/kaks_output.txt",append=TRUE)
  ##also write file with the specific Ks and Ka values
  write(line, file="Y:/Sophie/transcriptome_paper/Analysis_of_Promoters/kaks_output_separatevalues.txt",append=TRUE)
}

##make a list of all folders in file
##Here change the working directory to that containing the un-zipped "pep_cds_alignments" folder 

setwd("Y:/expression_browser/WGA/pep_cds_alignments/")
folders <- list.files("./")
for (i in folders) {
  dir = paste0("Y:/expression_browser/WGA/pep_cds_alignments/",i)
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
## Now filter out top 10, low 10, and mid 80 of movement
##

##here set directory to the relevant folder containing the correct gene subset files as named below
setwd("Y:/expression_browser/dataset_reports/20170911_all_reports/02.movement/HC_CS_no_stress/")
top_10_genes <- read.csv("HC_CS_no_stress_movement_top_10pc.txt")
low_10_genes <- read.csv("HC_CS_no_stress_movement_low_10pc.txt")
mid_80_genes <- read.csv("HC_CS_no_stress_movement_middle_80pc.txt")

setwd("Y:/Sophie/transcriptome_paper/Analysis_of_Promoters/")
kaks <- read.csv("kaks_output.txt",sep="\t")

##replace all "Inf" with 10 (highest number available; indicates that Ks was zero); replace NaN and - numbers with NA
library(schoolmath)
kaks[is.na(kaks)] <- NA
is.na(kaks) <- sapply(kaks, is.infinite)
kaks[is.na(kaks)] <- 10
kaks[kaks<0] <- NA

##obtain "average kaks" for the three pairwise comparisons as well
kaks$average <- rowMeans(kaks[,2:4])
write.table(kaks, "kaks_output_withaverage.txt",sep="\t", row.names = F)
kaks <- read.csv("kaks_output_withaverage.txt",sep="\t")


##filter kaks for triads in top 10, low 10, mid 80
kaks_top10 <- kaks[kaks$A %in% top_10_genes$x,c("Triad","kaksAB","kaksAD","kaksBD","average")]
kaks_low10 <- kaks[kaks$A %in% low_10_genes$x,c("Triad","kaksAB","kaksAD","kaksBD","average")]
kaks_mid80 <- kaks[kaks$A %in% mid_80_genes$x,c("Triad","kaksAB","kaksAD","kaksBD","average")]

##get average and se of kaks for the different sets
mean(kaks_low10$average, na.rm=TRUE) ##0.2108422
sd(kaks_low10$average, na.rm=TRUE) ##0.3539746
sd(kaks_low10$average, na.rm=TRUE)/sqrt(length(kaks_low10$average[!is.na(kaks_low10$average)])) ##0.009386902

mean(kaks_top10$average, na.rm=TRUE) ##0.3301444
sd(kaks_top10$average, na.rm=TRUE) ##0.418373
sd(kaks_top10$average, na.rm=TRUE)/sqrt(length(kaks_top10$average[!is.na(kaks_top10$average)])) ##0.01108686

mean(kaks_mid80$average, na.rm=TRUE) ##0.2577139
sd(kaks_mid80$average, na.rm=TRUE) ##0.3319627
sd(kaks_mid80$average, na.rm=TRUE)/sqrt(length(kaks_mid80$average[!is.na(kaks_mid80$average)])) ##0.003109798

summary(kaks_low10$average)
summary(kaks_top10$average)
summary(kaks_mid80$average)
#quick normality test-- defnitely not normal, so use Mann-Whitney
#shapiro.test(kaks_top10$kaksAB)

##Mann-Whitney
## compare all three kaks measures, and the average kaks
## all have the same default parameters as the first example (shown explicitly)
## note: due to the way in which numbers are "printed" in R, all p-values are capped at 
## <2.2e-16; when the output of the wilcox test is saved and the p-value manually accessed
## you can see the variation in the p-values, as expected. However, for this case, having 
## such small distinctions is immaterial, and sticking with the 2.2e-16 values should be fine.

wilcox.test(kaks_low10$average, kaks_top10$average, alternative="two.sided", paired=FALSE, exact=FALSE)
wilcox.test(kaks_low10$kaksAB, kaks_top10$kaksAB)
wilcox.test(kaks_low10$kaksAD, kaks_top10$kaksAD)
wilcox.test(kaks_low10$kaksBD, kaks_top10$kaksBD)
wilcox.test(kaks_low10$average, kaks_mid80$average)
wilcox.test(kaks_low10$kaksAB, kaks_mid80$kaksAB)
wilcox.test(kaks_low10$kaksAD, kaks_mid80$kaksAD)
wilcox.test(kaks_low10$kaksBD, kaks_mid80$kaksBD)
wilcox.test(kaks_top10$average, kaks_mid80$average)
wilcox.test(kaks_top10$kaksAB, kaks_mid80$kaksAB)
wilcox.test(kaks_top10$kaksAD, kaks_mid80$kaksAD)
wilcox.test(kaks_top10$kaksBD, kaks_mid80$kaksBD)

#prepare for graphs
library(tidyr)
melt_kaks_low10 <- gather(kaks_low10, variable, value, -Triad)
melt_kaks_top10 <- gather(kaks_top10, variable, value, -Triad)
melt_kaks_mid80 <- gather(kaks_mid80, variable, value, -Triad)
melt_kaks_low10$variable <- paste(melt_kaks_low10$variable,"_low10",sep="")
melt_kaks_top10$variable <- paste(melt_kaks_top10$variable,"_top10",sep="")
melt_kaks_mid80$variable <- paste(melt_kaks_mid80$variable,"_mid80",sep="")
melt_kaks<- rbind(melt_kaks_low10,melt_kaks_mid80, melt_kaks_top10)
melt_kaks$variable <- as.factor(melt_kaks$variable)

library(ggplot2)
p1 <- ggplot(melt_kaks[melt_kaks$variable %in% c("average_top10","average_low10","average_mid80"),], aes(x=value, fill=variable))+
  geom_density(alpha=0.4) +
  theme(axis.text.x = element_text(angle=90)) +
  scale_x_continuous(limits = c(0, 2)) +
  ggtitle("Average Ka/Ks")

p2 <- ggplot(melt_kaks[melt_kaks$variable %in% c("kaksAB_top10","kaksAB_low10","kaksAB_mid80"),], aes(x=value, fill=variable))+
  geom_density(alpha=0.4) +
  theme(axis.text.x = element_text(angle=90)) +
  scale_x_continuous(limits = c(0, 2)) +
  ggtitle("A by B pairwise Ka/Ks")

p3 <- ggplot(melt_kaks[melt_kaks$variable %in% c("kaksAD_top10","kaksAD_low10","kaksAD_mid80"),], aes(x=value, fill=variable))+
  geom_density(alpha=0.4) +
  theme(axis.text.x = element_text(angle=90)) +
  scale_x_continuous(limits = c(0, 2)) +
  ggtitle("A by D pairwise Ka/Ks")

p4 <- ggplot(melt_kaks[melt_kaks$variable %in% c("kaksBD_top10","kaksBD_low10","kaksBD_mid80"),], aes(x=value, fill=variable))+
  geom_density(alpha=0.4) +
  theme(axis.text.x = element_text(angle=90)) +
  scale_x_continuous(limits = c(0, 2)) +
  ggtitle("B by D pairwise Ka/Ks")

library(ggpubr)
ggarrange(p1,p2,p3,p4, ncol=2,nrow=2)

ggplot(melt_kaks, aes(x=variable, y=value, fill=variable))+
  geom_boxplot() +
  theme(axis.text.x = element_text(angle=90)) +
  scale_y_continuous(limits = c(0, 2)) +
  #stat_summary(fun.y=mean, geom="errorbar", linetype="dashed", aes(ymax=..y..,ymin=..y..), width=.75) ##add mean line to graphs


##
##
## Also do for 5-90-5 and 1-98-1 subsets of genes

setwd("Y:/expression_browser/dataset_reports/20170911_all_reports/02.movement/HC_CS_no_stress/")
top_5_genes <- read.csv("HC_CS_no_stress_movement_top_5pc.txt")
low_5_genes <- read.csv("HC_CS_no_stress_movement_low_5pc.txt")
mid_90_genes <- read.csv("HC_CS_no_stress_movement_middle_90pc.txt")
top_1_genes <- read.csv("HC_CS_no_stress_movement_top_1pc.txt")
low_1_genes <- read.csv("HC_CS_no_stress_movement_low_1pc.txt")
mid_98_genes <- read.csv("HC_CS_no_stress_movement_middle_98pc.txt")

setwd("Y:/Sophie/transcriptome_paper/Analysis_of_Promoters/")
kaks <- read.csv("kaks_output_withaverage.txt",sep="\t")
avekaks <- summary(kaks$average)
##filter kaks for triads in subsets
kaks_top1 <- kaks[kaks$A %in% top_1_genes$x,c("Triad","kaksAB","kaksAD","kaksBD","average")]
kaks_low1 <- kaks[kaks$A %in% low_1_genes$x,c("Triad","kaksAB","kaksAD","kaksBD","average")]
kaks_mid98 <- kaks[kaks$A %in% mid_98_genes$x,c("Triad","kaksAB","kaksAD","kaksBD","average")]

kaks_top5 <- kaks[kaks$A %in% top_5_genes$x,c("Triad","kaksAB","kaksAD","kaksBD","average")]
kaks_low5 <- kaks[kaks$A %in% low_5_genes$x,c("Triad","kaksAB","kaksAD","kaksBD","average")]
kaks_mid90 <- kaks[kaks$A %in% mid_90_genes$x,c("Triad","kaksAB","kaksAD","kaksBD","average")]

##get average and se of kaks for the different sets
mean(kaks_low5$average, na.rm=TRUE) ##0.2254781
sd(kaks_low5$average, na.rm=TRUE) ##0.4484071
sd(kaks_low5$average, na.rm=TRUE)/sqrt(length(kaks_low5$average[!is.na(kaks_low5$average)])) ##0.01639539

mean(kaks_top5$average, na.rm=TRUE) ##0.32293
sd(kaks_top5$average, na.rm=TRUE) ##0.29212
sd(kaks_top5$average, na.rm=TRUE)/sqrt(length(kaks_top5$average[!is.na(kaks_top5$average)])) ##0.01065961

mean(kaks_mid90$average, na.rm=TRUE) ##0.2632166
sd(kaks_mid90$average, na.rm=TRUE) ##0.340061
sd(kaks_mid90$average, na.rm=TRUE)/sqrt(length(kaks_mid90$average[!is.na(kaks_mid90$average)])) ##0.002925695

summary(kaks_low5$average)
summary(kaks_top5$average)
summary(kaks_mid90$average)

##wilcox-test
wilcox.test(kaks_low1$average, kaks_top1$average, alternative="two.sided", paired=FALSE, exact=FALSE)
wilcox.test(kaks_low1$kaksAB, kaks_top1$kaksAB)
wilcox.test(kaks_low1$kaksAD, kaks_top1$kaksAD)
wilcox.test(kaks_low1$kaksBD, kaks_top1$kaksBD)
wilcox.test(kaks_low1$average, kaks_mid98$average)
wilcox.test(kaks_low1$kaksAB, kaks_mid98$kaksAB)
wilcox.test(kaks_low1$kaksAD, kaks_mid98$kaksAD)
wilcox.test(kaks_low1$kaksBD, kaks_mid98$kaksBD)
wilcox.test(kaks_top1$average, kaks_mid98$average)
wilcox.test(kaks_top1$kaksAB, kaks_mid98$kaksAB)
wilcox.test(kaks_top1$kaksAD, kaks_mid98$kaksAD)
wilcox.test(kaks_top1$kaksBD, kaks_mid98$kaksBD)

wilcox.test(kaks_low5$average, kaks_top5$average, alternative="two.sided", paired=FALSE, exact=FALSE)
wilcox.test(kaks_low5$kaksAB, kaks_top5$kaksAB)
wilcox.test(kaks_low5$kaksAD, kaks_top5$kaksAD)
wilcox.test(kaks_low5$kaksBD, kaks_top5$kaksBD)
wilcox.test(kaks_low5$average, kaks_mid90$average)
wilcox.test(kaks_low5$kaksAB, kaks_mid90$kaksAB)
wilcox.test(kaks_low5$kaksAD, kaks_mid90$kaksAD)
wilcox.test(kaks_low5$kaksBD, kaks_mid90$kaksBD)
wilcox.test(kaks_top5$average, kaks_mid90$average)
wilcox.test(kaks_top5$kaksAB, kaks_mid90$kaksAB)
wilcox.test(kaks_top5$kaksAD, kaks_mid90$kaksAD)
wilcox.test(kaks_top5$kaksBD, kaks_mid90$kaksBD)

##prepare for plotting
library(tidyr)
melt_kaks_low1 <- gather(kaks_low1, variable, value, -Triad)
melt_kaks_top1 <- gather(kaks_top1, variable, value, -Triad)
melt_kaks_mid98 <- gather(kaks_mid98, variable, value, -Triad)
melt_kaks_low1$variable <- paste(melt_kaks_low1$variable,"_low1",sep="")
melt_kaks_top1$variable <- paste(melt_kaks_top1$variable,"_top1",sep="")
melt_kaks_mid98$variable <- paste(melt_kaks_mid98$variable,"_mid98",sep="")
melt_kaks_1<- rbind(melt_kaks_low1,melt_kaks_mid98, melt_kaks_top1)
melt_kaks_1$variable <- as.factor(melt_kaks_1$variable)

melt_kaks_low5 <- gather(kaks_low5, variable, value, -Triad)
melt_kaks_top5 <- gather(kaks_top5, variable, value, -Triad)
melt_kaks_mid90 <- gather(kaks_mid90, variable, value, -Triad)
melt_kaks_low5$variable <- paste(melt_kaks_low5$variable,"_low5",sep="")
melt_kaks_top5$variable <- paste(melt_kaks_top5$variable,"_top5",sep="")
melt_kaks_mid90$variable <- paste(melt_kaks_mid90$variable,"_mid90",sep="")
melt_kaks_5<- rbind(melt_kaks_low5,melt_kaks_mid90, melt_kaks_top5)
melt_kaks_5$variable <- as.factor(melt_kaks_5$variable)

##plots for 5-90-5
library(ggplot2)
p1 <- ggplot(melt_kaks_5[melt_kaks_5$variable %in% c("average_top5","average_low5","average_mid90"),], aes(x=value, fill=variable))+
  geom_density(alpha=0.4) +
  theme(axis.text.x = element_text(angle=90)) +
  scale_x_continuous(limits = c(0, 2)) +
  ggtitle("Average Ka/Ks")

p2 <- ggplot(melt_kaks_5[melt_kaks_5$variable %in% c("kaksAB_top5","kaksAB_low5","kaksAB_mid90"),], aes(x=value, fill=variable))+
  geom_density(alpha=0.4) +
  theme(axis.text.x = element_text(angle=90)) +
  scale_x_continuous(limits = c(0, 2)) +
  ggtitle("A by B pairwise Ka/Ks")

p3 <- ggplot(melt_kaks_5[melt_kaks_5$variable %in% c("kaksAD_top5","kaksAD_low5","kaksAD_mid90"),], aes(x=value, fill=variable))+
  geom_density(alpha=0.4) +
  theme(axis.text.x = element_text(angle=90)) +
  scale_x_continuous(limits = c(0, 2)) +
  ggtitle("A by D pairwise Ka/Ks")

p4 <- ggplot(melt_kaks_5[melt_kaks_5$variable %in% c("kaksBD_top5","kaksBD_low5","kaksBD_mid90"),], aes(x=value, fill=variable))+
  geom_density(alpha=0.4) +
  theme(axis.text.x = element_text(angle=90)) +
  scale_x_continuous(limits = c(0, 2)) +
  ggtitle("B by D pairwise Ka/Ks")

library(ggpubr)
ggarrange(p1,p2,p3,p4, ncol=2,nrow=2)

ggplot(melt_kaks_5, aes(x=variable, y=value, fill=variable))+
  geom_boxplot() +
  theme(axis.text.x = element_text(angle=90)) +
  scale_y_continuous(limits = c(0, 2))

##plots for 1-98-1
p1 <- ggplot(melt_kaks_1[melt_kaks_1$variable %in% c("average_top1","average_low1","average_mid98"),], aes(x=value, fill=variable))+
  geom_density(alpha=0.4) +
  theme(axis.text.x = element_text(angle=90)) +
  scale_x_continuous(limits = c(0, 2)) +
  ggtitle("Average Ka/Ks")

p2 <- ggplot(melt_kaks_1[melt_kaks_1$variable %in% c("kaksAB_top1","kaksAB_low1","kaksAB_mid98"),], aes(x=value, fill=variable))+
  geom_density(alpha=0.4) +
  theme(axis.text.x = element_text(angle=90)) +
  scale_x_continuous(limits = c(0, 2)) +
  ggtitle("A by B pairwise Ka/Ks")

p3 <- ggplot(melt_kaks_1[melt_kaks_1$variable %in% c("kaksAD_top1","kaksAD_low1","kaksAD_mid98"),], aes(x=value, fill=variable))+
  geom_density(alpha=0.4) +
  theme(axis.text.x = element_text(angle=90)) +
  scale_x_continuous(limits = c(0, 2)) +
  ggtitle("A by D pairwise Ka/Ks")

p4 <- ggplot(melt_kaks_1[melt_kaks_1$variable %in% c("kaksBD_top1","kaksBD_low1","kaksBD_mid98"),], aes(x=value, fill=variable))+
  geom_density(alpha=0.4) +
  theme(axis.text.x = element_text(angle=90)) +
  scale_x_continuous(limits = c(0, 2)) +
  ggtitle("B by D pairwise Ka/Ks")

library(ggpubr)
ggarrange(p1,p2,p3,p4, ncol=2,nrow=2)

ggplot(melt_kaks_1, aes(x=variable, y=value, fill=variable))+
  geom_boxplot() +
  theme(axis.text.x = element_text(angle=90)) +
  scale_y_continuous(limits = c(0, 2))

##
##
## Also do for 25-50-25 subsets of genes

##
## Now filter out top 10, low 10, and mid 80 of movement
##

setwd("Y:/expression_browser/dataset_reports/20170911_all_reports/02.movement/HC_CS_no_stress/")
top_25_genes <- read.csv("HC_CS_no_stress_movement_top_25pc.txt")
low_25_genes <- read.csv("HC_CS_no_stress_movement_low_25pc.txt")
mid_50_genes <- read.csv("HC_CS_no_stress_movement_middle_50pc.txt")

setwd("Y:/Sophie/transcriptome_paper/Analysis_of_Promoters/")
kaks <- read.csv("kaks_output.txt",sep="\t")

##replace all "Inf" with 10 (highest number available; indicates that Ks was zero); replace NaN and - numbers with NA
library(schoolmath)
kaks[is.na(kaks)] <- NA
is.na(kaks) <- sapply(kaks, is.infinite)
kaks[is.na(kaks)] <- 10
kaks[kaks<0] <- NA

##obtain "average kaks" for the three pairwise comparisons as well
kaks$average <- rowMeans(kaks[,2:4])
write.table(kaks, "kaks_output_withaverage.txt",sep="\t", row.names = F)
kaks <- read.csv("kaks_output_withaverage.txt",sep="\t")


##filter kaks for triads in top 25, low 25, mid 50
kaks_top25 <- kaks[kaks$A %in% top_25_genes$x,c("Triad","kaksAB","kaksAD","kaksBD","average")]
kaks_low25 <- kaks[kaks$A %in% low_25_genes$x,c("Triad","kaksAB","kaksAD","kaksBD","average")]
kaks_mid50 <- kaks[kaks$A %in% mid_50_genes$x,c("Triad","kaksAB","kaksAD","kaksBD","average")]

##get average and se of kaks for the different sets
mean(kaks_low25$average, na.rm=TRUE) ##0.2208852
sd(kaks_low25$average, na.rm=TRUE) #0.330918
sd(kaks_low25$average, na.rm=TRUE)/sqrt(length(kaks_low25$average[!is.na(kaks_low25$average)])) ##0.005545419

mean(kaks_top25$average, na.rm=TRUE) ##0.2986229
sd(kaks_top25$average, na.rm=TRUE) ##0.341176
sd(kaks_top25$average, na.rm=TRUE)/sqrt(length(kaks_top25$average[!is.na(kaks_top25$average)])) ##0.005721338

mean(kaks_mid50$average, na.rm=TRUE) ##0.2608251
sd(kaks_mid50$average, na.rm=TRUE) ##0.3512537
sd(kaks_mid50$average, na.rm=TRUE)/sqrt(length(kaks_mid50$average[!is.na(kaks_mid50$average)])) ##0.004161587

summary(kaks_low25$average)
summary(kaks_top25$average)
summary(kaks_mid50$average)
#quick normality test-- defnitely not normal, so use Mann-Whitney
#shapiro.test(kaks_top10$kaksAB)

##Mann-Whitney
## compare all three kaks measures, and the average kaks
## all have the same default parameters as the first example (shown explicitly)
## note: due to the way in which numbers are "printed" in R, all p-values are capped at 
## <2.2e-16; when the output of the wilcox test is saved and the p-value manually accessed
## you can see the variation in the p-values, as expected. However, for this case, having 
## such small distinctions is immaterial, and sticking with the 2.2e-16 values should be fine.

wilcox.test(kaks_low25$average, kaks_top25$average, alternative="two.sided", paired=FALSE, exact=FALSE)
wilcox.test(kaks_low25$kaksAB, kaks_top25$kaksAB)
wilcox.test(kaks_low25$kaksAD, kaks_top25$kaksAD)
wilcox.test(kaks_low25$kaksBD, kaks_top25$kaksBD)
wilcox.test(kaks_low25$average, kaks_mid50$average)
wilcox.test(kaks_low25$kaksAB, kaks_mid50$kaksAB)
wilcox.test(kaks_low25$kaksAD, kaks_mid50$kaksAD)
wilcox.test(kaks_low25$kaksBD, kaks_mid50$kaksBD)
wilcox.test(kaks_top25$average, kaks_mid50$average)
wilcox.test(kaks_top25$kaksAB, kaks_mid50$kaksAB)
wilcox.test(kaks_top25$kaksAD, kaks_mid50$kaksAD)
wilcox.test(kaks_top25$kaksBD, kaks_mid50$kaksBD)

#prepare for graphs
library(tidyr)
melt_kaks_low25 <- gather(kaks_low25, variable, value, -Triad)
melt_kaks_top25 <- gather(kaks_top25, variable, value, -Triad)
melt_kaks_mid50 <- gather(kaks_mid50, variable, value, -Triad)
melt_kaks_low25$variable <- paste(melt_kaks_low25$variable,"_low25",sep="")
melt_kaks_top25$variable <- paste(melt_kaks_top25$variable,"_top25",sep="")
melt_kaks_mid50$variable <- paste(melt_kaks_mid50$variable,"_mid50",sep="")
melt_kaks<- rbind(melt_kaks_low25,melt_kaks_mid50, melt_kaks_top25)
melt_kaks$variable <- as.factor(melt_kaks$variable)


##plots for 25-50-25
p1 <- ggplot(melt_kaks[melt_kaks$variable %in% c("average_top25","average_low25","average_mid50"),], aes(x=value, fill=variable))+
  geom_density(alpha=0.4) +
  theme(axis.text.x = element_text(angle=90)) +
  scale_x_continuous(limits = c(0, 2)) +
  ggtitle("Average Ka/Ks")

p2 <- ggplot(melt_kaks[melt_kaks$variable %in% c("kaksAB_top25","kaksAB_low25","kaksAB_mid50"),], aes(x=value, fill=variable))+
  geom_density(alpha=0.4) +
  theme(axis.text.x = element_text(angle=90)) +
  scale_x_continuous(limits = c(0, 2)) +
  ggtitle("A by B pairwise Ka/Ks")

p3 <- ggplot(melt_kaks[melt_kaks$variable %in% c("kaksAD_top25","kaksAD_low25","kaksAD_mid50"),], aes(x=value, fill=variable))+
  geom_density(alpha=0.4) +
  theme(axis.text.x = element_text(angle=90)) +
  scale_x_continuous(limits = c(0, 2)) +
  ggtitle("A by D pairwise Ka/Ks")

p4 <- ggplot(melt_kaks[melt_kaks$variable %in% c("kaksBD_top25","kaksBD_low25","kaksBD_mid50"),], aes(x=value, fill=variable))+
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