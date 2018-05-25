##
##
## Do Non-Syntenic Comparison for the development subset -- this requires the use of files produced
## by the "kaks_nonsyntenictriads_chinesespring_UPLOAD.R" script. Please run that script first.
##
##


#kaks_nonsyntenic_Dev <- kaks_nonsyntenic[kaks_nonsyntenic$A %in% nonsyn_Dev$x, c("Triad","kaksAB","kaksAD","kaksBD","average")]

##
## Now compare to the top 10, low 10, and middle 80 of the syntenic triads
##

##Here you are reading files made in "kaks_nonsyntenictriads_chinesespring_UPLOAD.R"
setwd("Y:/Sophie/transcriptome_paper/Analysis_of_Promoters/")
kaks_nonsyntenic <- read.csv("kaks_output_nonsyntenic_withaverage.txt",sep="\t")
kaks <- read.csv("kaks_output_withaverage.txt",sep="\t")

##here set directory to the relevant folder containing the correct gene subset files as named below
setwd("Y:/expression_browser/dataset_reports/20170911_all_reports/02.movement/HC_Development/")
top_10_genes <- read.csv("HC_Development_movement_top_10pc.txt")
low_10_genes <- read.csv("HC_Development_movement_low_10pc.txt")
mid_80_genes <- read.csv("HC_Development_movement_middle_80pc.txt")

setwd("Y:/Sophie/transcriptome_paper/non_syntenic_identities/NonSyntenicIdentities/")
nonsyn_Dev <- read.csv("NonSyn_A_homs_Development.txt")

##filter kaks for triads in top 10, low 10, mid 80
kaks_top10 <- kaks[kaks$A %in% top_10_genes$x,c("Triad","kaksAB","kaksAD","kaksBD","average")]
kaks_low10 <- kaks[kaks$A %in% low_10_genes$x,c("Triad","kaksAB","kaksAD","kaksBD","average")]
kaks_mid80 <- kaks[kaks$A %in% mid_80_genes$x,c("Triad","kaksAB","kaksAD","kaksBD","average")]

kaks_nonsyntenic <- kaks_nonsyntenic[kaks_nonsyntenic$A %in% nonsyn_Dev$x, c("Triad","kaksAB","kaksAD","kaksBD","average")]

##get statistics on kaks_nonsyntenic
mean(kaks_nonsyntenic$average, na.rm=TRUE) ##0.3899711
sd(kaks_nonsyntenic$average, na.rm=TRUE) ##0.409161
sd(kaks_nonsyntenic$average, na.rm=TRUE)/sqrt(length(kaks_nonsyntenic$average[!is.na(kaks_nonsyntenic$average)])) ##0.01364628
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

##then adjust p-value for multiple comparisons??
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
