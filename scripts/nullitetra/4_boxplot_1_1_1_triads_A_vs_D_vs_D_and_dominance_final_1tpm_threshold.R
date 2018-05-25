# Plot A, B, D expression in nullitetras, do they have a difference in % mismapping?
# Do genes which are dominance, balanced or suppressed have different % of mismapping in nullitetras?
# 28.2.2018

setwd("C:\\Users\\borrillp\\Documents\\28.2.18\\")
library(reshape2)

##### leaf #####
leaf_nullitetra_dom <- read.csv(file="leaf_triad_expr_in_nullis_no_zeroCS_with_dominance.csv")
# select only triads which are >1 tpm in CS
head(leaf_nullitetra_dom)
dim(leaf_nullitetra_dom)

leaf_nullitetra_dom <- leaf_nullitetra_dom[(leaf_nullitetra_dom$CS_tpmA + leaf_nullitetra_dom$CS_tpmB + leaf_nullitetra_dom$CS_tpmD) > 1,]
head(leaf_nullitetra_dom)
dim(leaf_nullitetra_dom)

t.test(leaf_nullitetra_dom$perc_nulliA_mismap, leaf_nullitetra_dom$perc_nulliB_mismap)
t.test(leaf_nullitetra_dom$perc_nulliA_mismap, leaf_nullitetra_dom$perc_nulliD_mismap)

# first make graphs not split by dominance category

leaf_nullitetra_long <- melt(leaf_nullitetra_dom[,c(1,14:16)],
                                 id.vars = c("group_id"))
head(leaf_nullitetra_long)
tail(leaf_nullitetra_long)

kruskal.test(value ~ variable, data=leaf_nullitetra_long)

### Dunn test
#install.packages("FSA")
library(FSA)
PT = dunnTest(value ~ variable,
              data=leaf_nullitetra_long,
              method="by")
PT

# plot graph 
library(ggplot2)
pdf(file="leaf_triad_expr_in_nullis_no_zeroCS_1tpm_threshold.pdf")
# plot boxplot without the outliers, zoom using coord cartesian (doesnt change the data just zooms)
ggplot(data=leaf_nullitetra_long, aes(x= variable,y=value)) + 
  geom_boxplot(outlier.shape = NA)  +
  coord_cartesian(ylim = c(0,12)) + 
  scale_y_continuous(breaks=seq(0, 12, 2)) + theme_bw() +
  ylab("Mismapping (%)") + xlab(NULL) +
  scale_x_discrete(breaks=c("perc_nulliA_mismap", "perc_nulliB_mismap", "perc_nulliD_mismap"),
                   labels=c("nulli A", "nulli B", "nulli D")) 
dev.off()

# calc summary stats
library(tidyr)
library(dplyr)
summarised_leaf <- leaf_nullitetra_long %>%
  group_by(variable) %>%
  summarise(mean = mean(value),
            Q1 = quantile(value,0.25),
            Q2 = quantile(value,0.5),
            Q3= quantile(value,0.75),
            percentile90 = quantile(value,0.9),
            percentile95 = quantile(value,0.95),
            num_triads = length(value))

as.data.frame(summarised_leaf)
write.csv(file="leaf_triad_expr_in_nullis_no_zeroCS_1tpm_threshold_stats.csv",summarised_leaf, row.names = F)


# convert to long format for ggplot2
head(leaf_nullitetra_dom)
leaf_nullitetra_dom_long <- melt(leaf_nullitetra_dom[,c(1,14:19)],
                                 id.vars = c("group_id","description","general_description"))
head(leaf_nullitetra_dom_long)


# now get A dom in nulliA, B dom in nulliB, D dom in nulliD
A.dom <- leaf_nullitetra_dom_long[leaf_nullitetra_dom_long$description == "A.dominant" & leaf_nullitetra_dom_long$variable == "perc_nulliA_mismap",]
head(A.dom)
B.dom <- leaf_nullitetra_dom_long[leaf_nullitetra_dom_long$description == "B.dominant" & leaf_nullitetra_dom_long$variable == "perc_nulliB_mismap",]
head(B.dom)
D.dom <- leaf_nullitetra_dom_long[leaf_nullitetra_dom_long$description == "D.dominant" & leaf_nullitetra_dom_long$variable == "perc_nulliD_mismap",]
head(D.dom)

dom <- rbind(A.dom, B.dom, D.dom)
head(dom)
tail(dom)
dim(dom)

# now get A supp in nulliA, B supp in nulliB, D supp in nulliD
A.supp <- leaf_nullitetra_dom_long[leaf_nullitetra_dom_long$description == "A.suppressed" & leaf_nullitetra_dom_long$variable == "perc_nulliA_mismap",]
head(A.supp)
B.supp <- leaf_nullitetra_dom_long[leaf_nullitetra_dom_long$description == "B.suppressed" & leaf_nullitetra_dom_long$variable == "perc_nulliB_mismap",]
head(B.supp)
D.supp <- leaf_nullitetra_dom_long[leaf_nullitetra_dom_long$description == "D.suppressed" & leaf_nullitetra_dom_long$variable == "perc_nulliD_mismap",]
head(D.supp)

supp <- rbind(A.supp, B.supp, D.supp)
head(supp)
tail(supp)
dim(supp)

# now get Central
head(leaf_nullitetra_dom)
central <- leaf_nullitetra_dom_long[leaf_nullitetra_dom_long$description == "Central" & leaf_nullitetra_dom_long$variable != "av_perc_nulli_mismap",]
head(central)
tail(central)
dim(central)

dom_supp_central <- rbind(dom, supp, central)
head(dom_supp_central)
tail(dom_supp_central)

# plot graph 
library(ggplot2)
pdf(file="leaf_triad_expr_in_nullis_no_zeroCS_dominance_summary.pdf")
ggplot(data=dom_supp_central, aes(x= general_description,y=value)) + 
  geom_boxplot(outlier.shape = NA)  +
  coord_cartesian(ylim = c(0,20)) + 
  scale_y_continuous(breaks=seq(0, 20, 2)) + theme_bw() +
  ylab("Mismapping (%)") + xlab(NULL)  
dev.off()

# calc summary stats
summarised_leaf_dom <- dom_supp_central %>%
  group_by(general_description) %>%
  summarise(mean = mean(value),
            Q1 = quantile(value,0.25),
            Q2 = quantile(value,0.5),
            Q3= quantile(value,0.75),
            percentile90 = quantile(value,0.9),
            percentile95 = quantile(value,0.95),
            num_genes = length(value))

as.data.frame(summarised_leaf_dom)

write.csv(file="leaf_triad_expr_in_nullis_no_zeroCS_dominance_summary.csv",summarised_leaf_dom, row.names = F)

kruskal.test(value ~ general_description, data=dom_supp_central)

### Dunn test
#install.packages("FSA")
library(FSA)
PT = dunnTest(value ~ general_description,
              data=dom_supp_central,
              method="by")
PT

# now want to also include dominant_non_dominant and suppressed_non_suppressed
A.dom.non.dom <- leaf_nullitetra_dom_long[leaf_nullitetra_dom_long$description == "A.dominant" & 
                                            (leaf_nullitetra_dom_long$variable == "perc_nulliB_mismap" | 
                                               leaf_nullitetra_dom_long$variable == "perc_nulliD_mismap"),]
head(A.dom.non.dom)
A.dom.non.dom$general_description <- "Dom_non_dom"
head(A.dom.non.dom)

B.dom.non.dom <- leaf_nullitetra_dom_long[leaf_nullitetra_dom_long$description == "B.dominant" &  
                                            (leaf_nullitetra_dom_long$variable == "perc_nulliA_mismap" | 
                                               leaf_nullitetra_dom_long$variable == "perc_nulliD_mismap"),]
head(B.dom.non.dom)
B.dom.non.dom$general_description <- "Dom_non_dom"
D.dom.non.dom <- leaf_nullitetra_dom_long[leaf_nullitetra_dom_long$description == "D.dominant" &  
                                            (leaf_nullitetra_dom_long$variable == "perc_nulliB_mismap" | 
                                               leaf_nullitetra_dom_long$variable == "perc_nulliA_mismap"),]
head(D.dom.non.dom)
D.dom.non.dom$general_description <- "Dom_non_dom"

# now get A supp in nulliA, B supp in nulliB, D supp in nulliD
A.supp.non.supp <- leaf_nullitetra_dom_long[leaf_nullitetra_dom_long$description == "A.suppressed" & 
                                              (leaf_nullitetra_dom_long$variable == "perc_nulliB_mismap" | 
                                                 leaf_nullitetra_dom_long$variable == "perc_nulliD_mismap"),]
head(A.supp.non.supp)
A.supp.non.supp$general_description <- "Supp_non_supp"
head(A.supp.non.supp)
B.supp.non.supp <- leaf_nullitetra_dom_long[leaf_nullitetra_dom_long$description == "B.suppressed" & 
                                              (leaf_nullitetra_dom_long$variable == "perc_nulliA_mismap" | 
                                                 leaf_nullitetra_dom_long$variable == "perc_nulliD_mismap"),]
head(B.supp.non.supp)
B.supp.non.supp$general_description <- "Supp_non_supp"
D.supp.non.supp <- leaf_nullitetra_dom_long[leaf_nullitetra_dom_long$description == "D.suppressed" & 
                                              (leaf_nullitetra_dom_long$variable == "perc_nulliB_mismap" | 
                                                 leaf_nullitetra_dom_long$variable == "perc_nulliA_mismap"),]
head(D.supp.non.supp)
D.supp.non.supp$general_description <- "Supp_non_supp"

# now add in the dom_non_dom and sup_non_sup
dom_sup_central_non_dom_non_sup <- rbind(dom_supp_central,A.dom.non.dom, B.dom.non.dom, D.dom.non.dom,
      A.supp.non.supp,B.supp.non.supp,D.supp.non.supp)
head(dom_sup_central_non_dom_non_sup)
tail(dom_sup_central_non_dom_non_sup)

# plot graph 
library(ggplot2)
pdf(file="leaf_triad_expr_in_nullis_no_zeroCS_dominance_summary_with_non_dom_non_sup.pdf")
ggplot(data=dom_sup_central_non_dom_non_sup, aes(x= general_description,y=value)) + 
  geom_boxplot(outlier.shape = NA)  +
  coord_cartesian(ylim = c(0,20)) + 
  scale_y_continuous(breaks=seq(0, 20, 2)) + theme_bw() +
  ylab("Mismapping (%)") + xlab(NULL)  
dev.off()

# calc summary stats
summarised_leaf_dom_non <- dom_sup_central_non_dom_non_sup %>%
  group_by(general_description) %>%
  summarise(mean = mean(value),
            Q1 = quantile(value,0.25),
            Q2 = quantile(value,0.5),
            Q3= quantile(value,0.75),
            percentile90 = quantile(value,0.9),
            percentile95 = quantile(value,0.95),
            num_genes = length(value))

as.data.frame(summarised_leaf_dom_non)

write.csv(file="leaf_triad_expr_in_nullis_no_zeroCS_dominance_summary_with_non_dom_non_sup.csv",summarised_leaf_dom_non, row.names = F)

kruskal.test(value ~ general_description, data=dom_sup_central_non_dom_non_sup)

### Dunn test
#install.packages("FSA")
library(FSA)
PT = dunnTest(value ~ general_description,
              data=dom_sup_central_non_dom_non_sup,
              method="by")
PT




#### now do the same for root ####

root_nullitetra_dom <- read.csv(file="root_triad_expr_in_nullis_no_zeroCS_with_dominance.csv")
# select only triads which are >1 tpm in CS
head(root_nullitetra_dom)
dim(root_nullitetra_dom)

root_nullitetra_dom <- root_nullitetra_dom[(root_nullitetra_dom$CS_tpmA + root_nullitetra_dom$CS_tpmB + root_nullitetra_dom$CS_tpmD) > 1,]
head(root_nullitetra_dom)
dim(root_nullitetra_dom)

t.test(root_nullitetra_dom$perc_nulliA_mismap, root_nullitetra_dom$perc_nulliB_mismap)
t.test(root_nullitetra_dom$perc_nulliA_mismap, root_nullitetra_dom$perc_nulliD_mismap)

# first make graphs not split by dominance category

root_nullitetra_long <- melt(root_nullitetra_dom[,c(1,14:16)],
                             id.vars = c("group_id"))
head(root_nullitetra_long)
tail(root_nullitetra_long)

kruskal.test(value ~ variable, data=root_nullitetra_long)

### Dunn test
#install.packages("FSA")
library(FSA)
PT = dunnTest(value ~ variable,
              data=root_nullitetra_long,
              method="by")
PT

# plot graph 
library(ggplot2)
pdf(file="root_triad_expr_in_nullis_no_zeroCS_1tpm_threshold.pdf")
# plot boxplot without the outliers, zoom using coord cartesian (doesnt change the data just zooms)
ggplot(data=root_nullitetra_long, aes(x= variable,y=value)) + 
  geom_boxplot(outlier.shape = NA)  +
  coord_cartesian(ylim = c(0,12)) + 
  scale_y_continuous(breaks=seq(0, 12, 2)) + theme_bw() +
  ylab("Mismapping (%)") + xlab(NULL) +
  scale_x_discrete(breaks=c("perc_nulliA_mismap", "perc_nulliB_mismap", "perc_nulliD_mismap"),
                   labels=c("nulli A", "nulli B", "nulli D")) 
dev.off()

# calc summary stats
library(tidyr)
library(dplyr)
summarised_root <- root_nullitetra_long %>%
  group_by(variable) %>%
  summarise(mean = mean(value),
            Q1 = quantile(value,0.25),
            Q2 = quantile(value,0.5),
            Q3= quantile(value,0.75),
            percentile90 = quantile(value,0.9),
            percentile95 = quantile(value,0.95),
            num_triads = length(value))

as.data.frame(summarised_root)
write.csv(file="root_triad_expr_in_nullis_no_zeroCS_1tpm_threshold_stats.csv",summarised_root, row.names = F)


# convert to long format for ggplot2
head(root_nullitetra_dom)
root_nullitetra_dom_long <- melt(root_nullitetra_dom[,c(1,14:19)],
                                 id.vars = c("group_id","description","general_description"))
head(root_nullitetra_dom_long)


# now get A dom in nulliA, B dom in nulliB, D dom in nulliD
A.dom <- root_nullitetra_dom_long[root_nullitetra_dom_long$description == "A.dominant" & root_nullitetra_dom_long$variable == "perc_nulliA_mismap",]
head(A.dom)
B.dom <- root_nullitetra_dom_long[root_nullitetra_dom_long$description == "B.dominant" & root_nullitetra_dom_long$variable == "perc_nulliB_mismap",]
head(B.dom)
D.dom <- root_nullitetra_dom_long[root_nullitetra_dom_long$description == "D.dominant" & root_nullitetra_dom_long$variable == "perc_nulliD_mismap",]
head(D.dom)

dom <- rbind(A.dom, B.dom, D.dom)
head(dom)
tail(dom)
dim(dom)

# now get A supp in nulliA, B supp in nulliB, D supp in nulliD
A.supp <- root_nullitetra_dom_long[root_nullitetra_dom_long$description == "A.suppressed" & root_nullitetra_dom_long$variable == "perc_nulliA_mismap",]
head(A.supp)
B.supp <- root_nullitetra_dom_long[root_nullitetra_dom_long$description == "B.suppressed" & root_nullitetra_dom_long$variable == "perc_nulliB_mismap",]
head(B.supp)
D.supp <- root_nullitetra_dom_long[root_nullitetra_dom_long$description == "D.suppressed" & root_nullitetra_dom_long$variable == "perc_nulliD_mismap",]
head(D.supp)

supp <- rbind(A.supp, B.supp, D.supp)
head(supp)
tail(supp)
dim(supp)

# now get Central
head(root_nullitetra_dom)
central <- root_nullitetra_dom_long[root_nullitetra_dom_long$description == "Central" & root_nullitetra_dom_long$variable != "av_perc_nulli_mismap",]
head(central)
tail(central)
dim(central)

dom_supp_central <- rbind(dom, supp, central)
head(dom_supp_central)
tail(dom_supp_central)

# plot graph 
library(ggplot2)
pdf(file="root_triad_expr_in_nullis_no_zeroCS_dominance_summary.pdf")
ggplot(data=dom_supp_central, aes(x= general_description,y=value)) + 
  geom_boxplot(outlier.shape = NA)  +
  coord_cartesian(ylim = c(0,20)) + 
  scale_y_continuous(breaks=seq(0, 20, 2)) + theme_bw() +
  ylab("Mismapping (%)") + xlab(NULL)  
dev.off()

# calc summary stats
summarised_root_dom <- dom_supp_central %>%
  group_by(general_description) %>%
  summarise(mean = mean(value),
            Q1 = quantile(value,0.25),
            Q2 = quantile(value,0.5),
            Q3= quantile(value,0.75),
            percentile90 = quantile(value,0.9),
            percentile95 = quantile(value,0.95),
            num_genes = length(value))

as.data.frame(summarised_root_dom)

write.csv(file="root_triad_expr_in_nullis_no_zeroCS_dominance_summary.csv",summarised_root_dom, row.names = F)

kruskal.test(value ~ general_description, data=dom_supp_central)

### Dunn test
#install.packages("FSA")
library(FSA)
PT = dunnTest(value ~ general_description,
              data=dom_supp_central,
              method="by")
PT

# now want to also include dominant_non_dominant and suppressed_non_suppressed
A.dom.non.dom <- root_nullitetra_dom_long[root_nullitetra_dom_long$description == "A.dominant" & 
                                            (root_nullitetra_dom_long$variable == "perc_nulliB_mismap" | 
                                               root_nullitetra_dom_long$variable == "perc_nulliD_mismap"),]
head(A.dom.non.dom)
A.dom.non.dom$general_description <- "Dom_non_dom"
head(A.dom.non.dom)

B.dom.non.dom <- root_nullitetra_dom_long[root_nullitetra_dom_long$description == "B.dominant" &  
                                            (root_nullitetra_dom_long$variable == "perc_nulliA_mismap" | 
                                               root_nullitetra_dom_long$variable == "perc_nulliD_mismap"),]
head(B.dom.non.dom)
B.dom.non.dom$general_description <- "Dom_non_dom"
D.dom.non.dom <- root_nullitetra_dom_long[root_nullitetra_dom_long$description == "D.dominant" &  
                                            (root_nullitetra_dom_long$variable == "perc_nulliB_mismap" | 
                                               root_nullitetra_dom_long$variable == "perc_nulliA_mismap"),]
head(D.dom.non.dom)
D.dom.non.dom$general_description <- "Dom_non_dom"

# now get A supp in nulliA, B supp in nulliB, D supp in nulliD
A.supp.non.supp <- root_nullitetra_dom_long[root_nullitetra_dom_long$description == "A.suppressed" & 
                                              (root_nullitetra_dom_long$variable == "perc_nulliB_mismap" | 
                                                 root_nullitetra_dom_long$variable == "perc_nulliD_mismap"),]
head(A.supp.non.supp)
A.supp.non.supp$general_description <- "Supp_non_supp"
head(A.supp.non.supp)
B.supp.non.supp <- root_nullitetra_dom_long[root_nullitetra_dom_long$description == "B.suppressed" & 
                                              (root_nullitetra_dom_long$variable == "perc_nulliA_mismap" | 
                                                 root_nullitetra_dom_long$variable == "perc_nulliD_mismap"),]
head(B.supp.non.supp)
B.supp.non.supp$general_description <- "Supp_non_supp"
D.supp.non.supp <- root_nullitetra_dom_long[root_nullitetra_dom_long$description == "D.suppressed" & 
                                              (root_nullitetra_dom_long$variable == "perc_nulliB_mismap" | 
                                                 root_nullitetra_dom_long$variable == "perc_nulliA_mismap"),]
head(D.supp.non.supp)
D.supp.non.supp$general_description <- "Supp_non_supp"

# now add in the dom_non_dom and sup_non_sup
dom_sup_central_non_dom_non_sup <- rbind(dom_supp_central,A.dom.non.dom, B.dom.non.dom, D.dom.non.dom,
                                         A.supp.non.supp,B.supp.non.supp,D.supp.non.supp)
head(dom_sup_central_non_dom_non_sup)
tail(dom_sup_central_non_dom_non_sup)

# plot graph 
library(ggplot2)
pdf(file="root_triad_expr_in_nullis_no_zeroCS_dominance_summary_with_non_dom_non_sup.pdf")
ggplot(data=dom_sup_central_non_dom_non_sup, aes(x= general_description,y=value)) + 
  geom_boxplot(outlier.shape = NA)  +
  coord_cartesian(ylim = c(0,20)) + 
  scale_y_continuous(breaks=seq(0, 20, 2)) + theme_bw() +
  ylab("Mismapping (%)") + xlab(NULL)  
dev.off()

# calc summary stats
summarised_root_dom_non <- dom_sup_central_non_dom_non_sup %>%
  group_by(general_description) %>%
  summarise(mean = mean(value),
            Q1 = quantile(value,0.25),
            Q2 = quantile(value,0.5),
            Q3= quantile(value,0.75),
            percentile90 = quantile(value,0.9),
            percentile95 = quantile(value,0.95),
            num_genes = length(value))

as.data.frame(summarised_root_dom_non)

write.csv(file="root_triad_expr_in_nullis_no_zeroCS_dominance_summary_with_non_dom_non_sup.csv",summarised_root_dom_non, row.names = F)

kruskal.test(value ~ general_description, data=dom_sup_central_non_dom_non_sup)

### Dunn test
#install.packages("FSA")
library(FSA)
PT = dunnTest(value ~ general_description,
              data=dom_sup_central_non_dom_non_sup,
              method="by")
PT

