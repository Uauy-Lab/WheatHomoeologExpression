# Aim is to plot the results from homoeologue allocation to modules
# 15.8.2017
# Philippa Borrill


list_dirs <- c("tissues\\grain","tissues\\leaf","tissues\\root","tissues\\spike")

list_names <- c("tissues\\grain" = "grain","tissues\\leaf" = "leaf","tissues\\root" = "root","tissues\\spike" = "spike")

out_dir <- "Y:\\expression_browser\\WGA\\WGCNA\\homoeologue_paralogue_all_networks\\same_near_far_for_manuscript\\"

library(reshape2)

for (dir in list_dirs) {
  
  setwd(paste0("Y:\\expression_browser\\WGA\\WGCNA\\",dir,"\\maxP0.05\\homoeologues_paralogues_numeric_threshold_correct"))
  
triad_1_1_1_synt <- read.csv(file="median_eigengene_distances_to_nearest_x_eigengenes_homoeologue_movements_1_1_1_syntenic.csv")
triad_1_1_1_non_synt <- read.csv(file="median_eigengene_distances_to_nearest_x_eigengenes_homoeologue_movements_1_1_1_non-syntenic.csv")
triad_random <- read.csv(file="median_eigengene_distances_to_nearest_x_eigengenes_homoeologue_movements_random.csv")

head(triad_1_1_1_synt)
head(triad_1_1_1_non_synt)
head(triad_random)
colnames(triad_1_1_1_non_synt) == colnames(triad_random)

# combine datasets
triad <- rbind(triad_1_1_1_synt,triad_1_1_1_non_synt,triad_random)
head(triad)
tail(triad)


triad_formatted <- triad[,c(1,2,3,7,9,10,12,13,14)]
head(triad_formatted)

colnames(triad_formatted)[4] <- "two_hom"
colnames(triad_formatted)[6] <- "perc_further_than_threshold"
colnames(triad_formatted)


df.all <- triad_formatted[,-1]
head(df.all)

# now calculate averages for different sets
library(plyr)
syntenic <- ddply(df.all[df.all$synteny == "syntenic",], 
                  c("threshold","identifier"),summarise, 
                  sum_two_hom= sum(two_hom),
                  sum_total_num_triads = sum(total_num_triads))
syntenic$set <- "syntenic_1_1_1"
head(syntenic)

non_syntenic <- ddply(df.all[df.all$synteny == "non-syntenic",], 
                     c("threshold","identifier"),summarise, 
                     sum_two_hom= sum(two_hom),
                     sum_total_num_triads = sum(total_num_triads))
non_syntenic$set <- "non-syntenic_1_1_1"
head(non_syntenic)


random <- ddply(df.all[df.all$type == "random" ,], 
                   c("threshold","identifier"),summarise, 
                   sum_two_hom= sum(two_hom),
                   sum_total_num_triads = sum(total_num_triads))
random$set <- "random"
head(random)

# now combine these averages into one df for plotting
combined.df <- rbind(syntenic, non_syntenic, random)
head(combined.df)

combined.df$under_threshold <- combined.df$sum_total_num_triads - combined.df$sum_two_hom 
head(combined.df)
tail(combined.df)



# now make into long format to plot in ggplot
combined.df$perc_over_threshold <- combined.df$sum_two_hom/combined.df$sum_total_num_triads
combined.df$perc_under_threshold <- combined.df$under_threshold/combined.df$sum_total_num_triads
head(combined.df)
is.numeric(combined.df$perc_over_threshold)

unique(combined.df$threshold)

write.csv(file=paste0(out_dir, list_names[dir], "_combined_results_all_thresholds.csv"), combined.df)

melted_combined.df <- melt(combined.df[,c(1,5,7,8)], id=c("threshold","set"))
head(melted_combined.df)
tail(melted_combined.df)
unique(melted_combined.df$variable)

is.numeric(melted_combined.df$threshold)
melted_combined.df$set <- as.character(melted_combined.df$set)


library(ggplot2)
plot_thresholds_2 <- ggplot(data=melted_combined.df, 
                          aes(x=threshold/max(threshold), y=value, colour=set, lty = variable)) + 
  geom_line(size=1) + ylab("Percentage of triads") + xlab("Percentage of distance to furthest eigengene")
pdf(file=paste0(out_dir, list_names[dir], "_percentages triads over percentage threshold.pdf"))
print(plot_thresholds_2)
dev.off()


### now want to make stacked bar charts ###
head(combined.df)

same_near_far <- combined.df[combined.df$identifier =="same",] # get percentage same (i.e. 100 %)
same_near_far

same_near_far$perc_over_near <- combined.df[combined.df$identifier == "0",7]  # get percentage further than same
same_near_far

same_near_far$perc_over_far <- combined.df[combined.df$identifier == "50",7] # get percentage further than nera (threshold 50)
same_near_far

same_near_far$perc_same <- same_near_far$perc_over_threshold - same_near_far$perc_over_near # perc same - perc further than same gives % same
same_near_far$perc_near <- same_near_far$perc_over_near - same_near_far$perc_over_far # perc over near - perc over far gives % near
same_near_far$perc_far <- same_near_far$perc_over_far # perc over far = perc far
same_near_far

same_near_far$num_same <- same_near_far$sum_total_num_triads * same_near_far$perc_same
same_near_far$num_near <- same_near_far$sum_total_num_triads * same_near_far$perc_near
same_near_far$num_far <- same_near_far$sum_total_num_triads * same_near_far$perc_far
same_near_far

same_near_far_csv <- same_near_far[,c(5,14:16,11:13,4)]
same_near_far_csv

write.csv(file=paste0(out_dir, list_names[dir], "_same_near_far.csv"),same_near_far_csv)

same_near_far_to_plot <- same_near_far[,c(5,11,12,13)]
same_near_far_to_plot

melted_same_near_far <- melt(same_near_far_to_plot, id.vars = "set")
head(melted_same_near_far)

melted_same_near_far$variable <- factor(melted_same_near_far$variable, levels= sort(levels(melted_same_near_far$variable)))
levels(melted_same_near_far$variable)

melted_same_near_far$set <- factor(melted_same_near_far$set, levels=c("syntenic_1_1_1","non-syntenic_1_1_1","random"))
levels(melted_same_near_far$set) 

library(ggplot2)
plot_bar_stacked_horiz <- ggplot(melted_same_near_far, aes(x=set, y=value, fill=variable)) + 
  geom_bar(position="stack", stat="identity")  +
  scale_fill_manual(values=c("#fdae61","#2c7bb6", "#abd9e9"),
                    name="Homoeologue module",
                    breaks=c("perc_far","perc_near","perc_same"),
                    labels=c("Far","Near","Same")) +
  ylab("Percentage of triads") + xlab(NULL) +theme_bw() + scale_y_continuous(expand=c(0,0)) +coord_flip()
plot_bar_stacked_horiz
pdf(file=paste0(out_dir, list_names[dir], "_bar_chart_50_horiz.pdf"), height=4)
print(plot_bar_stacked_horiz)
dev.off()


 
}


#### now plot one summary graph using "all" triads for each network ####


setwd(out_dir)

file<- "_same_near_far.csv"

grain <- read.csv(paste0("grain",file))
leaf <- read.csv(paste0("leaf",file))
spike <- read.csv(paste0("spike",file))
root <- read.csv(paste0("root",file))

# add column to say which network it came from
grain$network <- "grain"
leaf$network <- "leaf"
spike$network <- "spike"
root$network <- "root"

head(root)

all_data <- rbind(grain,leaf,spike,root)
head(all_data)

random_data <- all_data[all_data$set =="random" ,]
head(random_data)
dim(random_data)

random_summary <- apply(random_data[,c(6:8)], 2, mean)
head(random_summary)

random_summary_melted <- melt(random_summary)
random_summary_melted

random_summary_melted$network <- "random"
random_summary_melted$set <- "random"
random_summary_melted$variable <- rownames(random_summary_melted)
rownames(random_summary_melted) <- NULL
random_summary_melted

## make graph for syntenic:

summary_data <- all_data[all_data$set =="syntenic_1_1_1" |all_data$set =="non-syntenic_1_1_1"  ,]
dim(summary_data)
head(summary_data)

melted_summary_data <- melt(summary_data[,c(2,6,7,8,10)], id=c("network","set"))
head(melted_summary_data)
tail(melted_summary_data)
unique(melted_summary_data$variable)

melted_summary_data <- rbind(melted_summary_data, random_summary_melted)

melted_summary_data

melted_summary_data$variable <- factor(melted_summary_data$variable, levels= sort(levels(melted_summary_data$variable)))
levels(melted_summary_data$variable)

melted_summary_data$network <- factor(melted_summary_data$network, levels= rev(c("grain","spike","leaf","root","random")))
levels(melted_summary_data$network)



library(ggplot2)
plot_bar_stacked_horiz_all_synt <- ggplot(melted_summary_data[melted_summary_data$set != "non-syntenic_1_1_1",], aes(x=network, y=value, fill=variable)) + 
  geom_bar(position="stack", stat="identity")  +
  scale_fill_manual(values=c("#fdae61","#2c7bb6", "#abd9e9"),
                    name="Homoeologue module",
                    breaks=c("perc_far","perc_near","perc_same"),
                    labels=c("Far","Near","Same")) +
  ylab("Percentage of triads") + xlab(NULL) +theme_bw() + scale_y_continuous(expand=c(0,0)) +coord_flip()
plot_bar_stacked_horiz_all_synt
pdf(file=paste0(out_dir,  "all_networks_synt_1_1_1.pdf"), height=4)
print(plot_bar_stacked_horiz_all_synt)
dev.off()



plot_bar_stacked_horiz_all_non_synt <- ggplot(melted_summary_data[melted_summary_data$set != "syntenic_1_1_1",], aes(x=network, y=value, fill=variable)) + 
  geom_bar(position="stack", stat="identity")  +
  scale_fill_manual(values=c("#fdae61","#2c7bb6", "#abd9e9"),
                    name="Homoeologue module",
                    breaks=c("perc_far","perc_near","perc_same"),
                    labels=c("Far","Near","Same")) +
  ylab("Percentage of triads") + xlab(NULL) +theme_bw() + scale_y_continuous(expand=c(0,0)) +coord_flip()
plot_bar_stacked_horiz_all_non_synt
pdf(file=paste0(out_dir,  "all_networks_non_synt_1_1_1.pdf"), height=4)
print(plot_bar_stacked_horiz_all_non_synt)
dev.off()
