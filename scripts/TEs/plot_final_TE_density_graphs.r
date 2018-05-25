library(ggplot2)
library(reshape2)
setwd("W:/Jemima/companion_paper/TEs/")

rolling_counts <- read.table("W:/Jemima/companion_paper/TEs/TE_5kb_upstream_rolling_window_100window_10step.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
no_stress_triads <- read.table("W:/Jemima/companion_paper/HC_CS_no_stress_triads_with_five_cats.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
HC_CS_low10 <- read.table("W:/expression_browser/dataset_reports/20170911_all_reports/02.movement/HC_CS_no_stress/HC_CS_no_stress_movement_low_10pc.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
HC_CS_top10 <- read.table("W:/expression_browser/dataset_reports/20170911_all_reports/02.movement/HC_CS_no_stress/HC_CS_no_stress_movement_top_10pc.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
HC_CS_middle80 <- read.table("W:/expression_browser/dataset_reports/20170911_all_reports/02.movement/HC_CS_no_stress/HC_CS_no_stress_movement_middle_80pc.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)


#colours
group.colors<-c(A.dominant = "#579D1C", B.dominant = "#4B1F6F", D.dominant ="#FF950E", Central="#AAAAAA",Balanced="#AAAAAA", A.suppressed = "#b2E08a", B.suppressed = "#0eC7ff", D.suppressed ="#ffCF0e")
chr.colors<-c(A = "#579D1C", B = "#4B1F6F", D ="#FF950E")
movment.colors = c("Dynamic"="#e41a1c", "Middle 80"="#bdbdbd", "Stable"="#377eb8")
five.cat.colors <- c("non_suppressed" = "#a6611a", "suppressed" = "#dfc27d", "central" = "#AAAAAA", "non_dominant" = "#80cdc1", "dominant" = "#018571")


#plot for the five categories
no_stress_five <- no_stress_triads[,c("group_id", "gene", "five_cat")]
colnames(no_stress_five)[3] <- "description"

allocate_descriptions_mean <- function(data, description_cats, description){
	data_desc <- subset(data, gene %in% description_cats$gene)
	data_desc_means <- aggregate(data_desc$roll_count, list(data_desc$midpoint), mean)
	colnames(data_desc_means) <- c("midpoint", "pc_TE")
	data_desc_means$description <- description
	return(data_desc_means)
}

plot_subsets_colours <- function(rolling_counts, triad_list, plot_name, plot_colours){ 
	no_stress_only_roll <- subset(rolling_counts, gene %in% triad_list$gene)
	no_stress_rolling_means <- data.frame(midpoint = numeric(), pc_TE = numeric(), description = character())
	
	descriptions <- unique(triad_list$description)
	
	for (i in seq(1, length(descriptions))){
		description <- descriptions[i]
		description_cats <- triad_list[triad_list$description == description,]
		print(description)
		data_temp_mean <- allocate_descriptions_mean(data = no_stress_only_roll, description_cats = description_cats, description = description)
		no_stress_rolling_means <- rbind(no_stress_rolling_means, data_temp_mean)
		gc()
	}
	
	#plot these categories
	no_stress_plot <- ggplot(no_stress_rolling_means, aes(x=midpoint, y = pc_TE, fill = description, colour = description)) +
						geom_line(lwd = 1) +
						ggtitle(plot_name) +
						scale_colour_manual(values = plot_colours) +
						scale_fill_manual(values = plot_colours) +
						theme_bw()
	ggsave(no_stress_plot, file = paste(plot_name, ".pdf", sep = ""), height = 10, width = 10)
	return(no_stress_rolling_means)
}


no_stress_five_means <- plot_subsets_colours(rolling_counts = rolling_counts, triad_list = no_stress_five, plot_name = "TE_density_five_cats_HC_CS_no_stress", plot_colours = five.cat.colors)

plot_subset_rolling_only <- function(rolling_means, plot_name, plot_colours){
	plot <- ggplot(rolling_means, aes(x=midpoint, y = pc_TE, fill = description, colour = description)) +
						geom_line(lwd = 1) +
						ggtitle(plot_name) +
						scale_colour_manual(values = plot_colours) +
						scale_fill_manual(values = plot_colours) +
						theme_bw()
	ggsave(plot, file = paste(plot_name, ".pdf", sep = ""), height = 10, width = 10)
	return(plot)
}


five_means_roll_plot <- plot_subset_rolling_only(no_stress_five_means, plot_name = "TE_density_five_cats_HC_CS_no_stress", plot_colours = five.cat.colours)


#now do the p.value plots for the five cats
wilcox_midpoints_description <- function(test_data){
	test_data$description <- as.factor(test_data$description)
	
	midpoints <- unique(test_data$midpoint)
	
	Dynamic_Stable <- numeric()
	mid80_Stable <- numeric()
	Dynamic_mid80 <- numeric()
	
	for (i in seq(1, length(midpoints))){
		midpoint <- midpoints[i]
		midpoint_data <- test_data[test_data$midpoint == midpoint,]
		wilcox_midpoint <- pairwise.wilcox.test(x = midpoint_data$roll_count, g = midpoint_data$description, p.adjust.method = "BH")
		Dynamic_Stable[i] <- wilcox_midpoint$p.value["Stable", "Dynamic"]
		mid80_Stable[i] <- wilcox_midpoint$p.value["Stable", "Middle 80"]
		Dynamic_mid80[i] <- wilcox_midpoint$p.value["Middle 80", "Dynamic"]
		#top10_low10[i] <- midpoint
		#mid80_low80[i] <- midpoint
		#top10_mid80[i] <- midpoint
	}
	
	midpoints_wilcox <- data.frame(midpoints, Dynamic_Stable, mid80_Stable, Dynamic_mid80)
	return(midpoints_wilcox)
}

wilcox_midpoints_description_five_cat <- function(test_data){
	test_data$description_5cat <- as.factor(test_data$description_5cat)
	
	midpoints <- unique(test_data$midpoint)
	
	dom_central <- numeric()
	non.dom_central <- numeric()
	non.sup_central	<- numeric()
	sup_central <- numeric()
	non.dom_dom <- numeric()
	non.sup_dom	<- numeric()
	sup_dom <- numeric()
	non.sup_non.dom	<- numeric()
	sup_non.dom <- numeric()
	sup_non.sup <- numeric()
	
	for (i in seq(1, length(midpoints))){
		midpoint <- midpoints[i]
		midpoint_data <- test_data[test_data$midpoint == midpoint,]
		wilcox_midpoint <- pairwise.wilcox.test(x = midpoint_data$roll_count, g = midpoint_data$description_5cat, p.adjust.method = "BH")
		dom_central[i] <- wilcox_midpoint$p.value["dominant", "central"]
		non.dom_central[i] <- wilcox_midpoint$p.value["non_dominant", "central"]
		non.sup_central[i]	<- wilcox_midpoint$p.value["non_suppressed", "central"]
		sup_central[i] <- wilcox_midpoint$p.value["suppressed", "central"]
		non.dom_dom[i] <- wilcox_midpoint$p.value["non_dominant", "dominant"]
		non.sup_dom[i]	<- wilcox_midpoint$p.value["non_suppressed", "dominant"]
		sup_dom[i] <- wilcox_midpoint$p.value["suppressed", "dominant"]
		non.sup_non.dom[i]	<- wilcox_midpoint$p.value["non_suppressed", "non_dominant"]
		sup_non.dom[i] <- wilcox_midpoint$p.value["suppressed", "non_dominant"]
		sup_non.sup[i] <- wilcox_midpoint$p.value["suppressed", "non_suppressed"]
	}
	
	midpoints_wilcox <- data.frame(midpoints, dom_central, non.dom_central, non.sup_central, sup_central, non.dom_dom, non.sup_dom, sup_dom, non.sup_non.dom, sup_non.dom, sup_non.sup)
	return(midpoints_wilcox)
}

rolling_counts$description_5cat <- NA

five_cats <- unique(no_stress_five$description)

for (i in seq(1, length(five_cats))){
	five_cat <- five_cats[i]
	five_cat_genes <- no_stress_five[no_stress_five$description == five_cat, "gene"]
	
	rolling_counts[which(rolling_counts$gene %in% five_cat_genes ),"description_5cat"] <- five_cat
}

rolling_counts_fivecat <- rolling_counts[complete.cases(rolling_counts$description_5cat), ]

fivecat_rolling_wilcox_all <- wilcox_midpoints_description_five_cat(rolling_counts_fivecat)
fivecat_rolling_wilcox_all_melt <- melt(fivecat_rolling_wilcox, id = "midpoints")

pvalue_5_plot <- ggplot(fivecat_rolling_wilcox_melt, aes(x = midpoints, y = value, fill = variable, colour = variable)) +
				geom_line(lwd = 0.5) + 
				#coord_cartesian(ylim = c(0, 0.1)) +
				geom_hline(yintercept = 0.01) +
				scale_colour_manual(values = five.cat.colors) +
				scale_fill_manual(values = five.cat.colors) +
				theme_bw()
				
ggsave(pvalue_5_plot, file = "TE_density_five_cats_HC_CS_no_stress_wilcoxon_pvalue.pdf")

###### now repeat all of this for the movement categories
rolling_counts$description <- NA

rolling_counts[which(rolling_counts$gene %in% HC_CS_low10$x),"description"] <- "Stable"
rolling_counts[which(rolling_counts$gene %in% HC_CS_top10$x),"description"] <- "Dynamic"
rolling_counts[which(rolling_counts$gene %in% HC_CS_middle80$x),"description"] <- "Middle 80"

no_stress_movement <- unique(rolling_counts[,c("gene", "description")])
no_stress_movement <- no_stress_movement[complete.cases(no_stress_movement),]

no_stress_movement_means <- plot_subsets_colours(rolling_counts = rolling_counts, triad_list = no_stress_movement, plot_name = "TE_density_movement_HC_CS_no_stress", plot_colours = movment.colors)

no_stress_movement_means <- plot_subsets_colours(rolling_counts = rolling_counts, triad_list = no_stress_movement, plot_name = "TE_density_movement_HC_CS_no_stress", plot_colours = movment.colors)

#now do the p.value plots for movement
rolling_counts_movement <- rolling_counts[complete.cases(rolling_counts$description), ]

movement_rolling_wilcox_all <- wilcox_midpoints_description(rolling_counts_movement)
movement_rolling_wilcox_all_melt <- melt(movement_rolling_wilcox_all, id = "midpoints")

pvalue_move_plot <- ggplot(movement_rolling_wilcox_all_melt, aes(x = midpoints, y = value, fill = variable, colour = variable)) +
				geom_line(lwd = 0.5) + 
				#coord_cartesian(ylim = c(0, 0.1)) +
				geom_hline(yintercept = 0.01) +
				#scale_colour_manual(values = five.cat.colors) +
				#scale_fill_manual(values = five.cat.colors) +
				theme_bw()
				
ggsave(pvalue_move_plot, file = "TE_density_movement_HC_CS_no_stress_wilcoxon_pvalue.pdf", height = 10, width = 10)
pairwise_colours <- c("Dynamic_Stable"="#A329B3", "mid80_Stable" = "#959BE7", "Dynamic_mid80" = "#E79895")

pvalue_move_plot_log <- ggplot(movement_rolling_wilcox_all_melt, aes(x = midpoints, y = log10(value), fill = variable, colour = variable)) +
				geom_line(lwd = 0.5) + 
				geom_hline(yintercept = log10(0.01), linetype = "dashed") +
				scale_colour_manual(values = pairwise_colours) +
				theme_bw() +
				theme(legend.position = "none")
ggsave(pvalue_move_plot_log, file = "./plots_to_use/TE_density_movement_HC_CS_no_stress_wilcoxon_pvalue_log.pdf", height = 5, width = 10, units = "cm")

write.table(movement_rolling_wilcox_all, file = "TE_density_movement_CS_no_stress_wilcox_pvalues.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
				
