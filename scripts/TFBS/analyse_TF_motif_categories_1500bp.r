setwd("W:/Jemima/companion_paper/TF_motifs/motif_categories_analysis/1500bp_upstream")
library(ggplot2)

group.colors<-c(A.dominant = "#579D1C", B.dominant = "#4B1F6F", D.dominant ="#FF950E",
             Central="#AAAAAA",Balanced="#AAAAAA",
             A.suppressed = "#b2E08a", B.suppressed = "#0eC7ff", D.suppressed ="#ffCF0e")
chr.colors<-c(A = "#579D1C", B = "#4B1F6F", D ="#FF950E")


movment.colors = c("Dynamic"="#e41a1c", "Middle 80"="#bdbdbd", "Stable"="#377eb8")

five.cat.colors <- c("non_suppressed" = "#a6611a", "suppressed" = "#dfc27d", "central" = "#AAAAAA", "non_dominant" = "#80cdc1", "dominant" = "#018571")


unique_counts_path <- "W:/Jemima/companion_paper/out_2018_03_12_12_03/fimo_results/IWGSCv1.0_UTR.HC.firstCDS_consensus._1500bp_upstream_triads_only.fa/unique_hits_analysis/fimo_unique_TF_motif_counts_per_gene.txt"
no_stress_triads <- read.table("W:/Jemima/companion_paper/HC_CS_no_stress_triads_with_five_cats.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
unique_counts_path2 <- "W:/Jemima/companion_paper/out_2018_03_12_12_03/fimo_results/IWGSCv1.0_UTR.HC.firstCDS_consensus._1500bp_upstream_triads_only.fa/unique_hits_analysis/fimo_unique_TF_motifs_per_gene.txt"

unique_counts <- read.table(unique_counts_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
unique_counts2 <- read.table(unique_counts_path2, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

triads_motifs <- subset(triads, gene %in% unique(unique_counts$gene))

overall_plot <- ggplot(data = unique_counts, aes(unique_counts$unique_motifs)) +
					geom_density(fill = "gray", alpha = 0.3) +
					geom_vline(xintercept = median_motif, colour = "red") +
					theme_bw()
					


ggsave(overall_plot, file = "./plots_to_use/motif_count_all_genes_1500bp_to_use.pdf", height = 8, width = 8, units = "cm")


#plot the frequency of each individual motif in each promoter to see if presence/absence is relevant
motif_gene_counts <- read.table("W:/Jemima/companion_paper/TF_motifs/motif_gene_counts_1500bp_frequency.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
motif_gene_plot <- ggplot(motif_gene_counts, aes(x = no_of_counts, y = freq)) + 
							geom_bar(stat = 'identity') +
							theme_bw() +
							xlab("TFBS motif copy number") +
							ylab("Frequency")
							

ggsave(motif_gene_plot, file = "./plots_to_use/motif_gene_counts_1500bp_to_use.pdf", height = 8, width = 8, units = "cm")
ggsave(motif_gene_plot, file = "./plots_to_use/motif_gene_counts_1500bp_to_use.png", height = 8, width = 8, units = "cm", dpi = 600)

motif_gene_plot_log <- ggplot(motif_gene_counts, aes(x = no_of_counts, y = log10(freq))) + 
							geom_bar(stat = 'identity') +
							theme_bw() +
							xlab("TFBS motif occurence") +
							ylab("log10(Frequency)")
ggsave(motif_gene_plot_log, file = "./plots_to_use/motif_gene_counts_1500bp_to_use_log.png", height = 8, width = 8, units = "cm", dpi = 600)
							
							
#now plot without the 0 category
motif_gene_counts_gt0 <- motif_gene_counts[!motif_gene_counts$no_of_counts == 0,]

motif_gene_plot_gt0 <- ggplot(motif_gene_counts_gt0, aes(x = no_of_counts, y = freq)) + 
						geom_bar(stat = 'identity') +
						theme_bw() +
						xlab("TFBS motif copy number") +
						ylab("Frequency")
ggsave(motif_gene_plot_gt0, file = "./plots_to_use/motif_gene_counts_1500bp_no_0_counts.pdf", height = 8, width = 8, units = "cm")
ggsave(motif_gene_plot_gt0, file = "./plots_to_use/motif_gene_counts_1500bp_no_0_counts.png", height = 8, width = 8, units = "cm", dpi = 600)

#now if we split into the different genomes
no_stress_triad_counts <- merge(no_stress_triads, unique_counts, all.x = TRUE, all.y = FALSE)

					
genome_plot <- ggplot(data = no_stress_triad_counts , aes(x=no_stress_triad_counts $unique_motifs, fill = chr_group)) +
					geom_density(alpha = 0.3) +
					scale_fill_manual(values = chr.colors)
					
ggsave(genome_plot, file = "motif_count_no_stress_triads_1500bp_ABD.pdf")

#test the 3 genomes individually
unique_motif_medians <- aggregate(unique_motifs ~ chr_group, data = no_stress_triad_counts, FUN = "median", na.omit = TRUE)
no_stress_triad_counts$chr_group <- as.factor(no_stress_triad_counts$chr_group)
kruskal.test(unique_motifs ~ chr_group, data = no_stress_triad_counts)
pairwise.wilcox.test(x = no_stress_triad_counts$unique_motifs, g = no_stress_triad_counts$chr_group, p.adjust.method = "BH")

##now do the movement categories
HC_CS_low10 <- read.table("W:/expression_browser/dataset_reports/20170911_all_reports/02.movement/HC_CS_no_stress/HC_CS_no_stress_movement_low_10pc.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
HC_CS_top10 <- read.table("W:/expression_browser/dataset_reports/20170911_all_reports/02.movement/HC_CS_no_stress/HC_CS_no_stress_movement_top_10pc.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
HC_CS_middle80 <- read.table("W:/expression_browser/dataset_reports/20170911_all_reports/02.movement/HC_CS_no_stress/HC_CS_no_stress_movement_middle_80pc.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#plot the number of unique motifs based on these groups
HC_CS_low10$description <- "Stable"
HC_CS_top10$description <- "Dynamic"
HC_CS_middle80$description <- "Middle 80"


movement_cats <- unique(rbind(HC_CS_low10, HC_CS_top10, HC_CS_middle80))
colnames(movement_cats)[1] <- "gene"


unique_counts$description <- NA

unique_counts[which(unique_counts$gene %in% HC_CS_low10$x),"description"] <- "Stable"
unique_counts[which(unique_counts$gene %in% HC_CS_top10$x),"description"] <- "Dynamic"
unique_counts[which(unique_counts$gene %in% HC_CS_middle80$x),"description"] <- "Middle 80"

movement_counts <- unique_counts[complete.cases(unique_counts$description),]
#including 0s
movement_counts_all <- merge(movement_cats, unique_counts, all.x = TRUE, all.y = FALSE)
movement_counts_all[is.na(movement_counts_all$unique_motifs),"unique_motifs"] <- 0

movement_plot <- ggplot(data = movement_counts , aes(x=movement_counts$unique_motifs, fill = description)) +
					geom_density(alpha = 0.3) +
					scale_fill_manual(values = movment.colors)
					
ggsave(movement_plot, file = "unique_motif_count_1500bp_movement_cats_no0genes.pdf")

movement_plot0 <- ggplot(data = movement_counts_all , aes(x=movement_counts_all$unique_motifs, fill = description)) +
					geom_density(alpha = 0.3) +
					scale_fill_manual(values = movment.colors)
					
ggsave(movement_plot0, file = "unique_motif_count_1500bp_movement_cats.pdf")

#stats on the numbers for the movement cats
sink("movement_cats_unique_motifs_medians+stats.txt")
print("stats using only genes that have motifs called")
movement_counts_medians <- aggregate(unique_motifs ~ description, data = movement_counts, FUN = "median", na.omit = TRUE)
print(movement_counts_medians)
movement_counts$description <- as.factor(movement_counts$description)
print(kruskal.test(unique_motifs ~ description, data = movement_counts))
print(pairwise.wilcox.test(x = movement_counts$unique_motifs, g = movement_counts$description, p.adjust.method = "BH"))

print("stats using genes that have motifs called in addition to genes with no motifs")
print("stats using only genes that have motifs called")
movement_counts_all_medians <- aggregate(unique_motifs ~ description, data = movement_counts_all, FUN = "median", na.omit = TRUE)
print(movement_counts_all_medians)
movement_counts_all$description <- as.factor(movement_counts_all$description)
print(kruskal.test(unique_motifs ~ description, data = movement_counts_all))
print(pairwise.wilcox.test(x = movement_counts_all$unique_motifs, g = movement_counts_all$description, p.adjust.method = "BH"))

sink()
#now have a look at the simple categories
triad_groups <- unique(triads[,c("group_id", "gene")])
movement_triads <- merge(triad_groups, movement_cats, all.x = FALSE, all.y = TRUE)

movement_triads_groups <- unique(movement_triads[,c("group_id", "description")])
simple_categs <- read.table("W:/Jemima/companion_paper/TF_motifs/triad_categories_per_triad_stacked_simple_categories_no_none_1500bp.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
motif_categs <- read.table("/Jemima/companion_paper/TF_motifs/triad_categories_per_triad_stacked_1500bp.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
simple_categs_all <- read.table("W:/Jemima/companion_paper/TF_motifs/triad_categories_per_triad_stacked_simple_categories_1500bp.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

categs_move <- merge(movement_triads_groups, simple_categs)

categories <- unique(categs_move$simple_category)
## remember to normalise for the no none categories
motif_categ_sums <- aggregate(categs_move$count, list(categs_move$group_id), sum)
colnames(motif_categ_sums) <- c("group_id", "total_count")

categs_move_norm <- merge(categs_move, motif_categ_sums)
categs_move_norm$count_norm <- categs_move_norm$count/categs_move_norm$total_count
categs_move_norm$simple_category <- factor(categs_move_norm$simple_category, levels = c("all_same", "A_diff", "B_diff", "D_diff"))
categs_move_norm$description <- factor(categs_move_norm$description, levels = c("Stable", "Middle 80", "Dynamic"))


categs_all_move <- merge(movement_triads_groups, simple_categs_all)
categs_all_move$simple_category <- factor(categs_all_move$simple_category, levels = c("all_same", "A_diff", "B_diff", "D_diff"))

for (i in seq(1, length(categories[i]))){
	categ <- categories[i]
	categ_data <- categs_move_norm[categs_move_norm$simple_category == categ,]
	cat_plot <- ggplot(data = categ_data, aes(x=count_norm, fill = description)) +
				geom_density(alpha = 0.5) +
				scale_fill_manual(values = movment.colors) +
				theme_bw()
				
	cat_plot <- ggplot(data = categ_data, aes(y=count_norm, x = description, fill = description)) +
				geom_boxplot() +
				scale_fill_manual(values = movment.colors) +
				theme_bw()
				
}

#plot all on one graphs

all_cat_boxplot_norm <- ggplot(data = categs_move_norm, aes(y=count_norm, x = simple_category, fill = description)) +
				geom_boxplot() +
				scale_fill_manual(values = movment.colors) +
				theme_bw()
				
ggsave(all_cat_boxplot_norm, file = "movement_categories_simple_TFBS_cats_normalised_counts.pdf")
				
all_cat_boxplot <- ggplot(data = categs_move_norm, aes(y=count, x = simple_category, fill = description)) +
				geom_boxplot() +
				scale_fill_manual(values = movment.colors) +
				theme_bw()
ggsave(all_cat_boxplot, file = "movement_categories_simple_TFBS_cats_counts.pdf", height = 5, width = 5)
				
all_cat_boxplot_no_legend <- ggplot(data = categs_move_norm, aes(y=count, x = simple_category, fill = description)) +
				geom_boxplot() +
				scale_fill_manual(values = movment.colors) +
				theme_bw() +
				theme(legend.position = "none") +
				xlab("TFBS motif category") +
				ylab("Frequency")
				
				
ggsave(all_cat_boxplot_no_legend, file = "movement_categories_simple_TFBS_cats_counts_no_legend.pdf", height = 5, width = 5)
				
all_cat_boxplot_none <- ggplot(data = categs_all_move, aes(y=count, x = simple_category, fill = description)) +
				geom_boxplot() +
				scale_fill_manual(values = movment.colors) +
				theme_bw()

ggsave(all_cat_boxplot_none, file = "movement_categories_simple_TFBS_cats_counts_including_none.pdf")

#do statistics on each of the groups to test whether the dynamic/middle/stable vary

stats_for_movement_cats <- function(category_data, name){
	sink(paste(name, "_statistics.txt", sep = ""))
	categories <- unique(category_data$simple_category)
	for (i in seq(1, length(categories))){
		categ <- categories[i]
		test_data <- category_data[category_data$simple_category == categ, ]
		print(name)
		print(categ)
		print("medians")
		test_data_medians <- aggregate(count ~ description, data = test_data, FUN = "median", na.omit = TRUE)
		print(test_data_medians)
		test_data$description <- as.factor(test_data$description)
		print(kruskal.test(count ~ description, data = test_data))
		print(pairwise.wilcox.test(x = test_data$count, g = test_data$description, p.adjust.method = "BH"))
	}
	sink()
}

	
stats_for_movement_cats(category_data = categs_all_move, name = "movement_categories_simple_TFBS_cats_counts_including_none")
stats_for_movement_cats(category_data = categs_move_norm, name = "movement_categories_simple_TFBS_cats_counts_NOT_including_none")

#make the normalised table suitable for the stats function

norm_for_stats <- categs_move_norm[,c("description", "simple_category", "count_norm")]
colnames(norm_for_stats)[3] <- "count"

stats_for_movement_cats(category_data = norm_for_stats, name = "movement_categories_simple_TFBS_cats_normalised_NOT_including_none")
