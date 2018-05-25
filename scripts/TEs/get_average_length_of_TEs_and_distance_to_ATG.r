#### script includes all calculations. Final statistics used were including ALL triads in the the HC_no_stress categories

#read in the data
library(ggplot2)

group.colors<-c(A.dominant = "#579D1C", B.dominant = "#4B1F6F", D.dominant ="#FF950E",
             Central="#AAAAAA",Balanced="#AAAAAA",
             A.suppressed = "#b2E08a", B.suppressed = "#0eC7ff", D.suppressed ="#ffCF0e")
chr.colors<-c(A = "#579D1C", B = "#4B1F6F", D ="#FF950E")


movment.colors = c("Dynamic"="#e41a1c", "Middle 80"="#bdbdbd", "Stable"="#377eb8")

five.cat.colors <- c("non_suppressed" = "#a6611a", "suppressed" = "#dfc27d", "central" = "#AAAAAA", "non_dominant" = "#80cdc1", "dominant" = "#018571")


TE_111_all <- read.table("W:/Jemima/companion_paper/TEs/all_TE_position_size_relative_positions_ATG_start.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

triads_no_N <- read.table("W:/Jemima/companion_paper/genes_with_no_Ns_in_any_of_triad_4500bp_upstream.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#low 10 pc
HC_CS_low10 <- read.table("W:/expression_browser/dataset_reports/20170911_all_reports/02.movement/HC_CS_no_stress/HC_CS_no_stress_movement_low_10pc.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
HC_CS_top10 <- read.table("W:/expression_browser/dataset_reports/20170911_all_reports/02.movement/HC_CS_no_stress/HC_CS_no_stress_movement_top_10pc.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
HC_CS_middle80 <- read.table("W:/expression_browser/dataset_reports/20170911_all_reports/02.movement/HC_CS_no_stress/HC_CS_no_stress_movement_middle_80pc.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)


triads_no_N$description <- NA

triads_no_N[which(triads_no_N$gene_ids %in% HC_CS_low10$x),"description"] <- "HC_CS_no_stress_low10pc"
triads_no_N[which(triads_no_N$gene_ids %in% HC_CS_top10$x),"description"] <- "HC_CS_no_stress_top10pc"
triads_no_N[which(triads_no_N$gene_ids %in% HC_CS_middle80$x),"description"] <- "HC_CS_no_stress_middle80pc"

triads_move <- triads_no_N[complete.cases(triads_no_N),]

triads_TE <- subset(TE_111_all, geneId %in% triads_move$gene_ids)
colnames(triads_TE)[2] <- "gene"
colnames(triads_move)[1] <- "gene"

triads_move_TE <- merge(triads_move, triads_TE, all.x = TRUE, all.y = TRUE)
triads_move_TE_noNA <- triads_move_TE[complete.cases(triads_move_TE),]

triads_move_TE_means <- aggregate(TE_size ~ description, data = triads_move_TE, FUN = "mean")
triads_move_TE_medians <- aggregate(TE_size ~ description, data = triads_move_TE, FUN = "median")

#do the statistics now - overall kruskal test and then the wilcox pairwise test
triads_move_TE_noNA$description <- as.factor(triads_move_TE_noNA$description)
kruskal_HC_no_stress <- kruskal.test(TE_size ~ description, data = triads_move_TE_noNA)
#wilcox with benjamni hochbery adjust P-value

wilcox_no_stress <- pairwise.wilcox.test(x = triads_move_TE_noNA$TE_size, g = triads_move_TE_noNA$description, p.adjust.method = "BH")
movment.colors = c("HC_CS_no_stress_top10pc"="#e41a1c", "HC_CS_no_stress_middle80pc"="#bdbdbd", "HC_CS_no_stress_low10pc"="#377eb8")

size_plot <- ggplot(triads_move_TE, aes(x= TE_size, fill = description)) +
			geom_density(alpha = 0.5) +
			scale_fill_manual(values = movment.colors) + theme_bw() 
			xlim(0, 30000) +
			theme_bw()

ggsave(size_plot, file = "TE_size_distribution_movement.pdf")

size_boxplot <- ggplot(triads_move_TE, aes(x = description, y= TE_size, fill = description)) +
			geom_boxplot() +
			scale_fill_manual(values = movment.colors) +
			coord_cartesian(ylim = c(0, 10000))+
			theme_bw()

triads_move_TE$TE_mid_rel <- (triads_move_TE$TE_start_rel + triads_move_TE$TE_end_rel)/2

size_pos <- ggplot(triads_move_TE, aes(x = TE_mid_rel, y= TE_size, fill = description, colour = description)) +
			geom_point() +
			scale_fill_manual(values = movment.colors) +
			scale_colour_manual(values = movment.colors) +
			theme_bw()
			
size_start <- ggplot(triads_move_TE, aes(x = TE_start_rel, y= TE_size, fill = description, colour = description)) +
			geom_point() +
			scale_fill_manual(values = movment.colors) +
			scale_colour_manual(values = movment.colors) +
			theme_bw()
			
size_end <- ggplot(triads_move_TE[!triads_move_TE$description == "HC_CS_no_stress_middle80pc",], aes(x = TE_end_rel, y= TE_size, fill = description, colour = description)) +
			geom_point() +
			scale_fill_manual(values = movment.colors) +
			scale_colour_manual(values = movment.colors) +
			theme_bw()
			
ggsave(size_end, file = "TE_size_vs_TE_relative_end_position.pdf", width = 10, height = 7)

#### Ok, so we want a table with the median and standard deviation of TE size for: all triads, triads without Ns and at 5kb and 1.5 kb from the chinese_spring_no_stress
#triads_move_TE is currently just the 5kb no NA

colnames(triads_move_TE_medians)[2] <- "median_5kb_noN"
#get the 5kb noN stdev
triads_move_TE_medians$st_dev_5kb_noN <- aggregate(TE_size ~ description, data = triads_move_TE, FUN = "sd", na.action = na.omit)[,2]


triads_move_TE$closest_TE_point <- 0
triads_move_TE[(triads_move_TE$TE_start_rel > triads_move_TE$TE_end_rel), "closest_TE_point"] <- triads_move_TE$TE_start_rel

#get the median for just 1.5kb
triads_move_TE_1500 <- triads_move_TE[((triads_move_TE$TE_start_rel >= -1500) | (triads_move_TE$TE_end_rel >= -1500)),]

triads_move_TE_1500 <- triads_move_TE_1500[complete.cases(triads_move_TE_1500),]

triads_move_TE_medians$median_1.5kb_noN <- aggregate(TE_size ~ description, data = triads_move_TE_1500, FUN = "median", na.action = na.omit)[,2]
triads_move_TE_medians$st_dev_1.5kb_noN <- aggregate(TE_size ~ description, data = triads_move_TE_1500, FUN = "sd", na.action = na.omit)[,2]

#do the statistics for the 1.5kb - overall kruskal test and then the wilcox pairwise test
triads_move_TE_1500$description <- as.factor(triads_move_TE_1500$description)
kruskal_HC_no_stress_1500 <- kruskal.test(TE_size ~ description, data = triads_move_TE_1500)
#wilcox with benjamni hochbery adjust P-value

wilcox_no_stress_1500 <- pairwise.wilcox.test(x = triads_move_TE_1500$TE_size, g = triads_move_TE_1500$description, p.adjust.method = "BH")

#now perform the same analysis for all the triads
TE_111_all$description <- NA

TE_111_all[which(TE_111_all$geneId %in% HC_CS_low10$x),"description"] <- "HC_CS_no_stress_low10pc"
TE_111_all[which(TE_111_all$geneId %in% HC_CS_top10$x),"description"] <- "HC_CS_no_stress_top10pc"
TE_111_all[which(TE_111_all$geneId %in% HC_CS_middle80$x),"description"] <- "HC_CS_no_stress_middle80pc"

triads_move_TE_medians$median_5kb_all <- aggregate(TE_size ~ description, data = TE_111_all, FUN = "median", na.action = na.omit)[,2]
triads_move_TE_medians$stdev_5kb_all <- aggregate(TE_size ~ description, data = TE_111_all, FUN = "sd", na.action = na.omit)[,2]

## all triads_1500
TE_111_all_1500 <- TE_111_all[((TE_111_all$TE_start_rel >= -1500) | (TE_111_all$TE_end_rel >= -1500)),]

TE_111_all_1500 <- TE_111_all_1500[complete.cases(TE_111_all_1500),]

triads_move_TE_medians$median_1.5kb_all <- aggregate(TE_size ~ description, data = TE_111_all_1500, FUN = "median", na.action = na.omit)[,2]
triads_move_TE_medians$st_dev_1.5kb_all <- aggregate(TE_size ~ description, data = TE_111_all_1500, FUN = "sd", na.action = na.omit)[,2]

## statistics for all promoters
#do the statistics for the 1.5kb - overall kruskal test and then the wilcox pairwise test
TE_111_all$description <- as.factor(TE_111_all$description)
kruskal_HC_no_stress_all <- kruskal.test(TE_size ~ description, data = TE_111_all)
#wilcox with benjamni hochbery adjust P-value

wilcox_no_stress_all <- pairwise.wilcox.test(x = TE_111_all$TE_size, g = TE_111_all$description, p.adjust.method = "BH")

#do the statistics for the 1.5kb - overall kruskal test and then the wilcox pairwise test
TE_111_all_1500$description <- as.factor(TE_111_all_1500$description)
kruskal_HC_no_stress_all_1500 <- kruskal.test(TE_size ~ description, data = TE_111_all_1500)
#wilcox with benjamni hochberg adjust P-value

wilcox_no_stress_all_1500 <- pairwise.wilcox.test(x = TE_111_all_1500$TE_size, g = TE_111_all_1500$description, p.adjust.method = "BH")


write.table(triads_move_TE_medians, file = "TE_size_medians_and_stdev_movement_categories.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)

#### now calculate the sizes for the five dominance categories

triads_five <- read.table(file = "W:/Jemima/companion_paper/HC_CS_no_stress_triads_with_five_cats.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colnames(triads_five)[4] <- "geneId"
###first do all promoters including those with Ns
all_TE_size_five <- subset(TE_111_all, geneId %in% triads_five$geneId)
all_TE_size_five <- merge(all_TE_size_five, triads_five, by = "geneId")

get_medians_stats_five_cat <- function(all_TE_size_five){
	triads_dom_TE_medians <- aggregate(TE_size ~ five_cat, data = all_TE_size_five, FUN = "median", na.action = na.omit)
	colnames(triads_dom_TE_medians)[2] <- "median_5kb_all"
	triads_dom_TE_medians$stdev_5kb_all <- aggregate(TE_size ~ five_cat, data = all_TE_size_five, FUN = "sd", na.action = na.omit)[,2]
	
	##1500bp
	all_TE_size_five_1500 <- all_TE_size_five[((all_TE_size_five$TE_start_rel >= -1500) | (all_TE_size_five$TE_end_rel >= -1500)),]
	all_TE_size_five_1500 <- all_TE_size_five_1500[complete.cases(all_TE_size_five_1500),]
	
	triads_dom_TE_medians$median_1.5kb_all <- aggregate(TE_size ~ five_cat, data = all_TE_size_five_1500, FUN = "median", na.action = na.omit)[,2]
	triads_dom_TE_medians$stdev_1.5kb_all <- aggregate(TE_size ~ five_cat, data = all_TE_size_five_1500, FUN = "sd", na.action = na.omit)[,2]
	
	#stats with all promoters
	all_TE_size_five$five_cat <- as.factor(all_TE_size_five$five_cat)
	kruskal_all_TE_size_five <- kruskal.test(TE_size ~ five_cat, data = all_TE_size_five)
	print(kruskal_all_TE_size_five)
	#wilcox with benjamni hochberg adjust P-value
	wilcox_all_TE_size_five <- pairwise.wilcox.test(x = all_TE_size_five$TE_size, g = all_TE_size_five$five_cat, p.adjust.method = "BH")
	print(wilcox_all_TE_size_five)
	#do the statistics for the 1.5kb - overall kruskal test and then the wilcox pairwise test
	all_TE_size_five_1500$five_cat <- as.factor(all_TE_size_five_1500$five_cat)
	kruskal_all_TE_size_five_1500 <- kruskal.test(TE_size ~ five_cat, data = all_TE_size_five_1500)
	print(kruskal_all_TE_size_five_1500)
	#wilcox with benjamni hochberg adjust P-value
	wilcox_all_TE_size_five_1500 <- pairwise.wilcox.test(x = all_TE_size_five_1500$TE_size, g = all_TE_size_five_1500$five_cat, p.adjust.method = "BH")
	print(wilcox_all_TE_size_five_1500)
	return(triads_dom_TE_medians)
}

triads_dom_TE_medians2 <- get_medians_stats_five_cat(all_TE_size_five)

#do the same but without the Ns promoters
all_TE_size_five_noN <- subset(all_TE_size_five, geneId %in% triads_no_N$gene_ids)
triads_dom_TE_medians_noN <- get_medians_stats_five_cat(all_TE_size_five_noN)

colnames(triads_dom_TE_medians_noN) <- gsub("_all", "_noN", colnames(triads_dom_TE_medians_noN))

triads_dom_TE_medians_all <- merge(triads_dom_TE_medians2, triads_dom_TE_medians_noN)
write.table(triads_dom_TE_medians_all, file = "TE_size_medians_and_stdev_dom-supp_categories_HC_no_stress.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)


# now for both get distributions of the closest TE
# for the TE_111_all table, for each gene, find the closest start and end
gene_list <- unique(TE_111_all$geneId)
closest_TE <- numeric()

for (i in seq(1, length(gene_list))){
	gene <- gene_list[i]
	gene_TE <- TE_111_all[TE_111_all$geneId == gene, ]
	closest_TE[i] <- max(c(gene_TE$TE_start_rel, gene_TE$TE_end_rel))

}

closest_TE_df <- data.frame(gene_list, closest_TE)

#now plot the distributions of these for the 10-80-10 and also the five-categories
closest_TE_df$movement_cat <- NA
colnames(closest_TE_df)[1] <- "geneId"


closest_TE_df[which(closest_TE_df$geneId %in% HC_CS_low10$x),"movement_cat"] <- "Stable"
closest_TE_df[which(closest_TE_df$geneId %in% HC_CS_top10$x),"movement_cat"] <- "Dynamic"
closest_TE_df[which(closest_TE_df$geneId %in% HC_CS_middle80$x),"movement_cat"] <- "Middle 80"


closest_TE_df <- merge(closest_TE_df, triads_five, by = "geneId", all.x = TRUE, all.y = FALSE)
write.table(closest_TE_df, sep = ",", col.names = TRUE, row.names = FALSE, quote = 
movement_closest_TE <- ggplot(closest_TE_df[complete.cases(closest_TE_df$movement_cat),], aes(x=closest_TE, fill = movement_cat)) +
						geom_density(alpha = 0.5) +
						scale_fill_manual(values = movment.colors) +
						theme_bw()
						
movement_closest_TE_box <- ggplot(closest_TE_df[complete.cases(closest_TE_df$movement_cat),], aes(x=movement_cat, y = closest_TE, fill = movement_cat)) +
						geom_boxplot() +
						scale_fill_manual(values = movment.colors) +
						theme_bw()
						
fivecat_closest_TE <- ggplot(closest_TE_df[complete.cases(closest_TE_df$five_cat),], aes(x=closest_TE, fill = five_cat)) +
						geom_density(alpha = 0.5) +
						scale_fill_manual(values = five.cat.colors) +
						theme_bw()

fivecat_closest_TE_box <- ggplot(closest_TE_df[complete.cases(closest_TE_df$five_cat),], aes(x=five_cat, y=closest_TE, fill = five_cat)) +
						geom_boxplot() +
						scale_fill_manual(values = five.cat.colors) +
						theme_bw()
						
ggsave(movement_closest_TE, file = "closest_TE_density_movement_cats.pdf")
ggsave(movement_closest_TE_box, file = "closest_TE_boxplot_movement_cats.pdf")
ggsave(fivecat_closest_TE, file = "closest_TE_density_five_dom-supp_cats.pdf")
ggsave(fivecat_closest_TE_box, file = "closest_TE_boxplot_five_dom-supp_cats.pdf")

##get medians and statistics
sink("closest_TE_medians_and_stats_movement_and_five_cat.txt")
closest_TE_medians_five_cat <- aggregate(closest_TE ~ five_cat, data = closest_TE_df, FUN = "median", na.action = na.omit)
closest_TE_medians_five_cat$stdev <- aggregate(closest_TE ~ five_cat, data = closest_TE_df, FUN = "sd", na.action = na.omit)[,2]
print("five cats medians and stats\n")
print(closest_TE_medians_five_cat)
##statistics - kruskal followed by wilcox pairwise
closest_TE_five <- closest_TE_df[complete.cases(closest_TE_df$five_cat),]
closest_TE_five$five_cat <- as.factor(closest_TE_five$five_cat)

kruskal_closest_TE_five <- kruskal.test(closest_TE ~ five_cat, data = closest_TE_five)
print(kruskal_closest_TE_five)
#wilcox with benjamni hochberg adjust P-value
wilcox_closest_TE_five <- pairwise.wilcox.test(x = closest_TE_five$closest_TE, g = closest_TE_five$five_cat, p.adjust.method = "BH")
print(wilcox_closest_TE_five)

### calculate the same for the movement categories
##get medians and statistics
closest_TE_medians_movement <- aggregate(closest_TE ~ movement_cat, data = closest_TE_df, FUN = "median", na.action = na.omit)
closest_TE_medians_movement$stdev <- aggregate(closest_TE ~ movement_cat, data = closest_TE_df, FUN = "sd", na.action = na.omit)[,2]
print("movement cats medians and stats\n")
print(closest_TE_medians_movement)
##statistics - kruskal followed by wilcox pairwise
closest_TE_movement <- closest_TE_df[complete.cases(closest_TE_df$movement_cat),]
closest_TE_movement$movement_cat <- as.factor(closest_TE_movement$movement_cat)

kruskal_closest_TE_movement <- kruskal.test(closest_TE ~ movement_cat, data = closest_TE_movement)
print(kruskal_closest_TE_movement)
#wilcox with benjamni hochberg adjust P-value
wilcox_closest_TE_movement <- pairwise.wilcox.test(x = closest_TE_movement$closest_TE, g = closest_TE_movement$movement_cat, p.adjust.method = "BH")
print(wilcox_closest_TE_movement)
sink()


##### repeat with no Ns
closest_TE_noN <- subset(closest_TE_df, geneId %in% triads_no_N$gene_ids)
##get medians and statistics
sink("closest_TE_medians_and_stats_movement_and_five_cat_noN.txt")
closest_TE_medians_five_cat <- aggregate(closest_TE ~ five_cat, data = closest_TE_noN, FUN = "median", na.action = na.omit)
closest_TE_medians_five_cat$stdev <- aggregate(closest_TE ~ five_cat, data = closest_TE_noN, FUN = "sd", na.action = na.omit)[,2]
print("five cats medians and stats\n")
print(closest_TE_medians_five_cat)
##statistics - kruskal followed by wilcox pairwise
closest_TE_five <- closest_TE_noN[complete.cases(closest_TE_noN$five_cat),]
closest_TE_five$five_cat <- as.factor(closest_TE_five$five_cat)

kruskal_closest_TE_five <- kruskal.test(closest_TE ~ five_cat, data = closest_TE_five)
print(kruskal_closest_TE_five)
#wilcox with benjamni hochberg adjust P-value
wilcox_closest_TE_five <- pairwise.wilcox.test(x = closest_TE_five$closest_TE, g = closest_TE_five$five_cat, p.adjust.method = "BH")
print(wilcox_closest_TE_five)

### calculate the same for the movement categories
##get medians and statistics
closest_TE_medians_movement <- aggregate(closest_TE ~ movement_cat, data = closest_TE_noN, FUN = "median", na.action = na.omit)
closest_TE_medians_movement$stdev <- aggregate(closest_TE ~ movement_cat, data = closest_TE_noN, FUN = "sd", na.action = na.omit)[,2]
print("movement cats medians and stats\n")
print(closest_TE_medians_movement)
##statistics - kruskal followed by wilcox pairwise
closest_TE_movement <- closest_TE_noN[complete.cases(closest_TE_noN$movement_cat),]
closest_TE_movement$movement_cat <- as.factor(closest_TE_movement$movement_cat)

kruskal_closest_TE_movement <- kruskal.test(closest_TE ~ movement_cat, data = closest_TE_movement)
print(kruskal_closest_TE_movement)
#wilcox with benjamni hochberg adjust P-value
wilcox_closest_TE_movement <- pairwise.wilcox.test(x = closest_TE_movement$closest_TE, g = closest_TE_movement$movement_cat, p.adjust.method = "BH")
print(wilcox_closest_TE_movement)
sink()