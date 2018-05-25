#want to determine whether the number and type of TF motifs are different in the different triad categories
#first read in all the triad sets

all_111_triads <- "W:/Jemima/companion_paper/triads_1-1-1_HC_only.tsv"
all_111_triads_stacked <- "W:/Jemima/companion_paper/triads_1-1-1_HC_only_stacked.tsv"
fimo_unique_path <- "W:/Jemima/companion_paper/out_2018_03_12_12_03/fimo_results/IWGSCv1.0_UTR.HC.firstCDS_consensus._1500bp_upstream_triads_only.fa/unique_hits_analysis/fimo_unique_TF_motifs_per_gene.txt"

#then read in the promoter motif calling 
#start with 1.3 kb upstream, but will do 2kb upstream eventually as well

fimo_unique <- read.table(fimo_unique_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

triads <- read.table(all_111_triads, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

triads_stacked <- read.table(all_111_triads_stacked, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#extract only the fimo results in triads 
fimo_triads <- subset(fimo_unique, gene %in% unique(triads_stacked$gene))
#fimo_triads_ids <- merge(triads_stacked, fimo_triads, all.x = TRUE, all.y = FALSE)

#Now we need to do the breakdown of motif per triads

triad_ids <- unique(triads_stacked$group_id)
motif_ids <- unique(fimo_triads$motif_id)

motif_triad_categ <- data.frame(group_id = character(), motif_id = character(), motif_cat = character())

for (i in seq(1, length(triad_ids))){
	triad <- triad_ids[i]
	
	group_id <- rep(triad, times = length(motif_ids))
	motif_id <- character()
	motif_cat <- character()
		
	triad_genes <- triads_stacked[triads_stacked$group_id == triad,]
	#extract all the motifs
	motifs_triad <- subset(fimo_triads, gene %in% triad_genes$gene)
	motifs_triad_id <- merge(motifs_triad, triad_genes)
	for (j in seq(1, length(motif_ids))){
		motif <- motif_ids[j]
		triad_motif <- motifs_triad_id[motifs_triad_id$motif_id == motif,]
		motif_no <- nrow(triad_motif)
		if (motif_no == 0){
			motif_categ <- "none"
		} else if (motif_no == 1){
			motif_categ <- triad_motif[1,"genome"]
		} else {
			motif_categ <- paste(triad_motif$genome, collapse = "")
		}
		motif_id[j] <- motif
		motif_cat[j] <- motif_categ
	}
	triad_motif_cats <- data.frame(group_id, motif_id, motif_cat)
	motif_triad_categ <- rbind(motif_triad_categ, triad_motif_cats)
}

#correct for things that are the wrong way around
motif_triad_categ$motif_cat <- gsub("BA", "AB", motif_triad_categ$motif_cat)
motif_triad_categ$motif_cat <- gsub("DA", "AD", motif_triad_categ$motif_cat)
motif_triad_categ$motif_cat <- gsub("BDA", "ABD", motif_triad_categ$motif_cat)
motif_triad_categ$motif_cat <- gsub("ADB", "ABD", motif_triad_categ$motif_cat)
motif_triad_categ$motif_cat <- gsub("DBA", "ABD", motif_triad_categ$motif_cat)
motif_triad_categ$motif_cat <- gsub("DB", "BD", motif_triad_categ$motif_cat)
motif_triad_categ$motif_cat <- gsub("BAD", "ABD", motif_triad_categ$motif_cat)



motif_triad_categ$simple_category <- motif_triad_categ$motif_cat

motif_triad_categ$simple_category <- gsub("\\<ABD\\>", "all_same", motif_triad_categ$simple_category)
motif_triad_categ$simple_category <- gsub("\\<none\\>", "all_same", motif_triad_categ$simple_category)
motif_triad_categ$simple_category <- gsub("\\<BD\\>", "A_diff", motif_triad_categ$simple_category)
motif_triad_categ$simple_category <- gsub("\\<AD\\>", "B_diff", motif_triad_categ$simple_category)
motif_triad_categ$simple_category <- gsub("\\<AB\\>", "D_diff", motif_triad_categ$simple_category)
motif_triad_categ$simple_category <- gsub("\\<D\\>", "D_diff", motif_triad_categ$simple_category)
motif_triad_categ$simple_category <- gsub("\\<B\\>", "B_diff", motif_triad_categ$simple_category)
motif_triad_categ$simple_category <- gsub("\\<A\\>", "A_diff", motif_triad_categ$simple_category)

motif_triad_categ_dedup <- unique(motif_triad_categ)
write.table(motif_triad_categ_dedup, file = "W:/Jemima/companion_paper/TF_motifs/triad_categories_per_motif_stacked_1500bp.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


summarised_categories <- data.frame(table(motif_triad_categ$group_id, motif_triad_categ$motif_cat))
colnames(summarised_categories) <- c("group_id", "motif_cat", "count")


summarised_simple_categories <- data.frame(table(motif_triad_categ_dedup$group_id, motif_triad_categ_dedup$simple_category))
colnames(summarised_simple_categories) <- c("group_id", "simple_category", "count")

write.table(summarised_categories, file = "W:/Jemima/companion_paper/TF_motifs/triad_categories_per_triad_stacked_1500bp.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(summarised_simple_categories, file = "W:/Jemima/companion_paper/TF_motifs/triad_categories_per_triad_stacked_simple_categories_1500bp.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

#create a simple summary with no none values included

motif_triad_categ_dedup_no_none <- motif_triad_categ_dedup[!motif_triad_categ_dedup$motif_cat == "none",]
summarised_simple_categories_no_none <- data.frame(table(motif_triad_categ_dedup_no_none$group_id, motif_triad_categ_dedup_no_none$simple_category))
colnames(summarised_simple_categories_no_none) <- c("group_id", "simple_category", "count")
write.table(summarised_simple_categories_no_none, file = "W:/Jemima/companion_paper/TF_motifs/triad_categories_per_triad_stacked_simple_categories_no_none_1500bp.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
