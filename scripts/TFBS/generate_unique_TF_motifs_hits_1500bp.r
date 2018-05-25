library(stringr)

#read in the raw data

fimo_raw_dir <- "W:/Jemima/companion_paper/out_2018_03_12_12_03/fimo_results/IWGSCv1.0_UTR.HC.firstCDS_consensus._1500bp_upstream_triads_only.fa/unique_hits_analysis/"

fimo_raw <- paste(fimo_raw_dir, "fimo_unique.tsv", sep = "")

fimo <- read.table(fimo_raw, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
print("table_read_in")
#add column names
colnames(fimo) <- c("motif_id", "motif_alt_id", "sequence_name")

#extract the gene_id from the transcript id
fimo$gene <- str_split_fixed(fimo$sequence_name, "\\.", 3)[,1]

#get dataframe with just the gene_id and the motif_id
fimo_simple <- fimo[,c("gene", "motif_id")]

fimo_unique <- unique(fimo_simple)

#count the number of unique motifs per gene
fimo_unique_counts <- data.frame(table(fimo_unique$gene))

colnames(fimo_unique_counts) <- c("gene", "unique_motifs")

#write the unique table
out_prefix <- gsub(".tsv", "", fimo_raw)
print("out_prefix")
write.table(fimo_unique, file = paste(out_prefix, "_TF_motifs_per_gene.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

#write the unique counts table
write.table(fimo_unique_counts, file = paste(out_prefix, "_TF_motif_counts_per_gene.txt", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
 print("finished")