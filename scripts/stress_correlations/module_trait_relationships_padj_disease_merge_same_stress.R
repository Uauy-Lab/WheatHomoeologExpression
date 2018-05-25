# calculate module_trait relationships
# Philippa Borrill
# 28.7.2017

# Calculate module association to tissue, age, disease etc.

# 28.7.2017

setwd("Y:\\expression_browser\\WGA\\WGCNA\\abiotic_vs_disease\\disease\\maxP0.05\\")

eigengenes <- read.table(file="eigengenes.txt", header=T)
head(eigengenes)

# remove eigengene 0
eigengenes <- eigengenes[,2:ncol(eigengenes)]

# read meta-data
metadata_file <- "Y:\\expression_browser\\WGA\\WGCNA\\abiotic_vs_disease\\disease\\WGCNA_metadata_cluster_intermediate_stress.csv"
metadata <- read.csv(metadata_file, row.names = 1, sep=",")

dim(metadata)
head(metadata)
colnames(metadata)

metadata$Intermediate <- gsub("stripe rust control","mildew stripe rust control", (metadata$Intermediate))
metadata$Intermediate <- gsub("stripe rust 2","stripe rust", (metadata$Intermediate))
metadata$Intermediate <- gsub(" study2","",(metadata$Intermediate))
unique(metadata$Intermediate)

metadata_Sel <- metadata[,c(3,5:13,21)]
head(metadata_Sel)
metadata_Sel$Sample.IDs <- gsub("-",".",metadata_Sel$Sample.IDs)

# check rows are in same order in eigengenes and metadata
metadata_Sel$Sample.IDs == rownames(eigengenes)

colnames(metadata_Sel)

######### intermediate level stress #########

# now want to convert metadata to be binary

colnames(metadata_Sel)

metadata_intermed_stress <- metadata_Sel[,c(1,11)]
head(metadata_intermed_stress)
unique(metadata_intermed_stress$Intermediate)

# make new column for each intermed stress
metadata_intermed_stress$Fusarium_pseudograminearum <-  metadata_intermed_stress$Intermediate == "Fusarium pseudograminearum"
metadata_intermed_stress$Zymoseptoria_tritici <-  metadata_intermed_stress$Intermediate == "Zymoseptoria tritici"
metadata_intermed_stress$powdery_mildew <-  metadata_intermed_stress$Intermediate == "powdery mildew"
metadata_intermed_stress$stripe_rust <-  metadata_intermed_stress$Intermediate == "stripe rust"
metadata_intermed_stress$PAMP_chitin <-  metadata_intermed_stress$Intermediate == "PAMP chitin"
metadata_intermed_stress$PAMP_flg22 <-  metadata_intermed_stress$Intermediate == "PAMP flg22"
metadata_intermed_stress$Fusarium_pseudograminearum_control <-  metadata_intermed_stress$Intermediate == "Fusarium pseudograminearum control"
metadata_intermed_stress$Zymoseptoria_tritici_control <-  metadata_intermed_stress$Intermediate == "Zymoseptoria tritici control"
metadata_intermed_stress$mildew_stripe_rust2_control <-  metadata_intermed_stress$Intermediate == "mildew stripe rust control"
metadata_intermed_stress$PAMP_control <-  metadata_intermed_stress$Intermediate == "PAMP control"

head(metadata_intermed_stress,20)
tail(metadata_intermed_stress,20)

cols <- sapply(metadata_intermed_stress, is.logical)
metadata_intermed_stress[,cols] <- lapply(metadata_intermed_stress[,cols], as.numeric)
head(metadata_intermed_stress)
rownames(metadata_intermed_stress) <- metadata_intermed_stress$Sample.IDs

dim(metadata_intermed_stress)
metadata_logical <- metadata_intermed_stress[,3:ncol(metadata_intermed_stress)]
head(metadata_logical)
rownames(metadata_logical) ==rownames(eigengenes)
names(metadata_logical)

nSamples <- nrow(metadata_logical)
# now want to do a comparison between eigengenes and each column in metadata (i.e. intermed level stress)
library(WGCNA)
moduleTraitCor = cor(eigengenes, metadata_logical, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
head(moduleTraitPvalue)

# FDR adjust pvalues
moduleTraitPvalue <- matrix(p.adjust(as.vector(as.matrix(moduleTraitPvalue)), method="BY"),ncol=ncol(moduleTraitPvalue))
head(moduleTraitPvalue)
rownames(moduleTraitPvalue) <- rownames(moduleTraitCor)
colnames(moduleTraitPvalue) <- colnames(moduleTraitCor)
head(moduleTraitPvalue)

# now display eigengenes to metadata relationships
sizeGrWindow(15,6)


# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(8, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(metadata_logical),
               yLabels = names(eigengenes),
               ySymbols = names(eigengenes),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Intermediate level stress module-trait relationships"))



setwd("Y:\\expression_browser\\WGA\\WGCNA\\abiotic_vs_disease\\disease\\maxP0.05\\trait_module_correlations_merge_same_stress")

# now save as pdf and text matrix

write.csv(moduleTraitCor,file="Intermediate_level_stress_corr.csv")
write.csv(moduleTraitPvalue,file="Intermediate_level_stress_p-value.csv")

pdf(file="Intermediate-level-stress_corr_p-value_MEs.pdf", height=15, width=10)
par(mar = c(14, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(metadata_logical),
               yLabels = names(eigengenes),
               cex.lab.y = 0.5,
               ySymbols = names(eigengenes),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.25,
               zlim = c(-1,1),
               main = paste("Intermediate level stress module-trait relationships"))
dev.off()

# now make pdf colour by significance
library("RColorBrewer")

hmcols<-colorRampPalette(c("red","white"))(256)


pdf(file="Intermediate-level-stress_corr_p-value_MEs_coloured_by_padj_0.05.pdf", height=15, width=10)
par(mar = c(14, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitPvalue,
               xLabels = names(metadata_logical),
               yLabels = names(eigengenes),
               cex.lab.y = 0.5,
               ySymbols = names(eigengenes),
               colorLabels = FALSE,
               colors = hmcols,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.25,
               zlim = c(0,0.05),
               main = paste("Intermediate level stress module-trait relationships"))
dev.off()

pdf(file="Intermediate-level-stress_corr_p-value_MEs_coloured_by_padj_0.01.pdf", height=15, width=10)
par(mar = c(14, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitPvalue,
               xLabels = names(metadata_logical),
               yLabels = names(eigengenes),
               cex.lab.y = 0.5,
               ySymbols = names(eigengenes),
               colorLabels = FALSE,
               colors = hmcols,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.25,
               zlim = c(0,0.01),
               main = paste("Intermediate level stress module-trait relationships"))
dev.off()

pdf(file="Intermediate-level-stress_corr_p-value_MEs_coloured_by_padj_0.001.pdf", height=15, width=10)
par(mar = c(14, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitPvalue,
               xLabels = names(metadata_logical),
               yLabels = names(eigengenes),
               cex.lab.y = 0.5,
               ySymbols = names(eigengenes),
               colorLabels = FALSE,
               colors = hmcols,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.25,
               zlim = c(0,0.001),
               main = paste("Intermediate level stress module-trait relationships"))
dev.off()


