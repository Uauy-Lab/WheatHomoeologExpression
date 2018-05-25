# calculate hub genes and export for visualisation 

# 25.5.2017 # 25.8.2017

# check we can load the data
setwd("Y:\\expression_browser\\WGA\\WGCNA\\abiotic_vs_disease\\disease\\maxP0.05\\top_kME\\")

library(WGCNA)
options(stringsAsFactors = FALSE)

#load data from the 1st part of the analysis
lnames=load(file="Y:\\expression_browser\\WGA\\WGCNA\\abiotic_vs_disease\\disease\\filtered_data_ready_for_WGCNA_0.5tpm.Rdata")
# lnames contains the names of loaded variables
lnames

datExpr <- datExpr0
rm(datExpr0)

# Load network data saved in the second part.
# mergecutheight0.15

lnames=load(file="Y:\\expression_browser\\WGA\\WGCNA\\abiotic_vs_disease\\disease\\maxP0.05\\bwnet_network_components_mergeCutHeight0.15.Rdata")
lnames

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

## calculate hub genes (most highly connected = high kME values)
# read in module info
module_info <- read.csv(file="Y:\\expression_browser\\WGA\\WGCNA\\abiotic_vs_disease\\disease\\maxP0.05\\genes_with_modules_mergeCutHeight0.15_no_expr_values.csv", 
                        header=T, sep=" ")
head(module_info)
dim(module_info)

modulesA1 <- as.vector(module_info$bwnetModuleColors)
PCs <- moduleEigengenes((datExpr), colors=modulesA1)
ME <- PCs$eigengenes
colorsA1 <- names(table(modulesA1))
colorsA1

# now get KME values
geneModuleMembership <- signedKME(datExpr,ME)
colnames(geneModuleMembership) <- paste0("PC", colorsA1,".cor")
head(geneModuleMembership)

MMPvalue1 <- corPvalueStudent(as.matrix(geneModuleMembership), dim(datExpr)[[2]])
colnames(MMPvalue1) <- paste0("PC",colorsA1,".pval")
head(MMPvalue1)

Gene <- colnames(datExpr)
head(Gene)
kMEtable1 <- cbind(Gene, Gene, modulesA1)
head(kMEtable1)

for (i in 1: length(colorsA1)){
  kMEtable1 <- cbind(kMEtable1, geneModuleMembership[,i], MMPvalue1[,i] )
}

colnames(kMEtable1) <- c("PSID","Gene", "Module", sort(c(colnames(geneModuleMembership),colnames(MMPvalue1))))
head(kMEtable1)
dim(kMEtable1)
write.csv(kMEtable1, "kMEtable1.csv", row.names=F)

#find top ten hub genes per module
head(geneModuleMembership)

topGenesKME <- NULL

for (c in 1:length(colorsA1)) {
  kMErank1 <- rank(-geneModuleMembership[,c])
  topGenesKME <- cbind(topGenesKME,Gene[kMErank1<=10])
}


head(sort(geneModuleMembership$PCblack.cor, decreasing=T))

colnames(topGenesKME) <- colorsA1
topGenesKME

library("reshape2")
melted_topGenesKME <- melt(topGenesKME)
head(melted_topGenesKME)
melted_topGenesKME



write.csv(melted_topGenesKME_TF_modNo, file="top10_genes_kME_each_module.csv")

