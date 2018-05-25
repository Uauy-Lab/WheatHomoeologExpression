# aim is to run WGCNA on leaf samples from RNA-seq
# 31-08-2016
# Philippa Borrill

#Using the tutorials at https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/ as  guide

#### TUTORIAL 2c ########
# https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-blockwise.pdf 
# Automatic network construction and module detection with block-wise network construction
# My data set had >10,000 probes which is approx the max that an 8Gb memory computer can handle
# Therefore I need to split hte network into blocks to make it possible to compute
# I guess I could run it on the cluster if I wanted

# setup R correctly
# if in cluster
setwd("/nbi/Research-Groups/NBI/Cristobal-Uauy/expression_browser/WGA/WGCNA/tissues/leaf/")

library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads(3)

##########network analysis itself ##########
#load data from the 1st part of the analysis
lnames=load(file="filtered_data_ready_for_WGCNA_0.5tpm.RData")
# lnames contains the names of loaded variables
lnames
#check datExpr0 looks ok:
datExpr0[1:4,1:4]

# rename datExpr0 to match the tutorial
datExpr <- datExpr0
rm(datExpr0)
datExpr[1:4,1:4]
print("dimensions of datExpr after removing bad genes")
print(dim(datExpr))

# change into new folder to save results
setwd("/nbi/Research-Groups/NBI/Cristobal-Uauy/expression_browser/WGA/WGCNA/tissues/leaf/maxP0.05")


# This shows that for a signed hybrid network I should use a power of 12. 
power <- 12


##Run the model
# try maxBlockSize 46,000 rather than 10,000 i.e. 3 sets of probes together
bwnet = blockwiseModules(datExpr, maxBlockSize = 46000, 
                         power = power, networkType = "signed hybrid", TOMType = "unsigned", minModuleSize = 30,
                         corType="bicor", corOptions = "use = 'p', maxPOutliers = 0.05", 
                         reassignThreshold = 0, mergeCutHeight = 0.15,
                         numericLabels = TRUE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "Signed_hybrid_TOM-blockwise_maxP0.05",
                         verbose = 3)

# Examine module membership

# First let's see what the bwnet dataset contains
names(bwnet)
length(table(bwnet$blocks))
table(bwnet$blocks)

# get the modules colours
bwnetModuleColors <- labels2colors(bwnet$colors)

#get the module labels
bwnetModuleLabels <- bwnet$colors

# look at how many genes per module
table(bwnetModuleLabels)
write.csv(table(bwnetModuleLabels), file="bwnet_modules_mergeCutHeight_0.15.csv")

# plot histogram of module sizes
pdf(file="Number of genes in module mergeCutHeight_0.15.pdf")
par(mfrow=c(1,2))
hist(table(bwnetModuleLabels), xlim=c(1,10000), breaks=1000, xlab="Number of genes in Module", main="mergeCutHeight_0.15")
hist(table(bwnetModuleLabels), xlim=c(1,1000), breaks=1000, xlab="Number of genes in Module", main="Zoomed in mergeCutHeight_0.15")

dev.off()

# plot dendrograms of module membership

# as 1 pdf (1 page per block)
pdf(file="dendrogram_of_module_membership_0.15_mergeCutHeight.pdf", width=6, height=4)

for (i in 1:length(table(bwnet$blocks)))  {        # there are 5 blocks
plotDendroAndColors(bwnet$dendrograms[[i]], bwnetModuleColors[bwnet$blockGenes[[i]]],
                    "Module colors", main = paste0("Gene dendrogram and module colors in block ",i),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
}
dev.off()


# let's save this data
save(bwnet, file="bwnet_network_mergeCutHeight0.15.RData")


# also save the components separately (as the manual suggests)
bwnetMEs <- bwnet$MEs
bwnetdendrograms <- bwnet$dendrograms

save(bwnetModuleColors, bwnetModuleLabels, bwnetMEs, bwnetdendrograms, file = "bwnet_network_components_mergeCutHeight0.15.RData")

### now want to output genes with module membership
gene_expr <- cbind(t(datExpr),bwnetModuleLabels, bwnetModuleColors)
head(gene_expr)

# write genes with modules to file
write.csv(gene_expr,"genes_with_modules_mergeCutHeight0.15.csv")

######playing with clustering settings ########
# check we can load the data
setwd("/nbi/Research-Groups/NBI/Cristobal-Uauy/expression_browser/WGA/WGCNA/tissues/leaf/maxP0.05")

bwnet_test <- load(file="bwnet_network_mergeCutHeight0.15.RData")
bwnet_test


bwnet2 <- recutBlockwiseTrees(datExpr, goodSamples=bwnet$goodSamples, goodGenes =bwnet$goodGenes,
                    blocks = bwnet$blocks, TOMFiles = bwnet$TOMFiles, dendrograms = bwnet$dendrograms,
                    corType = "bicor", corOptions = "use = 'p', maxPOutliers = 0.05", networkType = "signed hybrid",
                    minModuleSize = 30, reassignThreshold = 0, 
                    mergeCutHeight = 0.10,  
                    detectCutHeight = 0.995, # changed from the default of 0.995
                    deepSplit = 2 , # changed from default of 2
                    numericLabels = TRUE,
                    verbose = 3)
                    
# First let's see what the bwnet dataset contains
names(bwnet2)

# get the modules colours
bwnet2ModuleColors <- labels2colors(bwnet2$colors)

#get the module labels
bwnet2ModuleLabels <- bwnet2$colors

# look at how many genes per module
table(bwnet2ModuleLabels)
write.csv(table(bwnet2ModuleLabels), file="bwnet2_modules_mergeCutHeight_0.10.csv")

# plot histogram of module sizes
pdf(file="Number of genes in module mergeCutHeight_0.10.pdf")
par(mfrow=c(1,2))
hist(table(bwnet2ModuleLabels), xlim=c(1,10000), breaks=1000, xlab="Number of genes in Module", main="mergeCutHeight_0.10")
hist(table(bwnet2ModuleLabels), xlim=c(1,1000), breaks=1000, xlab="Number of genes in Module", main="Zoomed in mergeCutHeight_0.10")

dev.off()

# plot dendrograms of module membership

# as 1 pdf (10 pages)
pdf(file="dendrogram_of_module_membership_0.10.pdf", width=6, height=4)
par(mfrow = c(5,2))
for (i in 1:length(table(bwnet$blocks)))  {        # there are 10 blocks
  plotDendroAndColors(bwnet$dendrograms[[i]], bwnet2ModuleColors[bwnet$blockGenes[[i]]],
                      "Module colors", main = paste0("Gene dendrogram and module colors in block ",i),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
}
dev.off()


# now let's save bwnet2 #NB that doesn't include the dendrogram- need that from bwnet original file
save(bwnet2, file="bwnet2_network_mergeCutHeight0.10.RData")
names(bwnet2)

# also save the components separately (as the manual suggests)
bwnetMEs <- bwnet2$MEs
bwnetdendrograms <- bwnet$dendrograms
# get the modules colours
bwnetModuleColors <- labels2colors(bwnet2$colors)

#get the module labels
bwnetModuleLabels <- bwnet2$colors

# 
save(bwnetModuleColors, bwnetModuleLabels, bwnetMEs, bwnetdendrograms, file = 
       "bwnet2_network_components_mergeCutHeight0.10.RData")


### now want to output genes with module membership
gene_expr <- cbind(t(datExpr),bwnetModuleLabels, bwnetModuleColors)
head(gene_expr)

# write genes with modules to file
write.csv(gene_expr,"genes_with_modules_mergeCutHeight0.10.csv")



