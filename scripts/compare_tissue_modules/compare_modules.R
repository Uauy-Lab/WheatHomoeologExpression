# compare the overlap between modules in different tissues

# Philippa Borrill
# 1.8.2017

# Found a package called "GeneOverlap" which can compare multiple lists 
# and calculate fishers exact test p-values and Jaccard index (0 means lists totally different, 1 lists identical)

# Do this for a test set first

# Compare grain and spike, in the end will need to do 6 comparisons (need to do pairwise)

setwd("Y:\\expression_browser\\WGA\\WGCNA\\tissues\\common_unique_modules\\genes_expr_all_networks")

source("https://bioconductor.org/biocLite.R")
biocLite(pkgs=c("GeneOverlap"))
library(GeneOverlap)

# first make a unified list where I only keep genes expressed in all 4 tissues 
# (NB. will need to go back and look at genes which were unique to one tissue which might be all in 1 module which is tissue specific)
grain_modules <- read.csv("Y:\\expression_browser\\WGA\\WGCNA\\tissues\\grain\\maxP0.05\\genes_with_modules_mergeCutHeight0.15_no_expr_values.csv", sep=" ")
spike_modules <- read.csv("Y:\\expression_browser\\WGA\\WGCNA\\tissues\\spike\\maxP0.05\\genes_with_modules_mergeCutHeight0.15_no_expr_values.csv", sep=" ")
root_modules <- read.csv("Y:\\expression_browser\\WGA\\WGCNA\\tissues\\root\\maxP0.05\\genes_with_modules_mergeCutHeight0.15_no_expr_values.csv", sep=" ")
leaf_modules <- read.csv("Y:\\expression_browser\\WGA\\WGCNA\\tissues\\leaf\\maxP0.05\\genes_with_modules_mergeCutHeight0.15_no_expr_values.csv", sep=" ")

head(grain_modules)
head(spike_modules)
head(root_modules)
head(leaf_modules)

modules <- merge(grain_modules, spike_modules, by ="X", suffixes = c(".grain",".spike"))
head(modules)
modules <- merge(modules, root_modules, by ="X")
head(modules)
modules <- merge(modules, leaf_modules, by ="X", suffixes = c(".root",".leaf"))
head(modules)

modules <- modules[,c(1,3,5,7,9)]
head(modules)
dim(modules)

write.csv(file="modules_for_all_4_tissues.csv", modules)

# only 59,777 genes were expressed in all tissues so can be compared


#### for just two modules and two tissues as a practise #######
grain_1 <- modules[modules$bwnetModuleLabels.grain =="1",][,1]
grain_2 <- modules[modules$bwnetModuleLabels.grain =="2",][,1]

spike_1 <- modules[modules$bwnetModuleLabels.spike =="1",][,1]
spike_2 <- modules[modules$bwnetModuleLabels.spike =="2",][,1]

grain_list <- list(grain_1,grain_2) # put both vectors into a list
names(grain_list) <- c("grain_1","grain_2")
sapply(grain_list,length) # check lengths of list

spike_list <- list(spike_1,spike_2)
names(spike_list) <- c("spike_1","spike_2")
sapply(spike_list,length) # check lengths of list

background_size <- nrow(modules) # need the number of genes in the background
background_size

# now just do for 1 list vs 1 list
go.obj <- newGeneOverlap(spike_list$spike_1, 
                         grain_list$grain_1,
                         genome.size = background_size)
go.obj

go.obj <- testGeneOverlap(go.obj)
go.obj
print(go.obj)

# now do for both lists against both lists
gom.obj <- newGOM(spike_list, grain_list, background_size)
drawHeatmap(gom.obj) # this shows the odds ratio as the colour

pvalue_matrix <- getMatrix(gom.obj, name="pval") # get pvalue as a table
Jaccard_matrix <- getMatrix(gom.obj, name="Jaccard") # get Jaccard as table

# we can also plot showing the Jaccard index as the colour
drawHeatmap(gom.obj, what="Jaccard", ncolused = 5, grid.col = "Blues", note.col = "black")

write.csv(file="pvalues_matrix.csv",pvalue_matrix )
write.csv(file="Jaccard_matrix.csv",Jaccard_matrix )


### now for all tissues and all modules #####

#### first grain vs spike ####
background_size <- nrow(modules) # need the number of genes in the background
background_size

# make grainList which is a list of character vectors containing gene names for each module
grainList <- list()
for (i in 0:(length(unique(modules$bwnetModuleLabels.grain))-1)) {
  print(i)
  grainList[[i+1]] <- assign(paste0("grain",i), modules[modules$bwnetModuleLabels.grain ==i,][,1])
}

names(grainList) <- c(paste0("grain",0:(length(unique(modules$bwnetModuleLabels.grain))-1)))
sapply(grainList,length) # check lengths of list


# make spikeList which is a list of character vectors containing gene names for each module
spikeList <- list()
for (i in 0:(length(unique(modules$bwnetModuleLabels.spike))-1)) {
  print(i)
  spikeList[[i+1]] <- assign(paste0("spike",i), modules[modules$bwnetModuleLabels.spike ==i,][,1])
}

names(spikeList) <- c(paste0("spike",0:(length(unique(modules$bwnetModuleLabels.spike))-1)))
sapply(spikeList,length) # check lengths of list

# now do grain vs spike comparison
# now do for both lists against both lists
gom.obj_grain_spike <- newGOM(grainList, spikeList, background_size)

pvalue_matrix_grain_spike <- getMatrix(gom.obj_grain_spike, name="pval") # get pvalue as a table
Jaccard_matrix_grain_spike <- getMatrix(gom.obj_grain_spike, name="Jaccard") # get Jaccard as table

# we can also plot showing the Jaccard index as the colour
pdf(file="grain_spike_module_comparison_fisher_test_coloured_by_Jaccard.pdf", width=50, height=50)
drawHeatmap(gom.obj_grain_spike, what="Jaccard", ncolused = 5, grid.col = "Blues", note.col = "black")
dev.off()

write.csv(file="pvalues_matrix_grain_spike.csv",pvalue_matrix_grain_spike )
write.csv(file="Jaccard_matrix_grain_spike.csv",Jaccard_matrix_grain_spike )


#### grain vs leaf ####
background_size <- nrow(modules) # need the number of genes in the background
background_size

# make grainList which is a list of character vectors containing gene names for each module
grainList <- list()
for (i in 0:(length(unique(modules$bwnetModuleLabels.grain))-1)) {
  print(i)
  grainList[[i+1]] <- assign(paste0("grain",i), modules[modules$bwnetModuleLabels.grain ==i,][,1])
}

names(grainList) <- c(paste0("grain",0:(length(unique(modules$bwnetModuleLabels.grain))-1)))
sapply(grainList,length) # check lengths of list


# make spikeList which is a list of character vectors containing gene names for each module
leafList <- list()
for (i in 0:(length(unique(modules$bwnetModuleLabels.leaf))-1)) {
  print(i)
  leafList[[i+1]] <- assign(paste0("leaf",i), modules[modules$bwnetModuleLabels.leaf ==i,][,1])
}

names(leafList) <- c(paste0("leaf",0:(length(unique(modules$bwnetModuleLabels.leaf))-1)))
sapply(leafList,length) # check lengths of list

# now do grain vs leaf comparison
# now do for both lists against both lists
gom.obj_grain_leaf <- newGOM(grainList, leafList, background_size)

pvalue_matrix_grain_leaf <- getMatrix(gom.obj_grain_leaf, name="pval") # get pvalue as a table
Jaccard_matrix_grain_leaf <- getMatrix(gom.obj_grain_leaf, name="Jaccard") # get Jaccard as table

# we can also plot showing the Jaccard index as the colour
pdf(file="grain_leaf_module_comparison_fisher_test_coloured_by_Jaccard.pdf", width=50, height=50)
drawHeatmap(gom.obj_grain_leaf, what="Jaccard", ncolused = 5, grid.col = "Blues", note.col = "black")
dev.off()

write.csv(file="pvalues_matrix_grain_leaf.csv",pvalue_matrix_grain_leaf )
write.csv(file="Jaccard_matrix_grain_leaf.csv",Jaccard_matrix_grain_leaf )

#### grain vs root ####
background_size <- nrow(modules) # need the number of genes in the background
background_size

# make grainList which is a list of character vectors containing gene names for each module
grainList <- list()
for (i in 0:(length(unique(modules$bwnetModuleLabels.grain))-1)) {
  print(i)
  grainList[[i+1]] <- assign(paste0("grain",i), modules[modules$bwnetModuleLabels.grain ==i,][,1])
}

names(grainList) <- c(paste0("grain",0:(length(unique(modules$bwnetModuleLabels.grain))-1)))
sapply(grainList,length) # check lengths of list


# make spikeList which is a list of character vectors containing gene names for each module
rootList <- list()
for (i in 0:(length(unique(modules$bwnetModuleLabels.root))-1)) {
  print(i)
  rootList[[i+1]] <- assign(paste0("root",i), modules[modules$bwnetModuleLabels.root ==i,][,1])
}

names(rootList) <- c(paste0("root",0:(length(unique(modules$bwnetModuleLabels.root))-1)))
sapply(rootList,length) # check lengths of list

# now do grain vs root comparison
# now do for both lists against both lists
gom.obj_grain_root <- newGOM(grainList, rootList, background_size)

pvalue_matrix_grain_root <- getMatrix(gom.obj_grain_root, name="pval") # get pvalue as a table
Jaccard_matrix_grain_root <- getMatrix(gom.obj_grain_root, name="Jaccard") # get Jaccard as table

# we can also plot showing the Jaccard index as the colour
pdf(file="grain_root_module_comparison_fisher_test_coloured_by_Jaccard.pdf", width=50, height=50)
drawHeatmap(gom.obj_grain_root, what="Jaccard", ncolused = 5, grid.col = "Blues", note.col = "black")
dev.off()

write.csv(file="pvalues_matrix_grain_root.csv",pvalue_matrix_grain_root )
write.csv(file="Jaccard_matrix_grain_root.csv",Jaccard_matrix_grain_root )

#### spike vs root ####

# now do spike vs root comparison
# now do for both lists against both lists
gom.obj_spike_root <- newGOM(spikeList, rootList, background_size)

pvalue_matrix_spike_root <- getMatrix(gom.obj_spike_root, name="pval") # get pvalue as a table
Jaccard_matrix_spike_root <- getMatrix(gom.obj_spike_root, name="Jaccard") # get Jaccard as table

# we can also plot showing the Jaccard index as the colour
pdf(file="spike_root_module_comparison_fisher_test_coloured_by_Jaccard.pdf", width=50, height=50)
drawHeatmap(gom.obj_spike_root, what="Jaccard", ncolused = 5, grid.col = "Blues", note.col = "black")
dev.off()

write.csv(file="pvalues_matrix_spike_root.csv",pvalue_matrix_spike_root )
write.csv(file="Jaccard_matrix_spike_root.csv",Jaccard_matrix_spike_root )

#### spike vs leaf ####

# now do spike vs leaf comparison
# now do for both lists against both lists
gom.obj_spike_leaf <- newGOM(spikeList, leafList, background_size)

pvalue_matrix_spike_leaf <- getMatrix(gom.obj_spike_leaf, name="pval") # get pvalue as a table
Jaccard_matrix_spike_leaf <- getMatrix(gom.obj_spike_leaf, name="Jaccard") # get Jaccard as table

# we can also plot showing the Jaccard index as the colour
pdf(file="spike_leaf_module_comparison_fisher_test_coloured_by_Jaccard.pdf", width=50, height=50)
drawHeatmap(gom.obj_spike_leaf, what="Jaccard", ncolused = 5, grid.col = "Blues", note.col = "black")
dev.off()

write.csv(file="pvalues_matrix_spike_leaf.csv",pvalue_matrix_spike_leaf )
write.csv(file="Jaccard_matrix_spike_leaf.csv",Jaccard_matrix_spike_leaf )

#### spike vs leaf ####

# now do root vs leaf comparison
# now do for both lists against both lists
gom.obj_root_leaf <- newGOM(rootList, leafList, background_size)

pvalue_matrix_root_leaf <- getMatrix(gom.obj_root_leaf, name="pval") # get pvalue as a table
Jaccard_matrix_root_leaf <- getMatrix(gom.obj_root_leaf, name="Jaccard") # get Jaccard as table

# we can also plot showing the Jaccard index as the colour
pdf(file="root_leaf_module_comparison_fisher_test_coloured_by_Jaccard.pdf", width=50, height=50)
drawHeatmap(gom.obj_root_leaf, what="Jaccard", ncolused = 5, grid.col = "Blues", note.col = "black")
dev.off()

write.csv(file="pvalues_matrix_root_leaf.csv",pvalue_matrix_root_leaf )
write.csv(file="Jaccard_matrix_root_leaf.csv",Jaccard_matrix_root_leaf )


## now want to FDR adjust the pvalues

# root_leaf
pvalue_root_leaf <- read.csv(file="pvalues_matrix_root_leaf.csv")
pvalue_root_leaf[1:4,1:4]
rownames(pvalue_root_leaf) <- pvalue_root_leaf$X
pvalue_root_leaf <- pvalue_root_leaf[,-1]
pvalue_root_leaf[1:4,1:4]
pvalue_matrix_root_leaf <- as.matrix(pvalue_root_leaf)
pvalue_matrix_root_leaf[1:4,1:4]
dim(pvalue_matrix_root_leaf)

padj_matrix_root_leaf <- p.adjust(pvalue_matrix_root_leaf, method= "BY")
length(padj_matrix_root_leaf)
padj_matrix_root_leaf <- matrix(padj_matrix_root_leaf, ncol=ncol(pvalue_matrix_root_leaf))
padj_matrix_root_leaf[1:4,1:4]
padj_matrix_root_leaf <- as.data.frame(padj_matrix_root_leaf)
colnames(padj_matrix_root_leaf) <- colnames(pvalue_matrix_root_leaf)
rownames(padj_matrix_root_leaf) <- rownames(pvalue_root_leaf)
padj_matrix_root_leaf[1:4,1:4]

write.csv(file="padj_matrix_root_leaf.csv",padj_matrix_root_leaf)

#grain_leaf
pvalue_grain_leaf <- read.csv(file="pvalues_matrix_grain_leaf.csv")
pvalue_grain_leaf[1:4,1:4]
rownames(pvalue_grain_leaf) <- pvalue_grain_leaf$X
pvalue_grain_leaf <- pvalue_grain_leaf[,-1]
pvalue_grain_leaf[1:4,1:4]
pvalue_matrix_grain_leaf <- as.matrix(pvalue_grain_leaf)
pvalue_matrix_grain_leaf[1:4,1:4]
dim(pvalue_matrix_grain_leaf)

padj_matrix_grain_leaf <- p.adjust(pvalue_matrix_grain_leaf, method= "BY")
length(padj_matrix_grain_leaf)
padj_matrix_grain_leaf <- matrix(padj_matrix_grain_leaf, ncol=ncol(pvalue_matrix_grain_leaf))
padj_matrix_grain_leaf[1:4,1:4]
padj_matrix_grain_leaf <- as.data.frame(padj_matrix_grain_leaf)
colnames(padj_matrix_grain_leaf) <- colnames(pvalue_matrix_grain_leaf)
rownames(padj_matrix_grain_leaf) <- rownames(pvalue_grain_leaf)
padj_matrix_grain_leaf[1:4,1:4]

write.csv(file="padj_matrix_grain_leaf.csv",padj_matrix_grain_leaf)

#grain_root
pvalue_grain_root <- read.csv(file="pvalues_matrix_grain_root.csv")
pvalue_grain_root[1:4,1:4]
rownames(pvalue_grain_root) <- pvalue_grain_root$X
pvalue_grain_root <- pvalue_grain_root[,-1]
pvalue_grain_root[1:4,1:4]
pvalue_matrix_grain_root <- as.matrix(pvalue_grain_root)
pvalue_matrix_grain_root[1:4,1:4]
dim(pvalue_matrix_grain_root)

padj_matrix_grain_root <- p.adjust(pvalue_matrix_grain_root, method= "BY")
length(padj_matrix_grain_root)
padj_matrix_grain_root <- matrix(padj_matrix_grain_root, ncol=ncol(pvalue_matrix_grain_root))
padj_matrix_grain_root[1:4,1:4]
padj_matrix_grain_root <- as.data.frame(padj_matrix_grain_root)
colnames(padj_matrix_grain_root) <- colnames(pvalue_matrix_grain_root)
rownames(padj_matrix_grain_root) <- rownames(pvalue_grain_root)
padj_matrix_grain_root[1:4,1:4]

write.csv(file="padj_matrix_grain_root.csv",padj_matrix_grain_root)

#grain_spike
pvalue_grain_spike <- read.csv(file="pvalues_matrix_grain_spike.csv")
pvalue_grain_spike[1:4,1:4]
rownames(pvalue_grain_spike) <- pvalue_grain_spike$X
pvalue_grain_spike <- pvalue_grain_spike[,-1]
pvalue_grain_spike[1:4,1:4]
pvalue_matrix_grain_spike <- as.matrix(pvalue_grain_spike)
pvalue_matrix_grain_spike[1:4,1:4]
dim(pvalue_matrix_grain_spike)

padj_matrix_grain_spike <- p.adjust(pvalue_matrix_grain_spike, method= "BY")
length(padj_matrix_grain_spike)
padj_matrix_grain_spike <- matrix(padj_matrix_grain_spike, ncol=ncol(pvalue_matrix_grain_spike))
padj_matrix_grain_spike[1:4,1:4]
padj_matrix_grain_spike <- as.data.frame(padj_matrix_grain_spike)
colnames(padj_matrix_grain_spike) <- colnames(pvalue_matrix_grain_spike)
rownames(padj_matrix_grain_spike) <- rownames(pvalue_grain_spike)
padj_matrix_grain_spike[1:4,1:4]

write.csv(file="padj_matrix_grain_spike.csv",padj_matrix_grain_spike)

#spike_leaf
pvalue_spike_leaf <- read.csv(file="pvalues_matrix_spike_leaf.csv")
pvalue_spike_leaf[1:4,1:4]
rownames(pvalue_spike_leaf) <- pvalue_spike_leaf$X
pvalue_spike_leaf <- pvalue_spike_leaf[,-1]
pvalue_spike_leaf[1:4,1:4]
pvalue_matrix_spike_leaf <- as.matrix(pvalue_spike_leaf)
pvalue_matrix_spike_leaf[1:4,1:4]
dim(pvalue_matrix_spike_leaf)

padj_matrix_spike_leaf <- p.adjust(pvalue_matrix_spike_leaf, method= "BY")
length(padj_matrix_spike_leaf)
padj_matrix_spike_leaf <- matrix(padj_matrix_spike_leaf, ncol=ncol(pvalue_matrix_spike_leaf))
padj_matrix_spike_leaf[1:4,1:4]
padj_matrix_spike_leaf <- as.data.frame(padj_matrix_spike_leaf)
colnames(padj_matrix_spike_leaf) <- colnames(pvalue_matrix_spike_leaf)
rownames(padj_matrix_spike_leaf) <- rownames(pvalue_spike_leaf)
padj_matrix_spike_leaf[1:4,1:4]

write.csv(file="padj_matrix_spike_leaf.csv",padj_matrix_spike_leaf)

#spike_root
pvalue_spike_root <- read.csv(file="pvalues_matrix_spike_root.csv")
pvalue_spike_root[1:4,1:4]
rownames(pvalue_spike_root) <- pvalue_spike_root$X
pvalue_spike_root <- pvalue_spike_root[,-1]
pvalue_spike_root[1:4,1:4]
pvalue_matrix_spike_root <- as.matrix(pvalue_spike_root)
pvalue_matrix_spike_root[1:4,1:4]
dim(pvalue_matrix_spike_root)

padj_matrix_spike_root <- p.adjust(pvalue_matrix_spike_root, method= "BY")
length(padj_matrix_spike_root)
padj_matrix_spike_root <- matrix(padj_matrix_spike_root, ncol=ncol(pvalue_matrix_spike_root))
padj_matrix_spike_root[1:4,1:4]
padj_matrix_spike_root <- as.data.frame(padj_matrix_spike_root)
colnames(padj_matrix_spike_root) <- colnames(pvalue_matrix_spike_root)
rownames(padj_matrix_spike_root) <- rownames(pvalue_spike_root)
padj_matrix_spike_root[1:4,1:4]

write.csv(file="padj_matrix_spike_root.csv",padj_matrix_spike_root)