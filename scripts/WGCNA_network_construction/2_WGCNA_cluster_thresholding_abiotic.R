# aim is to run WGCNA on grain samples from RNA-seq
# run on cluster
# 21-09-2016
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
setwd("/nbi/Research-Groups/NBI/Cristobal-Uauy/expression_browser/WGA/WGCNA/abiotic_vs_disease/abiotic")

library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads(10)

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
dim(datExpr)
print("dimensions of datExpr after removing bad genes")
print(dim(datExpr))

# change into new folder to save results
setwd("/nbi/Research-Groups/NBI/Cristobal-Uauy/expression_browser/WGA/WGCNA/abiotic_vs_disease/abiotic/maxP0.05")


# realise I should have run the pickSoftThreshold function with a signed network with bicor correlation therefore re-run
# Choose a set of soft-thresholding powers
powers = c(c(1:16), seq(from = 18, to=20, by=2))
# run soft-thresholding
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType ="signed hybrid", 
                        corFnc ="bicor", corOptions=list(maxPOutliers=0.05))


# Plot the results:
sizeGrWindow(9, 5)
pdf(file="soft-threshold_power_signed_hybrid_0.5tpm.pdf")
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
