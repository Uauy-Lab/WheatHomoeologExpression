# Aim is to run WGCNA on samples from RefSeqv1.0 analysis

# Philippa Borrill
#12-05-2017

#Steps will be:

#1: Summarise counts per gene (rather than transcript) using tximport for the studies which are to be included in the manuscript

#2: Summarise tpm per gene (rather than transcript) using tximport for the studies which are to be included in the manuscript

#3: Filter count data to only keep genes expressed >0.5 tpm in at least 1 % of samples

#4: Use variance stabilising normalisation from DESeq2 to normalise counts

#5: Load in metadata and plot dendrograms of sample clustering

#6: Run soft power threshold

#7: Run WGCNA (try in 1 block) - need to use a different function from usual which uses blocks

#############3: Filter count data to only keep genes expressed >0.5 tpm in at least 1 % of samples############

# In my previous script "WGCNA_tximport_summarise_counts_tpm_per_gene.R" I produce two tsv files per study with counts and tpm
# I now need to combine all the studies I want to use for WGCNA (i.e. all with Y publish excluding aneuploids, nullitetras and transgenics)
# I will use the script which Ricardo has written to load the data

##### Load the data #####

library("matrixStats")
library("DESeq2")
library("WGCNA")
library("dendextend")

setwd("/nbi/Research-Groups/NBI/Cristobal-Uauy/expression_browser/WGA/WGCNA/")

loadValuesFromExperiment<-function(metadata, folder, unit="tpm", values=c("Development")){
  metadata$Sample.IDs <- gsub("-",".",metadata$Sample.IDs)
#  print(metadata$Sample.IDs)
  v<-values[1]
  v<-gsub(" ","_",v)
  path<-paste0(folder,"/",v,"_",unit,".tsv")
  ret<-read.table(path, row.names = 1, header= TRUE)
#  print(v)
  for(i in 2:length(values)){
    v<-values[i]
    v<-gsub(" ","_",v)
 #   print(v)
    path<-paste0(folder,"/",v,"_",unit,".tsv")
    tmp<-read.table(path, row.names = 1, header= TRUE)
    ret<-cbind(ret,tmp)
  }
#  print(colnames(ret))
  md<-metadata[metadata$Sample.IDs%in%colnames(ret),]
  ret<-ret[,as.character(md$Sample.IDs),]
  list(ret,md)
}

folder<-"expressionValuesPerGene"

#metadata_file<- "/nbi/Research-Groups/NBI/Cristobal-Uauy/expression_browser/collaborators/Metadata_may04_FINAL_RR.txt"
metadata_file<- "/nbi/Research-Groups/NBI/Cristobal-Uauy/expression_browser/WGA/WGCNA/Metadata_june2_FINAL_RR_CS_spike.txt"
#metadata_file<-"Metadata_FINAL.txt"

metadata <- read.csv(metadata_file, row.names = 1, sep="\t")
metadata[1:10,1:4]
dim(metadata)
# filter out any samples we don't want to include

metadata <-  metadata[metadata[,"High.level.stress.disease"] != "aneuploidy",]
metadata <-  metadata[metadata[,"High.level.stress.disease"] != "transgenic",]
toFilter<-c("Sample_C18BCACXX_W-C-AR6",
            "HI.3309.007.Index_5.CS_125RNA_14d_Root8",
            "HI.3309.006.Index_7.CS_125RNA_14d_Root9",
            "ERR1141919")
studies_toFilter<-c("ERP009798", "cDNA RenSeq" )
metadata <-  metadata[metadata[,"Variety"] != "synthetic hexaploid (mother Triticum turgidum ssp dicoccon accession PI94655, father Aegilops tauschii ssp strangulata accession AS2404)",]
metadata <-  metadata[! metadata$Sample.IDs%in%toFilter,]
metadata <-  metadata[! metadata$study_title%in%studies_toFilter,]
metadata <-  metadata[metadata[,"Pubish.with.IWGSC..Y.N."] == "Y",]
metadata <-  metadata[metadata[,"High.level.stress.disease"] == "none",]
metadata <-  metadata[metadata[,"High.level.tissue"] == "leaves/shoots",]

# remove because it's an outlier
metadata <-  metadata[metadata[,"Sample.IDs"] != "ERR789099",]



metadata[1:10,1:4]
print("nrow of metadata")
print(nrow(metadata))

tpms  <-loadValuesFromExperiment(metadata, folder, unit="tpm",  values=unique(metadata$study_title))
counts<-loadValuesFromExperiment(metadata, folder, unit="count",values=unique(metadata$study_title))

metadata_used<-tpms[[2]]
tpms<-tpms[[1]]
counts<-counts[[1]]
print("nrow of metadata_used")
print(nrow(metadata_used))

tpms[1:4,1:4]
dim(tpms)
head(colnames(tpms))
head(colnames(counts))
head(rownames(tpms))
head(rownames(counts))

print("dimension of counts")
print(dim(counts))
setwd("/nbi/Research-Groups/NBI/Cristobal-Uauy/expression_browser/WGA/WGCNA/tissues/leaf")
write.csv(metadata_used, file="WGCNA_metadata_cluster.csv")

#3: Filter count data to only keep genes expressed >0.5 tpm in at least 3 samples
tpm_threshold <- 0.5 

# set what is number of samples we need over
perc_samples <- 2


#checked the filtering works as expected 15-05-2017
count_filt <- counts[rowCounts(as.matrix(tpms>tpm_threshold))>perc_samples,] # select only rows which have expr >0.5 tpm in >2 samples
nrow(counts)
print("dimensions of count_filt for tpm >0.5 in 3 samples")
print(dim(count_filt))
nrow(count_filt)

# now want to filter to only keep HC genes
highconf <- read.csv("/nbi/Research-Groups/NBI/Cristobal-Uauy/expression_browser/WGA/WGCNA/high_conf_genes.csv")
head(highconf)

count_filt <- count_filt[rownames(count_filt)%in%highconf$Gene,]
print("dimensions of count_filt for tpm >0.5 in 3 samples and high conf")
print(dim(count_filt))

#4: Use variance stabilising normalisation from DESeq2 to normalise counts
# round counts to integers (required for DESeq)
count_filt <- round(count_filt)
dim(count_filt)

# get data into DESeq form (convert to matrix)
counts_matrix <- as.matrix(count_filt)
counts_matrix[1:4,1:4]

#do variance stabilising transformation (see https://www.biostars.org/p/95788/)
vsd_blind <- varianceStabilizingTransformation(counts_matrix,blind=TRUE)
vc <- (vsd_blind)
vc[1:4,1:4]

# transform the data to be samples as rows and genes as columns (what WGCNA expects)
vsd2 <- t(vc)
vsd2[1:4,1:4]
print("dim of vsd2")
print(dim(vsd2))


#5: Load in metadata and plot dendrograms of sample clustering
# source the WGCNA package

options(stringsAsFactors = FALSE)

# if running in R-gui can use multi-threading
allowWGCNAThreads()

# copy vsd2 to a new dataframe called dataExpr0 to match the manual
datExpr0 <- vsd2

# Run this to check if there are gene outliers
gsg = goodSamplesGenes(datExpr0, verbose = 3)
print("are all genes ok?")
print(gsg$allOK)

#If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data with the following:
if (!gsg$allOK)
   {if (sum(!gsg$goodGenes)>0)
       printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse= ", ")));
       if (sum(!gsg$goodSamples)>0)
           printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse=", ")))
       datExpr0= datExpr0[gsg$goodSamples, gsg$goodGenes]
       }

print("dimensions of datExpr after removing bad genes")
print(dim(datExpr0))

#All my genes were OK

# Now check all samples look ok 
sampleTree = hclust(dist(datExpr0), method = "average");

pdf(file = paste0("Sample_Clustering_",tpm_threshold,"tpm.pdf"), width = 20, height = 4)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()

# If want to add trait data (e.g. chlorophyll level or age or protein content or moisture content need to do this now)
# should add in TISSUE + AGE info
# see 1c in https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-01-dataInput.pdf 


# save the data to use in the next step
save(datExpr0, metadata_used, file=paste0("filtered_data_ready_for_WGCNA_",tpm_threshold,"tpm.RData"))



# now plot dendrogram with colours according to tissue
metadata_selected <- metadata_used

# rename 
metadata_selected$High.level.tissue <- gsub("leaves/shoots", "leavesshoots", metadata_selected$High.level.tissue)


# convert sampleTree to dendrogram

dend <- as.dendrogram(sampleTree)

colorCodes <- c(roots= "brown", leavesshoots = "green", spike = "orange", grain= "purple")
labels_colors(dend) <- colorCodes[(metadata_selected$High.level.tissue[(order.dendrogram(dend))])]
labels_colors(dend)
head(metadata_selected$High.level.tissue)
head(metadata_selected$High.level.tissue[(order.dendrogram(dend))])

labels_colors(dend)

pdf(file = paste0("Sample_Clustering_coloured_by_tissue_",tpm_threshold,"tpm.pdf"), width = 60, height = 15)
par(cex = 0.6);
par(mar = c(6,4,2,0))
plot(dend, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()


# now plot with high level age as colour
metadata_selected$High.level.age
colorCodes <- c(seedling = "pale green", vegetative = "green3", reproductive= "olivedrab")
labels_colors(dend) <- colorCodes[(metadata_selected$High.level.age[(order.dendrogram(dend))])]
labels_colors(dend)
head(metadata_selected$High.level.age)
head(metadata_selected$High.level.age[(order.dendrogram(dend))])

pdf(file = paste0("Sample_Clustering_coloured_by_age_",tpm_threshold,"tpm.pdf"), width = 60, height = 15)
par(cex = 0.6);
par(mar = c(6,4,2,0))
plot(dend, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()



