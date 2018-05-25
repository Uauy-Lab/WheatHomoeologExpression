# Calculate how far apart the modules are for homoeologs
# 14.7.2017
# Philippa Borrill

list_dirs <- c("tissues\\grain","tissues\\leaf","tissues\\root","tissues\\spike")
list_names <- c("tissues\\grain" = "grain","tissues\\leaf" = "leaf","tissues\\root" = "root","tissues\\spike" = "spike")


for (dir in list_dirs) {
  

setwd(paste0("Y:\\expression_browser\\WGA\\WGCNA\\",dir,"\\maxP0.05\\homoeologues_paralogues_numeric_threshold_correct"))

########  get distances between modules and plot thresholds of nearest X modules #######

MEs <- read.table(file=paste0("Y:\\expression_browser\\WGA\\WGCNA\\",dir,"\\maxP0.05\\eigengenes.txt"), header=TRUE)
head(MEs)

# remove ME0
MEs <- MEs[,-1]
head(MEs)
distance <- dist(t(MEs), method = "euclidean")
head(distance)

distance_matrix <- as.matrix(distance)
head(distance_matrix)

## now want to choose different thresholds using e.g. 5 most closely related eigengenes

# need to calculate the distance to N closest eigengenes
head(distance_matrix)
min(distance_matrix)
names(distance_matrix[,1])[1]
mean(sort(distance_matrix[,3])[2:5])


# make empty for the for loop
eigengene_dist <- numeric()
eigengene_dist



dim(MEs)

eigengene_dist_result <- data.frame(matrix(ncol=(ncol(distance_matrix)-1),nrow =0))
eigengene_dist_result

for (i in 1:(nrow(distance_matrix))) {
  
  # for each ME calculate the distance to the 1 nearest, 2 nearest etc up to 37 nearest eigengenes
  for (j in 2:(ncol(distance_matrix))) {
    mean_dist <- mean(sort(distance_matrix[,i])[2:j])
    eigengene_dist <- append(eigengene_dist, mean_dist)
    eigengene_dist
  }
  eigengene_dist
  eigengene_dist_result[nrow(eigengene_dist_result)+1,] <- (eigengene_dist)
  eigengene_dist <- numeric()
}
colnames(eigengene_dist_result) <- paste0("nearest",1:(ncol(distance_matrix)-1))
rownames(eigengene_dist_result) <- paste0("ME",1:(ncol(distance_matrix)))
eigengene_dist_result[1:4,1:4]

pdf(file="distance_to_nearest_x_eigengenes.pdf",width=10)
boxplot(eigengene_dist_result,las=2)
dev.off()

######### now want to analyse how far homoeologues move ############
# read in homoeologue info table
hom.df <- read.csv(file="Y:\\expression_browser\\WGA\\data_tables\\wheat.homeolog_groups.release.nonTE.csv", header=T)
head(hom.df)

# only keep HC
hom.df <- hom.df[hom.df$HC.LC =="HC-only",]
head(hom.df)
dim(hom.df)

hom.df <- hom.df[,6:10]
head(hom.df)

# read in module info
module_info <- read.csv(file="Y:\\expression_browser\\WGA\\data_tables\\WGCNA_table.csv", 
                        header=T)
head(module_info)
module_info <- module_info[module_info$set == paste0("WGCNA_",list_names[dir]),1:3]
colnames(module_info) <- c("gene", "module", "bwnetModuleColors")
head(module_info)
dim(module_info)
module_info <- module_info[,1:2] # want to have 2 columns: gene and module
head(module_info)

####### for 1:1:1 syntenic ###########

# choose which homoeologues to look at 1:1:1 syntenic
head(hom.df)
hom.df_1_1_1 <- hom.df[hom.df$cardinality_abs == "1:1:1",]
dim(hom.df_1_1_1)
hom.df_1_1_1_synt <- hom.df_1_1_1[hom.df_1_1_1$synteny == "segmental homeologs",]
dim(hom.df_1_1_1_synt)
head(hom.df_1_1_1_synt)
# merge homoeologue info with module info

merged.df <- merge(hom.df_1_1_1_synt[,1:4],module_info, by.x= "A", by.y = "gene")
head(merged.df)

merged.df <- merge(merged.df,module_info, by.x= "B", by.y = "gene", suffixes = c(".A", ".B"))
head(merged.df)

merged.df <- merge(merged.df,module_info, by.x= "D", by.y = "gene")
head(merged.df)

colnames(merged.df) <- c("D","B","A", "cardinality_abs", "module.A", "module.B", "module.D")
merged.df <- cbind(merged.df[,4], merged.df[,3:1], merged.df[,5:7])
head(merged.df)
colnames(merged.df)[1] <- c("cardinality_abs")
head(merged.df)
merged.df$cardinality_abs <- merged.df$A
colnames(merged.df)[1] <- c("identifier")
head(merged.df)
dim(merged.df)

all_exp_triads <- nrow(merged.df)

# now remove module 0 (i.e. not in a module)
merged.df.filt <- merged.df[merged.df$module.A != "0" |merged.df$module.B != "0" | merged.df$module.D != "0", ]
dim(merged.df.filt)

triads_removed_all_hom_module0 <- nrow(merged.df.filt)

num_triads_1homoeologues_module0 <- nrow(merged.df[(merged.df$module.A == "0"  & merged.df$module.B != "0" & merged.df$module.D != "0")|
                                                     (merged.df$module.B == "0"  & merged.df$module.D != "0" & merged.df$module.A != "0")|
                                                     (merged.df$module.D == "0"  & merged.df$module.A != "0" & merged.df$module.B != "0"), ])
num_triads_1homoeologues_module0

num_triads_2homoeologues_module0 <- nrow(merged.df[(merged.df$module.A == "0"  & merged.df$module.B == "0" & merged.df$module.D != "0")|
                                                     (merged.df$module.B == "0"  & merged.df$module.D == "0" & merged.df$module.A != "0")|
                                                     (merged.df$module.A == "0"  & merged.df$module.D == "0" & merged.df$module.B != "0"), ])
num_triads_2homoeologues_module0

num_triads_3homoeologues_module0 <- nrow(merged.df[(merged.df$module.A == "0"  & merged.df$module.B == "0" & merged.df$module.D == "0")|
                                                     (merged.df$module.B == "0"  & merged.df$module.D == "0" & merged.df$module.A == "0")|
                                                     (merged.df$module.A == "0"  & merged.df$module.D == "0" & merged.df$module.B == "0"), ])
num_triads_3homoeologues_module0

# now remove any triad with a gene in module 0
merged.df.filt <- merged.df[merged.df$module.A != "0" & merged.df$module.B != "0" & merged.df$module.D != "0", ]
dim(merged.df.filt)

triads_removed_any_hom_module0 <- nrow(merged.df.filt)

num_triads_with_hom_in_module0 <- all_exp_triads - triads_removed_any_hom_module0
num_triads_with_hom_in_module0
num_triads_with_hom_in_module0 / all_exp_triads

#### now calculate distances between homoeologues

head(merged.df.filt)

# can extract value from distance matrix e.g. for module 5 and 4
distance_matrix[5,4]


# columns 5,6 and 7 contain distances between homoeologues 
distance_matrix[merged.df.filt[2,5],merged.df.filt[2,7]]

# now loop through all rows to calculate distances between each homoeologue

homoeologue_dist <- data.frame(identifier=character(),distAB=numeric(), distBD=numeric(), distAD=numeric(),  stringsAsFactors=F)
homoeologue_dist
for (k in 1:nrow(merged.df.filt)) {
  
  distAB <- distance_matrix[merged.df.filt[k,5],merged.df.filt[k,6]]
  distBD <- distance_matrix[merged.df.filt[k,6],merged.df.filt[k,7]]
  distAD <- distance_matrix[merged.df.filt[k,5],merged.df.filt[k,7]]
  identifier <- merged.df.filt[k,1]
  
  homoeologue_dist[nrow(homoeologue_dist)+1,] <- c(as.character(identifier),distAB, distBD, distAD)
  
}

head(homoeologue_dist)
homoeologue_dist$distAB <- as.numeric(homoeologue_dist$distAB)
homoeologue_dist$distBD <- as.numeric(homoeologue_dist$distBD)
homoeologue_dist$distAD <- as.numeric(homoeologue_dist$distAD)
head(homoeologue_dist)

write.csv(file="homoeologue_distances_inc_same_module_1_1_1_syntenic.csv",homoeologue_dist)

# now remove any triad with a distance of 0 across all three homoeologues
homoeologue_dist.filt <- homoeologue_dist[homoeologue_dist$distAB != "0" | homoeologue_dist$distBD != "0" | homoeologue_dist$distAD != "0", ]
dim(homoeologue_dist.filt)
head(homoeologue_dist.filt)

write.csv(file="homoeologue_distances_different_module_1_1_1_syntenic.csv",homoeologue_dist.filt)


### now want to see how may triads move over the median distance for each category ###

all_thresholds_mvmt <- data.frame(identifier=character(), threshold=numeric(),AB=numeric(), BD=numeric(), AD=numeric(), AB_or_BD_or_AD=numeric(), AB_and_BD_and_AD=numeric(),  stringsAsFactors=F)
all_thresholds_mvmt

# make vector with thresholds
# 100 % threshold is 1.874862
max_threshold <- 1.874862

medians <- seq(from=0, to=max_threshold, length.out=101)
names(medians) <- seq(from=0, to=100, length.out=101)
medians

# add first value which is negative so you can also count same module homoeologues
first_value <- -0.1
names(first_value) <- "same"
medians <- c(first_value,medians)
medians
is.vector(medians)

# need a counter in the loop
counter <- 0

for (threshold in medians) {
  
  AB_move_further <- nrow(homoeologue_dist[homoeologue_dist$distAB > threshold,])
  BD_move_further <- nrow(homoeologue_dist[homoeologue_dist$distBD > threshold,])
  AD_move_further <- nrow(homoeologue_dist[homoeologue_dist$distAD > threshold,])
  
  # either AB BD or AD is larger
  either_AB_or_BD_or_AD_move_further <- nrow(homoeologue_dist[homoeologue_dist$distAB > threshold | homoeologue_dist$distBD > threshold | homoeologue_dist$distAD > threshold  ,])
  
  #  AB BD AND AD are larger
  all_AB_BD_AD_move_further <- nrow(homoeologue_dist[homoeologue_dist$distAB > threshold & homoeologue_dist$distBD > threshold & homoeologue_dist$distAD > threshold  ,])
  
  counter <- counter +1
  
  identifier <- names(medians)[counter]
  
  all_thresholds_mvmt[nrow(all_thresholds_mvmt)+1,] <- c(as.character(identifier),threshold, 
                                                         AB_move_further,BD_move_further,AD_move_further,
                                                         either_AB_or_BD_or_AD_move_further, all_AB_BD_AD_move_further)
  
}

all_thresholds_mvmt
num_triads_homoeologues_alloc <- nrow(merged.df.filt) # number of triads where all three homoeologues were in a module

all_thresholds_mvmt$total_num_triads <- nrow(homoeologue_dist)
all_thresholds_mvmt$two_hom_further_than_threshold <- as.numeric(all_thresholds_mvmt$AB_or_BD_or_AD)/num_triads_homoeologues_alloc
all_thresholds_mvmt$three_hom_further_than_threshold <- as.numeric(all_thresholds_mvmt$AB_and_BD_and_AD)/num_triads_homoeologues_alloc
all_thresholds_mvmt$type <- "1:1:1"
all_thresholds_mvmt$synteny <- "syntenic"
all_thresholds_mvmt$set <- paste(all_thresholds_mvmt$type,all_thresholds_mvmt$synteny, sep = " ")

write.csv(file="median_eigengene_distances_to_nearest_x_eigengenes_homoeologue_movements_1_1_1_syntenic.csv", all_thresholds_mvmt)

####### for 1:1:1 non-syntenic ###########

head(hom.df)
hom.df_1_1_1 <- hom.df[hom.df$cardinality_abs == "1:1:1",]
dim(hom.df_1_1_1)
hom.df_1_1_1_nonsynt <- hom.df_1_1_1[hom.df_1_1_1$synteny != "segmental homeologs",]
dim(hom.df_1_1_1_nonsynt)
head(hom.df_1_1_1_nonsynt)
# merge homoeologue info with module info

merged.df <- merge(hom.df_1_1_1_nonsynt[,1:4],module_info, by.x= "A", by.y = "gene")
head(merged.df)

merged.df <- merge(merged.df,module_info, by.x= "B", by.y = "gene", suffixes = c(".A", ".B"))
head(merged.df)

merged.df <- merge(merged.df,module_info, by.x= "D", by.y = "gene")
head(merged.df)

colnames(merged.df) <- c("D","B","A", "cardinality_abs", "module.A", "module.B", "module.D")
merged.df <- cbind(merged.df[,4], merged.df[,3:1], merged.df[,5:7])
head(merged.df)
colnames(merged.df)[1] <- c("cardinality_abs")
head(merged.df)
merged.df$cardinality_abs <- merged.df$A
colnames(merged.df)[1] <- c("identifier")
head(merged.df)
dim(merged.df)

# now remove module 0 (i.e. not in a module)
merged.df.filt <- merged.df[merged.df$module.A != "0" |merged.df$module.B != "0" | merged.df$module.D != "0", ]
dim(merged.df.filt)

# now remove any triad with a gene in module 0
merged.df.filt <- merged.df[merged.df$module.A != "0" & merged.df$module.B != "0" & merged.df$module.D != "0", ]
dim(merged.df.filt)

#### now calculate distances between homoeologues

head(merged.df.filt)

# can extract value from distance matrix e.g. for module 5 and 4
distance_matrix[5,4]

# columns 5,6 and 7 contain distances between homoeologues 
distance_matrix[merged.df.filt[2,5],merged.df.filt[2,7]]

# now loop through all rows to calculate distances between each homoeologue

homoeologue_dist <- data.frame(identifier=character(),distAB=numeric(), distBD=numeric(), distAD=numeric(),  stringsAsFactors=F)
homoeologue_dist
for (k in 1:nrow(merged.df.filt)) {
  
  distAB <- distance_matrix[merged.df.filt[k,5],merged.df.filt[k,6]]
  distBD <- distance_matrix[merged.df.filt[k,6],merged.df.filt[k,7]]
  distAD <- distance_matrix[merged.df.filt[k,5],merged.df.filt[k,7]]
  identifier <- merged.df.filt[k,1]
  
  homoeologue_dist[nrow(homoeologue_dist)+1,] <- c(as.character(identifier),distAB, distBD, distAD)
  
}

head(homoeologue_dist)
homoeologue_dist$distAB <- as.numeric(homoeologue_dist$distAB)
homoeologue_dist$distBD <- as.numeric(homoeologue_dist$distBD)
homoeologue_dist$distAD <- as.numeric(homoeologue_dist$distAD)
head(homoeologue_dist)

write.csv(file="homoeologue_distances_inc_same_module_1_1_1_non-syntenic.csv",homoeologue_dist)


# now remove any triad with a distance of 0 across all three homoeologues
homoeologue_dist.filt <- homoeologue_dist[homoeologue_dist$distAB != "0" | homoeologue_dist$distBD != "0" | homoeologue_dist$distAD != "0", ]
dim(homoeologue_dist.filt)
head(homoeologue_dist.filt)


write.csv(file="homoeologue_distances_different_module_1_1_1_non-syntenic.csv",homoeologue_dist.filt)

### now want to see how may triads move over the median distance for each category ###

all_thresholds_mvmt <- data.frame(identifier=character(), threshold=numeric(),AB=numeric(), BD=numeric(), AD=numeric(), AB_or_BD_or_AD=numeric(), AB_and_BD_and_AD=numeric(),  stringsAsFactors=F)
all_thresholds_mvmt

# need a counter in the loop
counter <- 0

for (threshold in medians) {
  
  AB_move_further <- nrow(homoeologue_dist[homoeologue_dist$distAB > threshold,])
  BD_move_further <- nrow(homoeologue_dist[homoeologue_dist$distBD > threshold,])
  AD_move_further <- nrow(homoeologue_dist[homoeologue_dist$distAD > threshold,])
  
  # either AB BD or AD is larger
  either_AB_or_BD_or_AD_move_further <- nrow(homoeologue_dist[homoeologue_dist$distAB > threshold | homoeologue_dist$distBD > threshold | homoeologue_dist$distAD > threshold  ,])
  
  #  AB BD AND AD are larger
  all_AB_BD_AD_move_further <- nrow(homoeologue_dist[homoeologue_dist$distAB > threshold & homoeologue_dist$distBD > threshold & homoeologue_dist$distAD > threshold  ,])
  
  counter <- counter +1
  
  identifier <- names(medians)[counter]
  
  all_thresholds_mvmt[nrow(all_thresholds_mvmt)+1,] <- c(as.character(identifier),threshold, 
                                                         AB_move_further,BD_move_further,AD_move_further,
                                                         either_AB_or_BD_or_AD_move_further, all_AB_BD_AD_move_further)
  
}

all_thresholds_mvmt
num_triads_homoeologues_alloc <- nrow(merged.df.filt) # number of triads where all three homoeologues were in a module

all_thresholds_mvmt$total_num_triads <- nrow(homoeologue_dist)
all_thresholds_mvmt$two_hom_further_than_threshold <- as.numeric(all_thresholds_mvmt$AB_or_BD_or_AD)/num_triads_homoeologues_alloc
all_thresholds_mvmt$three_hom_further_than_threshold <- as.numeric(all_thresholds_mvmt$AB_and_BD_and_AD)/num_triads_homoeologues_alloc
all_thresholds_mvmt$type <- "1:1:1"
all_thresholds_mvmt$synteny <- "non-syntenic"
all_thresholds_mvmt$set <- paste(all_thresholds_mvmt$type,all_thresholds_mvmt$synteny, sep = " ")

write.csv(file="median_eigengene_distances_to_nearest_x_eigengenes_homoeologue_movements_1_1_1_non-syntenic.csv", all_thresholds_mvmt)


####### for random triads ###########

head(module_info)
dim(module_info)

module_info_no_0 <- module_info[module_info$module != 0,]
dim(module_info_no_0)

module_info_A <- module_info_no_0[(grepl("A",module_info_no_0$gene)),]
head(module_info_A)

module_info_B <- module_info_no_0[(grepl("B",module_info_no_0$gene)),]
head(module_info_B)

module_info_D <- module_info_no_0[(grepl("D",module_info_no_0$gene)),]
head(module_info_D)

# make a random set of triads
set.seed(0)
random_triads <- (data.frame(A = module_info_A[sample(nrow(module_info_A), 1000), 1], 
                             B =  module_info_B[sample(nrow(module_info_B), 1000),1 ], 
                             D = module_info_D[sample(nrow(module_info_D), 1000), 1]))

head(random_triads)
random_triads$cardinality_abs <- "1:1:1"
head(random_triads)

# merge homoeologue info with module info

merged.df <- merge(random_triads[,1:4],module_info, by.x= "A", by.y = "gene")
head(merged.df)

merged.df <- merge(merged.df,module_info, by.x= "B", by.y = "gene", suffixes = c(".A", ".B"))
head(merged.df)

merged.df <- merge(merged.df,module_info, by.x= "D", by.y = "gene")
head(merged.df)

colnames(merged.df)[7] <-"module.D"

merged.df <- cbind(merged.df[,4], merged.df[,3:1], merged.df[,5:7])
head(merged.df)
colnames(merged.df)[1] <- c("cardinality_abs")
head(merged.df)
merged.df$cardinality_abs <- merged.df$A 
colnames(merged.df)[1] <- c("identifier")
head(merged.df)
dim(merged.df)

# now remove module 0 (i.e. not in a module)
merged.df.filt <- merged.df[merged.df$module.A != "0" | merged.df$module.B != "0", ]
dim(merged.df.filt)

# now remove any triad with a gene in module 0
merged.df.filt <- merged.df[merged.df$module.A != "0" & merged.df$module.B != "0", ]
dim(merged.df.filt)

#### now calculate distances between homoeologues

head(merged.df.filt)

# can extract value from distance matrix e.g. for module 5 and 4
distance_matrix[5,4]

# columns 5,6 and 7 contain distances between homoeologues 
distance_matrix[merged.df.filt[2,5],merged.df.filt[2,7]]

# now loop through all rows to calculate distances between each homoeologue

homoeologue_dist <- data.frame(identifier=character(),distAB=numeric(), distBD=numeric(), distAD=numeric(),  stringsAsFactors=F)
homoeologue_dist
for (k in 1:nrow(merged.df.filt)) {
  
  distAB <- distance_matrix[merged.df.filt[k,5],merged.df.filt[k,6]]
  distBD <- distance_matrix[merged.df.filt[k,6],merged.df.filt[k,7]]
  distAD <- distance_matrix[merged.df.filt[k,5],merged.df.filt[k,7]]
  identifier <- merged.df.filt[k,1]
  
  homoeologue_dist[nrow(homoeologue_dist)+1,] <- c(as.character(identifier),distAB, distBD, distAD)
  
}

head(homoeologue_dist)
homoeologue_dist$distAB <- as.numeric(homoeologue_dist$distAB)
homoeologue_dist$distBD <- as.numeric(homoeologue_dist$distBD)
homoeologue_dist$distAD <- as.numeric(homoeologue_dist$distAD)
head(homoeologue_dist)

write.csv(file="homoeologue_distances_inc_same_module_random.csv",homoeologue_dist)

# now remove any triad with a distance of 0 across all three homoeologues
homoeologue_dist.filt <- homoeologue_dist[homoeologue_dist$distAB != "0" | homoeologue_dist$distBD != "0" | homoeologue_dist$distAD != "0", ]
dim(homoeologue_dist.filt)
head(homoeologue_dist.filt)

pdf(file="distances_between_homoeologues_different_modules_random.pdf")
boxplot(homoeologue_dist.filt[,2:4])
text(y = boxplot.stats(homoeologue_dist.filt$distAB)$stats, labels = round(boxplot.stats(homoeologue_dist.filt$distAB)$stats, digits = 4), x = 1.25)
text(y = boxplot.stats(homoeologue_dist.filt$distBD)$stats, labels = round(boxplot.stats(homoeologue_dist.filt$distBD)$stats, digits = 4), x = 2.25)
text(y = boxplot.stats(homoeologue_dist.filt$distAD)$stats, labels = round(boxplot.stats(homoeologue_dist.filt$distAD)$stats, digits = 4), x = 3.25)
dev.off()

write.csv(file="homoeologue_distances_different_module_random.csv",homoeologue_dist.filt)


### now want to see how may triads move over the median distance for each category ###

all_thresholds_mvmt <- data.frame(identifier=character(), threshold=numeric(),AB=numeric(), BD=numeric(), AD=numeric(), AB_or_BD_or_AD=numeric(), AB_and_BD_and_AD=numeric(),  stringsAsFactors=F)
all_thresholds_mvmt

# need a counter in the loop
counter <- 0

for (threshold in medians) {
  
  AB_move_further <- nrow(homoeologue_dist[homoeologue_dist$distAB > threshold,])
  BD_move_further <- nrow(homoeologue_dist[homoeologue_dist$distBD > threshold,])
  AD_move_further <- nrow(homoeologue_dist[homoeologue_dist$distAD > threshold,])
  
  # either AB BD or AD is larger
  either_AB_or_BD_or_AD_move_further <- nrow(homoeologue_dist[homoeologue_dist$distAB > threshold | homoeologue_dist$distBD > threshold | homoeologue_dist$distAD > threshold  ,])
  
  #  AB BD AND AD are larger
  all_AB_BD_AD_move_further <- nrow(homoeologue_dist[homoeologue_dist$distAB > threshold & homoeologue_dist$distBD > threshold & homoeologue_dist$distAD > threshold  ,])
  
  counter <- counter +1
  
  identifier <- names(medians)[counter]
  
  all_thresholds_mvmt[nrow(all_thresholds_mvmt)+1,] <- c(as.character(identifier),threshold, 
                                                         AB_move_further,BD_move_further,AD_move_further,
                                                         either_AB_or_BD_or_AD_move_further, all_AB_BD_AD_move_further)
  
}

all_thresholds_mvmt
num_triads_homoeologues_alloc <- nrow(merged.df.filt) # number of triads where all three homoeologues were in a module

all_thresholds_mvmt$total_num_triads <- nrow(homoeologue_dist)
all_thresholds_mvmt$two_hom_further_than_threshold <- as.numeric(all_thresholds_mvmt$AB_or_BD_or_AD)/num_triads_homoeologues_alloc
all_thresholds_mvmt$three_hom_further_than_threshold <- as.numeric(all_thresholds_mvmt$AB_and_BD_and_AD)/num_triads_homoeologues_alloc
all_thresholds_mvmt$type <- "random"
all_thresholds_mvmt$synteny <- "random"
all_thresholds_mvmt$set <- paste(all_thresholds_mvmt$type,all_thresholds_mvmt$synteny, sep = " ")

write.csv(file="median_eigengene_distances_to_nearest_x_eigengenes_homoeologue_movements_random.csv", all_thresholds_mvmt)

}
