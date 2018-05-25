# Add dominance category to 1:1:1 triads from nullitetras
# 21.2.2018

setwd("Y:\\expression_browser\\WGA\\nullitetras\\1_1_1_triads")

leaf_nullitetra <- read.csv(file="leaf_triad_expr_in_nullis_no_zeroCS.csv")
head(leaf_nullitetra)
# add columns with % A, B and D in CS
leaf_nullitetra$CS_tpmA_perc <- leaf_nullitetra$CS_tpmA/(leaf_nullitetra$CS_tpmA +leaf_nullitetra$CS_tpmB +leaf_nullitetra$CS_tpmD)
leaf_nullitetra$CS_tpmB_perc <- leaf_nullitetra$CS_tpmB/(leaf_nullitetra$CS_tpmA +leaf_nullitetra$CS_tpmB +leaf_nullitetra$CS_tpmD)
leaf_nullitetra$CS_tpmD_perc <- leaf_nullitetra$CS_tpmD/(leaf_nullitetra$CS_tpmA +leaf_nullitetra$CS_tpmB +leaf_nullitetra$CS_tpmD)

leaf_nullitetra <- leaf_nullitetra[,c(1:6,15:17,7:14)]
head(leaf_nullitetra)

leaf_dom <- read.csv(file="LeafCategories.csv")
head(leaf_dom)

leaf_dom_to_merge <- unique(leaf_dom[,1:3])
head(leaf_dom_to_merge)
dim(leaf_dom)
dim(leaf_dom_to_merge)

leaf_nullitetra_dom <- merge(leaf_nullitetra, leaf_dom_to_merge, by="group_id")
head(leaf_nullitetra_dom)

dim(leaf_nullitetra)
dim(leaf_nullitetra_dom)

# save the data 
write.csv(file="leaf_triad_expr_in_nullis_no_zeroCS_with_dominance.csv", leaf_nullitetra_dom, row.names = F)




### now do the same for root ###

root_nullitetra <- read.csv(file="root_triad_expr_in_nullis_no_zeroCS.csv")
head(root_nullitetra)
# add columns with % A, B and D in CS
root_nullitetra$CS_tpmA_perc <- root_nullitetra$CS_tpmA/(root_nullitetra$CS_tpmA +root_nullitetra$CS_tpmB +root_nullitetra$CS_tpmD)
root_nullitetra$CS_tpmB_perc <- root_nullitetra$CS_tpmB/(root_nullitetra$CS_tpmA +root_nullitetra$CS_tpmB +root_nullitetra$CS_tpmD)
root_nullitetra$CS_tpmD_perc <- root_nullitetra$CS_tpmD/(root_nullitetra$CS_tpmA +root_nullitetra$CS_tpmB +root_nullitetra$CS_tpmD)

root_nullitetra <- root_nullitetra[,c(1:6,15:17,7:14)]
head(root_nullitetra)

root_dom <- read.csv(file="RootCategories.csv")
head(root_dom)

root_dom_to_merge <- unique(root_dom[,1:3])
head(root_dom_to_merge)
dim(root_dom)
dim(root_dom_to_merge)

root_nullitetra_dom <- merge(root_nullitetra, root_dom_to_merge, by="group_id")
head(root_nullitetra_dom)

dim(root_nullitetra)
dim(root_nullitetra_dom)

# save the data 
write.csv(file="root_triad_expr_in_nullis_no_zeroCS_with_dominance.csv", root_nullitetra_dom, row.names = F)

