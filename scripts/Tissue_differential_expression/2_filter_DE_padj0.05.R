# Aim is to filter DE results for all intermediate vs all intermediate for padj <0.05 (to reduce file size)

setwd("Y:\\expression_browser\\WGA\\Bayer_analysis\\DE\\")

files <- list.files(pattern = "_results_padj_0.05.csv$")
files

#f <- files[1]
#f


for (f in files) {
  res <- read.csv(file=f)
  head(res)
  res_NA_rm <- na.omit(res)
  dim(res_NA_rm)
  head(res_NA_rm)
  tail(res_NA_rm)
  res_filt <- res_NA_rm[res_NA_rm$padj <0.05 ,]
  dim(res)
  dim(res_filt)
  head(res_filt)
  tail(res_filt)
  write.csv(file=paste0("padj0.05_DE_for_daniel\\",f), res_filt, row.names=F)
}
