sexCallerData <- read.delim("sex_caller_output_asd_Jan28.txt", 
                            header = TRUE, 
                            col.names = c("Sample_name", 
                                          "ChrX_reads",
                                          "ChrY_reads",
                                          "ChrY:ChrX_ratio",
                                          "ChrY:ChrX_percent",
                                          "Sex",
                                          "Sample_info_sex")) 


sexCluster <- kmeans(x = sexCallerData$ChrY.ChrX_ratio,
                     centers = 2)
sexCluster$centers
sexCluster$cluster

# Male has larger ratio 
# Find larger ratio and set to male
# Index of sexCluster$centers for male = 1
maleIndex <- which(sexCluster$centers == max(sexCluster$centers))
maleIndex

predictedSex <- character()
for (idx in sexCluster$cluster) {
  if (idx == maleIndex) {
    predictedSex <- c(predictedSex, "M")
  } else {
    predictedSex <- c(predictedSex, "F")
  }
}
predictedSex


