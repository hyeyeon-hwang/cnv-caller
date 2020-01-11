#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DNAcopy")

# Read in CNV file for now
# Will integrate into CNV pipeline later

library(DNAcopy)
library(magrittr) # pipe
library(tibble) # as_tibble(), add_column()

inputCnvFile = "Dec11_cnv_bins_full.txt"
sampleNames = c("chr", "start", "end",
                "cJLKD", "c6978", "c6980", "c7015", "c7016",
                "e7005", "e7006", "e7007", "e7008")
data <- read.delim(inputCnvFile, header = TRUE, col.names = sampleNames)

mat <- data[,-c(1:3)] %>% as.matrix()

CNA.object <- CNA(genomdat = mat,
                  chrom = data$chr,
                  maploc = data$start,
                  data.type = "logratio",
                  sampleid = colnames(data)[-c(1:3)],
                  presorted = TRUE) 
CNA.object
smoothed.CNA.object <- smooth.CNA(CNA.object)
segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1) # < 10 min

