# Title     : TMM and RPKM Normalization
# Objective : Normalize RNAseq counts using EdgeR's TMM/RPKM
# Created by: Rafael
# Created on: 2020-07-13
library(edgeR, limma)

# method to grab gene lengths from biomaRt (takes a vector of gene ids and returns a dataframe containing gene_id and lengths)
getGeneLengths <- function (geneList) {
  if(!file.exists(lengthFilePath)){
    library("biomaRt")
    humanDatabase <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
    geneAttributes <- getBM(attributes =
                            c("ensembl_gene_id", "start_position", "end_position"),
                          filters = "ensembl_gene_id", values = geneList, mart = humanDatabase, uniqueRows = FALSE)
    geneLengths <- data.frame("gene_id" = geneAttributes$ensembl_gene_id, "length" = geneAttributes$end_position - geneAttributes$start_position, stringsAsFactors = FALSE)
    geneLengths <- geneLengths[order(geneLengths$gene_id),]
    write.table(geneLengths, file = lengthFilePath, sep = "\t", row.names = FALSE)
  }
  lengths <- read.delim(lengthFilePath, check.names = FALSE, stringsAsFactors = FALSE)
}

checkGeneLengths <- function (normData, lengths){
  n <- normData
  unused <- NULL
  unusedNum <- NULL

  for (i in seq(1, nrow(n))){
    if(!(noquote(n[i, 1]) %in% lengths$gene_id)){
      unused <- c(unused, noquote(n[i, 1]))
      unusedNum <- c(unusedNum, i)
    }
  }

  for (i in unusedNum){
    index <- i - (match(i, unusedNum) - 1)
    n <- n[-index,]
  }

  write.table(unused, file = paste0(gsub('.txt', '', outFilePath), "_unused.txt"), row.names = FALSE)
  n
}

# reading command line arguments (norm method, input, output)
args <- commandArgs(trailingOnly = TRUE)
normMethod <- args[1]
inFilePath <- args[2]
outFilePath <- args[3]
lengthFilePath <- paste0("lengths/", substr(inFilePath, 8, nchar(inFilePath) - 4), "_lengths.txt") # is there a better way to do this?

# reading in the raw data in text file to a chart and adding it to a DGEList
rawData <- read.delim(inFilePath, check.names = FALSE, stringsAsFactors = FALSE)
normData <- rawData[order(rawData[1]),]

# check which column data begins in
firstCol <- 3 # default
for (i in seq_len(ncol(rawData))) {
  if (sapply(rawData[i], is.numeric)){
    firstCol <- i
    break
  }
}

filtered <- filterByExpr(normData[firstCol:ncol(normData)])
normData <- normData[filtered,]

if(normMethod == 'rpkm'){
  lengths <- getGeneLengths(normData[1])
  normData <- checkGeneLengths(normData, lengths)
}

data <- DGEList(counts = normData[,firstCol:ncol(normData)], genes = normData[,1:firstCol - 1])
data <- calcNormFactors(data, method = "TMM")

# running the normalization
if(normMethod == "tmm"){
  normCounts <- cpm(data)
} else if(normMethod == "rpkm"){
  normCounts <- rpkm(data, gene.length = lengths$length)
}

# applying the new data to the original chart, saving to text file
normData[,firstCol:ncol(normData)] <- normCounts
write.table(normData, file = outFilePath, sep = "\t", row.names = FALSE)