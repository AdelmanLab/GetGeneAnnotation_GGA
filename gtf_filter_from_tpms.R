#!/usr/bin/env Rscript

args = as.character(commandArgs(trailingOnly=TRUE))
print(args)
# args to give
# 1) TPM threshold
# 2) number of samples needing to pass filter
# 3) GTF to filter
# 4) output name for gtf
# 5) transcript info file
# 6) clean Chr file


tpmThresh <- as.numeric(args[1])
sampThresh <- as.numeric(args[2])
gtfOrig <- read.table(args[3],sep="\t", stringsAsFactors = F, quote="")
outputPrefix <- args[4]
txinfo <- read.table(args[5],sep="\t", stringsAsFactors = F, header=F)[,1:8]
colnames(txinfo) <- c("aTxID","Chr","Start","Stop","Strand","GeneID","GeneName","GeneType")
#note, added an a to "aTxID" to help with sorting later (so that header stays at top)

#read in clean Chr list
cleanChr <- read.table(args[6],sep="\t",header=F)$V1

#read in abundFiles, make TPM matrix, and call transcripts passing filter
abunFiles <- dir(Sys.glob(file.path(paste0(args[4],"/*"))), pattern = "abundance.tsv", full.names = T)
abunList <- vector("list", length(abunFiles))
for (i in 1:length(abunFiles)){
  abunList[[i]] <- read.table(abunFiles[i], header = T, sep = "\t", row.names = 1)
}
names(abunList) <- sapply(strsplit(abunFiles,"/"), "[", 1)
tpmMat <- matrix(unlist(sapply(abunList, "[", "tpm")), nrow = nrow(abunList[[1]]), ncol = length(abunList), 
                 dimnames = list(rownames(abunList[[1]]), names(abunList)))


#then merge txinfo with tpmMat

tpmAvg <- data.frame(aTxID=rownames(tpmMat),AvgTPM=rowMeans(tpmMat))
tpmMatInfo <- merge(txinfo,tpmAvg)

#remove weird Chrs
tpmMatInfo <- tpmMatInfo[tpmMatInfo$Chr %in% cleanChr,]

write.table(tpmMatInfo, paste0(outputPrefix,"_avg_tpms_all_tx.txt"), quote = F, sep = "\t", row.names = F, col.names = T)

#pull out tx ids from gtf and write transcripts passing filters to file
tx2keep <- rownames(tpmMat)[rowSums(tpmMat > tpmThresh) >= sampThresh & rownames(tpmMat) %in% tpmMatInfo$aTxID]

tx_id <- sapply(strsplit(sapply(strsplit(sapply(strsplit(gtfOrig$V9, "transcript_id"), "[", 2), ";"), "[", 1), "\""), "[", 2)
write.table(gtfOrig[tx_id %in% tx2keep, ], paste0(outputPrefix,"_kallisto-filt.gtf"), quote = F, sep = "\t", row.names = F, col.names = F)

write.table(tpmMatInfo[tpmMatInfo$aTxID %in% tx2keep,], paste0(outputPrefix,"_avg_tpms_kallisto-filt-tx.txt"), quote = F, sep = "\t", row.names = F, col.names = T)





   
