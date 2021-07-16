#!/usr/bin/env Rscript


args = as.character(commandArgs(trailingOnly=TRUE))
print(args)
# args to give
# 1) Dom TSS/TES matrix file
# 2) basic GTF annotations
# 3) output name for GTF

#read in dominant TSS-TES matrix
dat <- read.table(args[1],sep="\t", stringsAsFactors = F, header=T)
#read in basic GTF
gtfOrig <- read.table(args[2],sep="\t", stringsAsFactors = F, quote="")
#output gtf filename
outputFileName <- args[3]


#get list of transcript IDs to keep, with/without filtering for those associated with stable transcripts
tx2keep_all <- unlist(apply(dat,1,function(x) strsplit(as.character(x[22]),",") ))
gene2keep_all <- dat$GeneID

#note, this step is redundant because the matrix is now only containing stable transcripts.
tx2keep_stable <- unlist(apply(dat[dat$TSScallOnlyUnstableTx==0,],1,function(x) strsplit(as.character(x[22]),",") ))
gene2keep_stable <- dat$GeneID[dat$TSScallOnlyUnstableTx==0]



#pull out tx ids from gtf and write transcripts passing filters to file
tx_id <- sapply(strsplit(sapply(strsplit(sapply(strsplit(gtfOrig$V9, "transcript_id"), "[", 2), ";"), "[", 1), "\""), "[", 2)
gene_id <- sapply(strsplit(sapply(strsplit(sapply(strsplit(gtfOrig$V9, ";"), "[", 1), " "), "[", 2), "\""), "[", 2)




# write GTFs to file (only output the stable transcripts)

#write.table(gtfOrig[tx_id %in% tx2keep_all, ],paste0(outputPrefix,".annotated.transcripts.gtf"), quote = F, sep = "\t", row.names = F, col.names = F)

write.table(gtfOrig[tx_id %in% tx2keep_stable & gene_id %in% gene2keep_stable,],
            outputFileName, quote = F, sep = "\t", row.names = F, col.names = F)

#summary(tx_id %in% tx2keep_stable & gtfOrig$V3=="transcript")
#summary(tx_id %in% tx2keep_stable & gtfOrig$V3=="transcript" & gene_id %in% gene2keep_stable)

#print messsages with transcript numbers

#print(paste0(sum(gtfOrig$V3=="transcript"), " Trancsripts in basic GTF"))

#print(paste0(nrow(dat), " Dom TSS in basic GTF"))
#print(paste0(sum(tx_id %in% tx2keep_all & gtfOrig$V3=="transcript"), " Trancsripts stemming from dominant TSS and ending in dominant TES"))

#print(paste0(nrow(dat[stable_tx_tss,]), " Stable transcript associated Dom TSS in basic GTF"))
#print(paste0(sum(tx_id %in% tx2keep_stable & gtfOrig$V3=="transcript" & 
#                  gene_id %in% gene2keep_stable), " Stable trancsripts stemming from dominant TSS and ending in dominant TES"))







   
