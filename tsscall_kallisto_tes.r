#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
print(args)
# args to give
# 1) Kallisto counts file (will use the average across all samples to call TSSs)
# 2) Final TSS call output
# 3) output prefix
# 4) Filtered GTF


TPMfilter <- 0

#read in basic gtf
   gtf <- read.table(args[4], sep = "\t",stringsAsFactors = F, quote="")
#pull out transcript IDs and append to gtf as 10th column
   txIDcol <- grep("transcript_id",unlist(strsplit(head(gtf[gtf$V3 != "gene",9],1),";")))
   tx_ids <- data.frame(txID=apply(gtf[gtf$V3 != "gene",],1,function(x) gsub("transcript_id ","",unlist(strsplit(as.character(x[9]),split=c(";")))[txIDcol])))
   tx_ids$txID <- gsub("\"","",gsub(" \"","",as.character(tx_ids$txID)))
   gtf$txID <- as.character(rep(NA,nrow(gtf)))
   gtf$txID[gtf$V3 != "gene"] <- as.character(gsub(" ","",gsub("\"","",tx_ids$txID)))
#pull out gene IDs and append to gtf as 10th column, if no gene ids in your gtf then this will just come along for the ride
   gene_ids <- data.frame(geneID=apply(gtf,1,function(x) gsub("gene_id ","",unlist(strsplit(as.character(x[9]),split=c(";")))[1])))
   gene_ids$geneID <- gsub("\"","",gsub(" \"","",as.character(gene_ids$geneID)))
   gtf$geneID <- gene_ids$geneID
 

#read in tsscall file
tsscall <- read.table(args[2],sep="\t",header=T,stringsAsFactors = F)
rownames(tsscall) <- tsscall$TSS.ID

s <- strsplit(tsscall$Transcripts, split = "[;]")
tss_id_and_tx_ids <- data.frame(TxID = unlist(s),TSS.ID = rep(tsscall$TSS.ID, sapply(s, length)))



#read in kallisto counts file
dat <- read.table(args[1],sep="\t",header=T,stringsAsFactors = F)


#Add in TSS and TES columns (flipping strand as approppiate)
dat$TSS <- dat$Start
dat$TSS[dat$Strand=="-"] <- dat$Stop[dat$Strand=="-"]


dat$TES <- dat$Start
dat$TES[dat$Strand=="+"] <- dat$Stop[dat$Strand=="+"]


dat_tsscall <- merge(dat,tss_id_and_tx_ids,by=1)

#get avg counts per transcript file, create columns for aggregating overlapping TSS or TES for each gene (chr;tss;geneID)
outcountsfile <- data.frame(aTxID=dat_tsscall$aTxID,
                            Chr=dat_tsscall$Chr,
                            StartCor=dat_tsscall$Start,
                            StopCor=dat_tsscall$Stop,
                            TSS=dat_tsscall$TSS,
                            TES=dat_tsscall$TES,
                            Strand=dat_tsscall$Strand,
                            GeneID=dat_tsscall$GeneID,
                            GeneName=dat_tsscall$GeneName,
                            GeneType=dat_tsscall$GeneType,
                            AvgTPM=dat_tsscall$AvgTPM,
                            chrTSS=paste0(dat_tsscall$Chr,";",dat_tsscall$TSS,";",dat_tsscall$TSS.ID),
                            chrTES=paste0(dat_tsscall$Chr,";",dat_tsscall$TES,";",dat_tsscall$TSS.ID),
                            TSS.ID=dat_tsscall$TSS.ID)


#for TESs (same gene) at same site, sum TPMs and take the one with highest TPM as "representative transcript"
outcountsfileTES_TPMsum <- merge(aggregate(as.character(outcountsfile$aTxID),by=list(TES_id=outcountsfile$chrTES),toString),aggregate(outcountsfile$AvgTPM,by=list(TES_id=outcountsfile$chrTES),sum),by=1)
colnames(outcountsfileTES_TPMsum)[2:3] <- c("aTxIDs","TES_TPM")

#for each TES take the one with the highest TPM (as the representative txID)
rep.TES <- outcountsfile[order(outcountsfile$AvgTPM,decreasing=T),]
rep.TES <- rep.TES[!duplicated(rep.TES$chrTES),]
colnames(rep.TES)[1] <- "aTxID_TES_rep"

TES_summed <- merge(rep.TES[,c(13,1:10,14)],outcountsfileTES_TPMsum,by=1)
TES_summed$aTxIDs <- gsub(" ","",TES_summed$aTxIDs)


#then for each gene, make TES clusters

#to do this write a function to cluster TESs within 100bp (note that if there are multiple closely spaced TES then a TES cluster could extend > 100bp)
cluster_TES <- function(tss){
   tss.dat <- TES_summed[TES_summed$TSS.ID == tss,]
   tss.dat <- tss.dat[order(tss.dat$TES),]
   tss.TES <- tss.dat$TES
   if (length(tss.TES) > 1){
  
  cluster_filt <- abs(tss.TES[-length(tss.TES)] - tss.TES[-1]) < 100
  cluster_num <- 1
  cluster_status <- cluster_filt[1]
  clusters <- rep(paste0(tss,"_","Cluster1"),length(tss.TES))
  for (i in 1:(length(tss.TES)-1)){
    if(cluster_filt[i] == TRUE){
      clusters[i+1] <- paste0(tss,"_","Cluster",cluster_num)
    } else {
      cluster_num <- cluster_num +1
      clusters[i+1] <- paste0(tss,"_","Cluster",cluster_num)
    }
  }
} else {clusters <- paste0(tss,"_","Cluster1") }

tss.dat$TES_cluster <- clusters

#for each TES cluster, get the transcript IDs as a comma separated list
clusterTPMsIDs <- merge(aggregate(tss.dat$TES_TPM,by=list(TES_cluster=clusters),sum),aggregate(tss.dat$aTxIDs,by=list(TES_cluster=clusters),toString),by=1)
colnames(clusterTPMsIDs)[2:3] <- c("ClusterTPM","ClusterTxIDs")
clusterTPMsIDs$aTxIClusterTxIDs <- gsub(" ","",clusterTPMsIDs$ClusterTxIDs)

tss.cluster.dat <- merge(tss.dat,clusterTPMsIDs,by.x=15,by.y=1)
tss.cluster.dat <- tss.cluster.dat[order(tss.cluster.dat$TES_TPM,decreasing=T),]
tss.cluster.dat <- tss.cluster.dat[!duplicated(tss.cluster.dat$TES_cluster),-c(14:15,17)]
tss.cluster.dat$Percent_TES_TPMs <- round(100*tss.cluster.dat$ClusterTPM / sum(tss.cluster.dat$ClusterTPM),3)

return(tss.cluster.dat)    
}

#create empty dataframe to put in TES clusters
tsscall.with.tes.temp.ordered <- data.frame(TES_cluster=factor(),
                             chrTES=character(),
                             aTxID_TES_rep=character(),
                             Chr=character(),
                             StartCor=numeric(),
                             StopCor=numeric(),
                             TSS=numeric(),
                             TES=numeric(),
                             Strand=character(),
                             GeneID=character(),
                             GeneName=character(),
                             GeneType=character(),
                             TSS.ID=character(),
                             ClusterTPM=numeric(),
                             aTxIClusterTxIDs=character(),
                             Percent_TES_TPMs=numeric())

#run loop for through unique geneIDs, using the function I wrote above to cluster TESs
for (tssID in unique(TES_summed$TSS.ID)){
   tsscall.with.tes.temp.ordered <- rbind(tsscall.with.tes.temp.ordered,cluster_TES(tssID))
}
colnames(tsscall.with.tes.temp.ordered)[3] <- "TES_dom_tx_id"
colnames(tsscall.with.tes.temp.ordered)[15] <- "TES_all_tx_ids"
colnames(tsscall.with.tes.temp.ordered)[10] <- "TES_dom_geneID"

#calculate TSS to TES distance
tsscall.with.tes.temp.ordered$Distance_TSS_to_TES <- tsscall.with.tes.temp.ordered$TES - tsscall.with.tes.temp.ordered$TSS
tsscall.with.tes.temp.ordered$Distance_TSS_to_TES[tsscall.with.tes.temp.ordered$Strand=="-"] <- tsscall.with.tes.temp.ordered$TSS[tsscall.with.tes.temp.ordered$Strand=="-"] - tsscall.with.tes.temp.ordered$TES[tsscall.with.tes.temp.ordered$Strand=="-"]

#Remove genes with TSS downstream of the TES
Cluster_dist_positive <- tsscall.with.tes.temp.ordered[tsscall.with.tes.temp.ordered$Distance_TSS_to_TES > 0,]


#order clustered file by TSS.ID, TES TPM, and gene size
Cluster_ordered <- Cluster_dist_positive[order(Cluster_dist_positive$TSS.ID,Cluster_dist_positive$ClusterTPM,Cluster_dist_positive$Distance_TSS_to_TES,decreasing=T),]


#write all TES clusters to file
write.table(Cluster_ordered,paste0(args[3],"_TES_cluster_counts.txt"),sep="\t",quote=F,row.names=F)


#then order the clusters by clusterTPM and then use duplicated() to select the dominant TES for each gene
DominantTES <- Cluster_ordered[!duplicated(Cluster_ordered$TSS.ID),]


         
           
#merge with TSScall file, a handful of weird conflict driven dup gene ids, sort by TSScall reads and take the strongest TSS for duplicated gene IDs 
tsscall.with.tes.temp <- merge(tsscall,DominantTES[,c(13,10,8,3,14:17)],by.x=1,by.y=1)           
tsscall.with.tes.temp.ordered <- tsscall.with.tes.temp[order(tsscall.with.tes.temp$Reads,decreasing=T),]

#remove genes with TSS downstream of TES
tsscall.with.tes.temp.ordered$Distance_TSS_to_TES <- tsscall.with.tes.temp.ordered$TES - tsscall.with.tes.temp.ordered$Position
tsscall.with.tes.temp.ordered$Distance_TSS_to_TES[tsscall.with.tes.temp.ordered$Strand=="-"] <- tsscall.with.tes.temp.ordered$Position[tsscall.with.tes.temp.ordered$Strand=="-"] - tsscall.with.tes.temp.ordered$TES[tsscall.with.tes.temp.ordered$Strand=="-"]

#Remove genes with TSS downstream of the TES
tsscall.with.tes.ordered.pos <- tsscall.with.tes.temp.ordered[tsscall.with.tes.temp.ordered$Distance_TSS_to_TES > 0,]





tsscall.with.tes <- tsscall.with.tes.ordered.pos[!duplicated(tsscall.with.tes.ordered.pos$TES_dom_geneID),]




#write merged file to test
write.table(tsscall.with.tes,paste0(args[3],"_FINAL_Dominant_obsTSS_and_TES_perGeneID_List.txt"),sep="\t",quote=F,row.names=F )      
                                
#save weird conflict cases to separate case
weird.condlict.geneids <- tsscall.with.tes.temp.ordered$TES_dom_geneID[duplicated(tsscall.with.tes.temp.ordered$TES_dom_geneID)]

#write merged file to test
write.table(tsscall.with.tes.temp.ordered[tsscall.with.tes.temp.ordered$TES_dom_geneID %in% weird.condlict.geneids,],paste0(args[3],"_Weird_conflict_duplicates_obsTSS_and_TES_perGeneID_List.txt"),sep="\t",quote=F,row.names=F )      
 
     

