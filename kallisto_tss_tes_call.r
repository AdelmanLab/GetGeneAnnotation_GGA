#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
print(args)
# args to give
# 1) Kallisto counts file (use the filtered one herer)
# 2) filtered GTF
# 3) output prefix
# 4) Chrs to keep

#args <- c("mESC_RNAseq_test_avg_tpms_kallisto-filt-tx.txt","mESC_RNAseq_test_kallisto-filt.gtf","test_output")
                   
#read in clean chromosomes
cleanChr <- read.table(args[4],sep="\t",header=F)$V1

#read in filtered gtf
   gtf <- read.table(args[2], sep = "\t",stringsAsFactors = F, quote="")
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
 

#read in kallisto counts file
dat <- read.table(args[1],sep="\t",header=T,stringsAsFactors = F)

#Add in TSS and TES columns (flipping strand as approppiate)
dat$TSS <- dat$Start
dat$TSS[dat$Strand=="-"] <- dat$Stop[dat$Strand=="-"]

dat$TES <- dat$Start
dat$TES[dat$Strand=="+"] <- dat$Stop[dat$Strand=="+"]

#get avg counts per transcript file, create columns for aggregating overlapping TSS or TES for each gene (chr;tss;geneID)
outcountsfile <- data.frame(aTxID=dat$aTxID,
                            Chr=dat$Chr,
                            StartCor=dat$Start,
                            StopCor=dat$Stop,
                            TSS=dat$TSS,
                            TES=dat$TES,
                            Strand=dat$Strand,
                            GeneID=dat$GeneID,
                            GeneName=dat$GeneName,
                            GeneType=dat$GeneType,
                            AvgTPM=dat$AvgTPM,
                            chrTSS=paste0(dat$Chr,";",dat$TSS,";",dat$GeneID),
                            chrTES=paste0(dat$Chr,";",dat$TES,";",dat$GeneID))

#remove transcripts with a TPM of 0
outcountsfile <- outcountsfile[outcountsfile$AvgTPM > 0,]

#for TSSs (same gene) at same site, sum TPMs and take the one with highest TPM as "representative transcript"
outcountsfileTSS_TPMsum <- merge(aggregate(as.character(outcountsfile$aTxID),by=list(TSS_id=outcountsfile$chrTSS),toString),aggregate(outcountsfile$AvgTPM,by=list(TSS_id=outcountsfile$chrTSS),sum),by=1)
colnames(outcountsfileTSS_TPMsum)[2:3] <- c("aTxIDs","TSS_TPM")

#for each TSS take the one with the highest TPM (as the representative txID)
rep.TSS <- outcountsfile[order(outcountsfile$AvgTPM,decreasing=T),]
rep.TSS <- rep.TSS[!duplicated(rep.TSS$chrTSS),]
colnames(rep.TSS)[1] <- "aTxID_TSS_rep"

TSS_summed <- merge(rep.TSS[,c(12,1:10,13)],outcountsfileTSS_TPMsum,by=1)
TSS_summed$aTxIDs <- gsub(" ","",TSS_summed$aTxIDs)


#then for each gene, make TSS clusters

#to do this write a function to cluster TSSs within 100bp (note that if there are multiple closely spaced TSS then a TSS cluster could extend > 100bp)
cluster_tss <- function(gene){
   gene.dat <- TSS_summed[TSS_summed$GeneID == gene,]
   gene.dat <- gene.dat[order(gene.dat$TSS),]
   gene.tss <- gene.dat$TSS
   if (length(gene.tss) > 1){
  
  cluster_filt <- abs(gene.tss[-length(gene.tss)] - gene.tss[-1]) < 100
  cluster_num <- 1
  cluster_status <- cluster_filt[1]
  clusters <- rep(paste0(gene,"_","Cluster1"),length(gene.tss))
  for (i in 1:(length(gene.tss)-1)){
    if(cluster_filt[i] == TRUE){
      clusters[i+1] <- paste0(gene,"_","Cluster",cluster_num)
    } else {
      cluster_num <- cluster_num +1
      clusters[i+1] <- paste0(gene,"_","Cluster",cluster_num)
    }
  }
} else {clusters <- paste0(gene,"_","Cluster1") }

gene.dat$TSS_cluster <- clusters

#for each TSS cluster, get the transcript IDs as a comma separated list
clusterTPMsIDs <- merge(aggregate(gene.dat$TSS_TPM,by=list(TSS_cluster=clusters),sum),aggregate(gene.dat$aTxIDs,by=list(TSS_cluster=clusters),toString),by=1)
colnames(clusterTPMsIDs)[2:3] <- c("ClusterTPM","ClusterTxIDs")
clusterTPMsIDs$aTxIClusterTxIDs <- gsub(" ","",clusterTPMsIDs$ClusterTxIDs)

gene.cluster.dat <- merge(gene.dat,clusterTPMsIDs,by.x=15,by.y=1)
gene.cluster.dat <- gene.cluster.dat[order(gene.cluster.dat$TSS_TPM,decreasing=T),]
gene.cluster.dat <- gene.cluster.dat[!duplicated(gene.cluster.dat$TSS_cluster),-c(14:15,17)]

return(gene.cluster.dat)    
}

#create empty dataframe to put in TSS clusters
Cluster_summed <- data.frame(TSS_cluster=factor(),chrTSS=character(),
                             aTxID_TSS_rep=character(),
                             Chr=character(),
                             StartCor=numeric(),
                             StopCor=numeric(),
                             TSS=numeric(),
                             TES=numeric(),
                             Strand=character(),
                             GeneID=character(),
                             GeneName=character(),
                             GeneType=character(),
                             chrTES=character(),
                             ClusterTPM=numeric(),
                             aTxIClusterTxIDs=character())

#run loop for through unique geneIDs, using the function I wrote above to cluster TSSs
for (geneID in unique(TSS_summed$GeneID)){
   Cluster_summed <- rbind(Cluster_summed,cluster_tss(geneID))
}

#then order the clusters by clusterTPM and then use duplicated() to select the dominant TSS for each gene
Cluster_ordered <- Cluster_summed[order(Cluster_summed$ClusterTPM,decreasing=T),]
DominantTSS <- Cluster_ordered[!duplicated(Cluster_ordered$GeneID),]


##################################### TES section

#now repeat this for dominant TES

s <- strsplit(DominantTSS$aTxIClusterTxIDs, split = "[,]")
dom_tss_txs <- unlist(s)

#remove transcripts not from dominant TSS
outcountsfile_tes <- outcountsfile[outcountsfile$aTxID %in% dom_tss_txs,]


#for TESs (same gene) at same site, sum TPMs and take the one with highest TPM as "representative transcript"
outcountsfileTES_TPMsum <- merge(aggregate(as.character(outcountsfile_tes$aTxID),by=list(TES_id=outcountsfile_tes$chrTES),toString),aggregate(outcountsfile_tes$AvgTPM,by=list(TES_id=outcountsfile_tes$chrTES),sum),by=1)
colnames(outcountsfileTES_TPMsum)[2:3] <- c("aTxIDs","TES_TPM")


#for each TES take the one with the highest TPM (as the representative txID)
rep.TES <- outcountsfile_tes[order(outcountsfile_tes$AvgTPM,decreasing=T),]
rep.TES <- rep.TES[!duplicated(rep.TES$chrTES),]
colnames(rep.TES)[1] <- "aTxID_TES_rep"

TES_summed <- merge(rep.TES[,c(13,1:10,13)],outcountsfileTES_TPMsum,by=1)
TES_summed$aTxIDs <- gsub(" ","",TES_summed$aTxIDs)


#then for each gene, make TES clusters

#to do this write a function to cluster TESs within 100bp (note that if there are multiple closely spaced TES then a TES cluster could extend > 100bp)
cluster_TES <- function(gene){
   gene.dat <- TES_summed[TES_summed$GeneID == gene,]
   gene.dat <- gene.dat[order(gene.dat$TES),]
   gene.TES <- gene.dat$TES
   if (length(gene.TES) > 1){
  
  cluster_filt <- abs(gene.TES[-length(gene.TES)] - gene.TES[-1]) < 100
  cluster_num <- 1
  cluster_status <- cluster_filt[1]
  clusters <- rep(paste0(gene,"_","Cluster1"),length(gene.TES))
  for (i in 1:(length(gene.TES)-1)){
    if(cluster_filt[i] == TRUE){
      clusters[i+1] <- paste0(gene,"_","Cluster",cluster_num)
    } else {
      cluster_num <- cluster_num +1
      clusters[i+1] <- paste0(gene,"_","Cluster",cluster_num)
    }
  }
} else {clusters <- paste0(gene,"_","Cluster1") }

gene.dat$TES_cluster <- clusters

#for each TES cluster, get the transcript IDs as a comma separated list
clusterTPMsIDs <- merge(aggregate(gene.dat$TES_TPM,by=list(TES_cluster=clusters),sum),aggregate(gene.dat$aTxIDs,by=list(TES_cluster=clusters),toString),by=1)
colnames(clusterTPMsIDs)[2:3] <- c("ClusterTPM","ClusterTxIDs")
clusterTPMsIDs$aTxIClusterTxIDs <- gsub(" ","",clusterTPMsIDs$ClusterTxIDs)

gene.cluster.dat <- merge(gene.dat,clusterTPMsIDs,by.x=15,by.y=1)
gene.cluster.dat <- gene.cluster.dat[order(gene.cluster.dat$TES_TPM,decreasing=T),]
gene.cluster.dat <- gene.cluster.dat[!duplicated(gene.cluster.dat$TES_cluster),-c(14:15,17)]

return(gene.cluster.dat)    
}



#create empty dataframe to put in TES clusters
TES_Cluster_summed <- data.frame(TES_cluster=factor(),chrTES=character(),
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
                             chrTES=character(),
                             ClusterTPM=numeric(),
                             aTxIClusterTxIDs=character())

#run loop for through unique geneIDs, using the function I wrote above to cluster TESs
for (geneID in unique(TES_summed$GeneID)){
   TES_Cluster_summed <- rbind(TES_Cluster_summed,cluster_TES(geneID))
}



## merge tss and tes?
# merge dominant TSS/TES files? I guess I don't really need to, but doing it all the same. 
Dominant_tss_tes_clusters <- merge(DominantTSS[,c(10,11,12,4,7,9,15)],TES_Cluster_summed[,c(10,3,8,14,15)],by=1)
colnames(Dominant_tss_tes_clusters)[c(7:8,11)] <- c("TSS.associated.Tx.IDs","Dom.Tx.ID","TES.associated.Tx.IDs")



#calculate TSS to TES distance
Dominant_tss_tes_clusters$Distance_TSS_to_TES <- Dominant_tss_tes_clusters$TES - Dominant_tss_tes_clusters$TSS
Dominant_tss_tes_clusters$Distance_TSS_to_TES[Dominant_tss_tes_clusters$Strand=="-"] <- Dominant_tss_tes_clusters$TSS[Dominant_tss_tes_clusters$Strand=="-"] - Dominant_tss_tes_clusters$TES[Dominant_tss_tes_clusters$Strand=="-"]

#Remove genes with TSS downstream of the TES (should be 0 here, this is really for the TSScall version)
Cluster_dist_positive <- Dominant_tss_tes_clusters[Dominant_tss_tes_clusters$Distance_TSS_to_TES > 0,]


#order clustered file by TSS.ID, TES TPM, and gene size
Cluster_ordered <- Cluster_dist_positive[order(Cluster_dist_positive$GeneID,Cluster_dist_positive$ClusterTPM,Cluster_dist_positive$Distance_TSS_to_TES,decreasing=T),]


#write all TES clusters to file
write.table(Cluster_ordered,paste0(args[3],"_kallisto_domTSS_TES_cluster_counts.txt"),sep="\t",quote=F,row.names=F)


#then order the clusters by clusterTPM and then use duplicated() to select the dominant TES for each gene
DominantTES <- Cluster_ordered[!duplicated(Cluster_ordered$GeneID),]


# filter to only keep the cleaned up chromosomes
Dominant_tss_tes <- DominantTES[DominantTES$Chr %in% cleanChr,]


Dominant_tss_tes_bed <- data.frame(Chr=Dominant_tss_tes$Chr,
                                   Start=Dominant_tss_tes$TSS,
                                   Stop=Dominant_tss_tes$TES,
                                   Name=Dominant_tss_tes$GeneID,
                                   Score=rep(0,nrow(Dominant_tss_tes)),
                                   Strand=Dominant_tss_tes$Strand)
Dominant_tss_tes_bed$Start[Dominant_tss_tes$Strand=="-"] <- Dominant_tss_tes$TES[Dominant_tss_tes$Strand=="-"]
Dominant_tss_tes_bed$Stop[Dominant_tss_tes$Strand=="-"] <- Dominant_tss_tes$TSS[Dominant_tss_tes$Strand=="-"]


write.table(Dominant_tss_tes,paste0(args[3],"_Kallisto_called_dominant_TSS_to_TES_clusters.txt"),sep="\t",quote=F,row.names=F,col.names=T)      


tss_tes_makeheatmap <- data.frame(Name=Dominant_tss_tes_bed$Name,
                                  anchor=Dominant_tss_tes$TSS,
                                  Chr=Dominant_tss_tes_bed$Chr,
                                  Start=Dominant_tss_tes_bed$Start,
                                  Stop=Dominant_tss_tes_bed$Stop,
                                  Strand=Dominant_tss_tes_bed$Strand)

write.table(Dominant_tss_tes_bed,paste0(args[3],"_Kallisto_called_dominant_TSS_to_TES.bed"),sep="\t",quote=F,row.names=F,col.names=F)      
write.table(tss_tes_makeheatmap,paste0(args[3],"_Kallisto_called_dominant_TSS_to_TES_formakeheatmap.txt"),sep="\t",quote=F,row.names=F,col.names=F)      


                   
                    
#get dominant transcript GTF             
dominant.tx.gtf <- gtf[gtf$txID %in% Dominant_tss_tes$Dom.Tx.ID,]
                              
#write dominant GTF to txt file
write.table(dominant.tx.gtf[,1:9],paste0(args[3],"_Kallisto_called_dominant_transcript.per.gene.gtf"),sep="\t",quote=F,row.names=F,col.names=F)
    
    
#get dominant TSS to dominant TES transcripts
#get list of transcript IDs to keep, with/without filtering for those associated with stable transcripts
tx2keep_all <- unlist(apply(data.frame(TxIDs=Dominant_tss_tes$TES.associated.Tx.IDs),
                      1,function(x) strsplit(as.character(x[1]),",") ))
dominant.clusters.tx.gtf <- gtf[gtf$txID %in% tx2keep_all,]

write.table(dominant.clusters.tx.gtf[,1:9],paste0(args[3],"_Kallisto_called_dominant_TSS_to_TES_cluster.transcripts.gtf"),sep="\t",quote=F,row.names=F,col.names=F)


### retun all non-redundant TSSs with TPM scores

#mESC_test_avg_tpms_all_tx.txt



