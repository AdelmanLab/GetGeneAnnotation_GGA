#!/usr/bin/env Rscript


args = commandArgs(trailingOnly=TRUE)
print(args)
# args to give
# 1) Kallisto TSS/TES calls (use the filtered one here)
# 2) TSScall final output
# 3) Kallisto Tx counts all
# 4) output prefix
# 5) TPM threshold


#read in Kallisto tx counts, extract gene info
tx.info <- read.table(args[3],sep="\t",header=T,as.is=T)
gene.info <- tx.info[!duplicated(tx.info$GeneID),c(6:8,2,5)]
                   
#read in kallisto called dominant TSS/TESs, append kallisto. in front of colnames
kal <- read.table(args[1],sep="\t",header=T,as.is=T)[,c(1:3,4:5,9,10,6,8,7)]
colnames(kal) <- paste0("kallisto.",colnames(kal))

#read in TSScall called dominant TSSs (with kal called TESs), append tsscall in front of colnames
tsscall <- read.table(args[2],sep="\t",header=T,as.is=T)
colnames(tsscall) <- paste0("tsscall.",colnames(tsscall))

#merge kallisto and TSScall calls by TES_dom_geneID 
tss.call.gene.id.col <- grep("TES_dom_geneID",colnames(tsscall))
kal.tsscall <- merge(kal,tsscall,by.x=1,by.y=tss.call.gene.id.col,all=T)
colnames(kal.tsscall)[1] <- "GeneID"

#merge with gene info
geneinfo.kal.tsscall <- merge(gene.info,kal.tsscall[,c(1,5:7,9:10,11:ncol(kal.tsscall))],by=1)

#read in TSScall threshold
tpm <- as.numeric(args[5])

#get distance from TSScall to Kallisto called dominant TSSs
geneinfo.kal.tsscall$Distance.kallisto.to.tsscall.TSS <- geneinfo.kal.tsscall$kallisto.TSS - geneinfo.kal.tsscall$tsscall.Position
geneinfo.kal.tsscall$Distance.kallisto.to.tsscall.TSS[geneinfo.kal.tsscall$Strand=="-"] <- (geneinfo.kal.tsscall$kallisto.TSS - geneinfo.kal.tsscall$tsscall.Position)[geneinfo.kal.tsscall$Strand=="-"]

#caall whether TSSs are in kallisto, TSScall or both
#for TSScall only, use the TES cluster TPM to filter, this means that we're using the dominant TSS-TES associated transcripts
tsscall.only.filt <- geneinfo.kal.tsscall$tsscall.ClusterTPM <= tpm & !is.na(geneinfo.kal.tsscall$tsscall.ClusterTPM)
kal.only.filt <- is.na(geneinfo.kal.tsscall$tsscall.Position)
both.filt <- geneinfo.kal.tsscall$tsscall.ClusterTPM > tpm & !is.na(geneinfo.kal.tsscall$tsscall.ClusterTPM)

geneinfo.kal.tsscall$TSScallOnly <- rep(0,nrow(geneinfo.kal.tsscall))
geneinfo.kal.tsscall$TSScallOnly[tsscall.only.filt] <- rep(1,sum(tsscall.only.filt))

geneinfo.kal.tsscall$KallistoOnly <- rep(0,nrow(geneinfo.kal.tsscall))
geneinfo.kal.tsscall$KallistoOnly[kal.only.filt] <- rep(1,sum(kal.only.filt))

geneinfo.kal.tsscall$Both <- rep(0,nrow(geneinfo.kal.tsscall))
geneinfo.kal.tsscall$Both[both.filt] <- rep(1,sum(both.filt))

tsscallMergeOutput <- geneinfo.kal.tsscall[!is.na(geneinfo.kal.tsscall$tsscall.Position),]

#output master speadsheet for all TSScalls

write.table(tsscallMergeOutput,paste0(args[4],"_Merged.kallisto.TSScall.TSS.TES.calls.txt"),sep="\t",quote=F,row.names=F)


#output simplified speadsheet for all TSScalls

simple_df <- data.frame(GeneID=tsscallMergeOutput$GeneID,
                        GeneName=tsscallMergeOutput$GeneName,
                        GeneType=tsscallMergeOutput$GeneType,
                        TSS_ID=tsscallMergeOutput$tsscall.TSS.ID,
                        Chr=tsscallMergeOutput$Chr,
                        Strand=tsscallMergeOutput$Strand,
                        TSScallOnlyUnstableTx=tsscallMergeOutput$TSScallOnly,
                        DominantTSS=tsscallMergeOutput$tsscall.Position,
                        DominantTES=tsscallMergeOutput$tsscall.TES,
                        TSScallReads=tsscallMergeOutput$tsscall.Reads,
                        Distance.kallisto.to.TSScall.TSS=tsscallMergeOutput$Distance.kallisto.to.tsscall.TSS,
                        Divergent=tsscallMergeOutput$tsscall.Divergent.,
                        Divergent.partner=tsscallMergeOutput$tsscall.Divergent.partner,
                        Divergent.distance=tsscallMergeOutput$tsscall.Divergent.distance,
                        Convergent=tsscallMergeOutput$tsscall.Convergent.,
                        Convergent.partner=tsscallMergeOutput$tsscall.Convergent.partner,
                        Convergent.distance=tsscallMergeOutput$tsscall.Convergent.distance,
                        MultipleGeneID=tsscallMergeOutput$tsscall.Multiple_Gene_ID,
                        AlternativeGeneID=tsscallMergeOutput$tsscall.Alternative_Gene.ID,
                        TESdominantTxID=tsscallMergeOutput$tsscall.TES_dom_tx_id,
                        TESclusterTPM=tsscallMergeOutput$tsscall.ClusterTPM,
                        TESallTxIDs=tsscallMergeOutput$tsscall.TES_all_tx_ids)


#write final matrix for just stable obs dominant TSSs
write.table(simple_df[simple_df$TSScallOnlyUnstableTx==0,],paste0(args[4],"_Dominant.TSS.TES.calls_WORKING_FILE.txt"),sep="\t",quote=F,row.names=F)

#Also output matrix including the unstable TSSs
write.table(simple_df,paste0(args[4],"_Dominant.TSS.TES.calls_WORKING_FILE.including.unstable.txt"),sep="\t",quote=F,row.names=F)



#### make BED and makeheatmap files


# write dominant TSS to TES for TESs passing or failing TPM filter
tss_tes_bed <- data.frame(Chr=simple_df$Chr,
                          Start=simple_df$DominantTSS,
                          Stop=simple_df$DominantTES,
                          GeneID=simple_df$GeneID,
                          Score=rep(0,nrow(simple_df)),
                          Strand=simple_df$Strand)
tss_tes_bed$Start[simple_df$Strand=="-"] <- simple_df$DominantTES[simple_df$Strand=="-"]
tss_tes_bed$Stop[simple_df$Strand=="-"] <- simple_df$DominantTSS[simple_df$Strand=="-"]

#make BED file 0-based
tss_tes_bed$Start <- tss_tes_bed$Start -1

tss_tes_makeheatmap <- data.frame(Name=tss_tes_bed$GeneID,
                                  anchor=simple_df$DominantTSS,
                                  Chr=simple_df$Chr,
                                  Start=tss_tes_bed$Start,
                                  Stop=tss_tes_bed$Stop,
                                  Strand=tss_tes_bed$Strand)


#write bed/makeheatmap for stable only files
write.table(tss_tes_bed[simple_df$TSScallOnlyUnstableTx==0,],paste0(args[4],"_Dominant.TSS.TES.calls.bed"),sep="\t",quote=F,row.names=F,col.names=F)      
write.table(tss_tes_makeheatmap[simple_df$TSScallOnlyUnstableTx==0,],paste0(args[4],"_Dominant.TSS.TES.calls_formakeheatmap.txt"),sep="\t",quote=F,row.names=F,col.names=F)      

#write bed/makeheatmap for files including unstable
write.table(tss_tes_bed,paste0(args[4],"_Dominant.TSS.TES.calls.including.unstable.bed"),sep="\t",quote=F,row.names=F,col.names=F)      
write.table(tss_tes_makeheatmap,paste0(args[4],"_Dominant.TSS.TES.calls_formakeheatmap.including.unstable.txt"),sep="\t",quote=F,row.names=F,col.names=F)      

#write bed/makeheatmap for unstable only
write.table(tss_tes_bed[simple_df$TSScallOnlyUnstableTx==1,],paste0(args[4],"_Dominant.TSS.calls.toLast.anno.TES.only.unstable.bed"),sep="\t",quote=F,row.names=F,col.names=F)      
write.table(tss_tes_makeheatmap[simple_df$TSScallOnlyUnstableTx==1,],paste0(args[4],"_Dominant.TSS.calls.toLast.anno.TES_formakeheatmap.only.unstable.txt"),sep="\t",quote=F,row.names=F,col.names=F)      

                   


#Also output the kallisto only genes

conflict.genes <- unique(c(unlist(strsplit(tsscall$tsscall.Alternative_Gene.ID[tsscall$tsscall.Conflict_GeneID_Yes1==1],";")),
                    tsscall$tsscall.Unique_GeneID[tsscall$tsscall.Conflict_GeneID_Yes1==1]))
conflict.genes.in.kal.only <- kal.only.filt & !(geneinfo.kal.tsscall$GeneID %in% conflict.genes)

write.table(geneinfo.kal.tsscall[conflict.genes.in.kal.only,],paste0(args[4],"_Merged.KALLISTO.ONLY.TSS.TES.calls.txt"),sep="\t",quote=F,row.names=F)



