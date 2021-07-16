#!/usr/bin/env Rscript

args = as.character(commandArgs(trailingOnly=TRUE))
print(args)
# args to give
# 1) input TSS classify file
# 2) output prefix
# 3) R library directory path
# 4) clean chromosome list
# 5) TPM threshold

#3 functions from CAM, first define these and then run them

#load stringr package
  libDir <- args[3]
  library(stringr,lib.loc=libDir)


# read in clean chromosome list
cleanChr <- read.table(args[4],sep="\t",header=F)$V1


#remove nuTSS annotations 
#Function input Output_TSScall_detail_file_TSSclassify_TSSfilter
rvm_nuTSS <- function(TSScall){
  TSScall <- TSScall[TSScall$Chromosome %in% cleanChr,]
  TSScall$TSS.ID_type <- substr(TSScall$TSS.ID, 1 , 3)            #In the TSS.ID column, TSSs are labeled as obsTSS_#, nuTSS_# and annoTSS_#. This line will extract the first three characters in the TSS.ID Column for each row
  TSScall_obs <- subset(TSScall, TSScall$TSS.ID_type == "obs")     #subset obsTSS
  TSScall_ann <- subset(TSScall, TSScall$TSS.ID_type == "ann")     #subset obsTSS
  TSScall_nuTSS <- subset(TSScall, TSScall$TSS.ID_type == "nuT")
  
  write.table(TSScall_obs, paste0(args[2],"_TSScall_Output_obsTSS_List.txt"), sep="\t", row.names = FALSE, quote = FALSE)     #Export List with only obs 
  write.table(TSScall_ann, paste0(args[2],"_TSScall_Output_annoTSS_List.txt"), sep="\t", row.names = FALSE, quote = FALSE)     #Export List with only annotated TSS. To be used for Dominant calling with Kalisto 
  write.table(TSScall_nuTSS, paste0(args[2],"_TSScall_Output_nuTSS_List_IGNORE.txt"), sep="\t", row.names = FALSE, quote = FALSE)     #Export List with nuTSS
  
  print(nrow(TSScall_obs))
  print(nrow(TSScall_ann))
  print(nrow(TSScall_nuTSS))
  
  return(TSScall_obs)
  
}

### Unique IDs
Unique_ID_Info <- function(TSScall_Filtered){
  
  #Flag TSSs with multiple ENSUMG IDs in the Gene.ID column
  value = c(";") 
  TSScall_Filtered$Multiple_Gene_ID <- grepl(value, TSScall_Filtered$Gene.ID)   # If a semicolon is in the string, TRUE 
  
  Single_ID <- subset(TSScall_Filtered, TSScall_Filtered$Multiple_Gene_ID == FALSE) #Just 1 Gene.ID listed 
  Single_ID$Unique_GeneID <- Single_ID$Gene.ID
  
  loop_temp <- subset(TSScall_Filtered, TSScall_Filtered$Multiple_Gene_ID == TRUE) #TSSs with multiple GeneIDs
  loop_temp$Unique_GeneID <- NA
  i = 1
  

  
  while (i <= nrow(loop_temp)) {
    IDs <- str_split(loop_temp$Gene.ID[i], value, n = Inf, simplify = TRUE) #Break the single string into a matric whos size is dependent on the number of ENSUMGS 
    IDs <- t(IDs) #Transpose the matrix 
    Unique_ID <- unique(IDs)  #Get a list of unique IDs 
    
    #If there is only 1 unique Gene.ID. Then add the Unique ID to the Unique_GeneID Column
    if (isTRUE(nrow(Unique_ID) == 1)){
      loop_temp$Unique_GeneID[i] <- Unique_ID
      
      #Flag it, to be dealt with later 
    } else if (isTRUE(nrow(Unique_ID) > 1)){
      loop_temp$Unique_GeneID[i] <- "TBD"
      
    }
    i = i+1
    
  }
  
  #Merge Annnotations with just 1 Gene.ID listed, and those that had multiple ENSUMG listed that only corresponded to 1 unique Gene ID
  Unique_Filtered <- subset(loop_temp, loop_temp$Unique_GeneID != "TBD")
  Single_ID_V2 <- rbind(Single_ID, Unique_Filtered)
  
  #Clean Up TSSs with multiple unique ENSUMG.IDs. Use Geoff's script to expand the file so its 1 gene per line.  
  tssTable <- subset(loop_temp, loop_temp$Unique_GeneID == "TBD")
  
  genesPerTss <- lapply(strsplit(as.character(tssTable$Gene.ID), ";"), unique)
  expandedTss <- tssTable[0, ]
  expandedTss$Gene.ID <- character()
  
  rowCounter = 1
  for (i in 1:length(genesPerTss)){
    for (j in 1:length(genesPerTss[[i]])){
      expandedTss[rowCounter, ] <- tssTable[i, ]
      expandedTss[rowCounter, "Gene.ID"] <- genesPerTss[[i]][j]
      rowCounter <- rowCounter + 1
    }
  }
  
  #Deduplicate The list by Gene.ID. If a Gene.ID has two associated TSSs, Retain the TSS with the most reads 
  List <- unique(expandedTss$Gene.ID)
  i = 1
  
  while (i <= length(List)) {
    
    Working <- subset(expandedTss, expandedTss$Gene.ID == List[i])
    
    #If its unique, don't do anything 
    if (isTRUE(nrow(Working) > 1)){
      Working <- Working[order(-Working$Reads) , ]
      Working <- Working[1,]
    }
    
    if (isTRUE(i == 1)){
      Conflict_Filtered <- Working
      
    } else {
      Conflict_Filtered <- rbind(Conflict_Filtered, Working)
    }
    i = i + 1
  }
  
  Conflict_Filtered$Conflict_GeneID_Yes1 <- 1
  Conflict_Filtered$Unique_GeneID <- Conflict_Filtered$Gene.ID
  write.table(Conflict_Filtered, paste0(args[2],"_Conflict_Filtered.txt"), sep = "\t", row.names = F, col.names = T, quote = F)
  
  
  #Merge Together Conflict Filtered and SingleID_V2 Lists 
  Single_ID_V2$Conflict_GeneID_Yes1 <- 0 
  Final <- rbind(Single_ID_V2, Conflict_Filtered)
  
  write.table(Final, paste0(args[2],"_TSScall_Output_obsTSS_List_Filtered_For_Dom_TSScalling.txt"), sep = "\t", row.names = F, col.names = T, quote = F) #Output will be used for Dominant TSS calling
  return(Final)
}

#Generate a List of Dominant obsTSSs per Gene Model 
Dominant_TSScall <- function(Unique_ID_Info){
  
  #Generate a unique list of Gene.IDs
  GeneID <- unique(Unique_ID_Info$Unique_GeneID)
  
  #Upload Kallisto Dominant TSS to TES list 
  Kallisto <- read.delim(paste0(args[2], "_Kallisto_called_dominant_TSS_to_TES_clusters.txt"), sep="\t", as.is = T)
  
  #Condition for while loop and files to hold dominant and nondominant TSS lists
  i = 1
  Dominant_TSS_List <- subset(Unique_ID_Info, Unique_ID_Info$TSS.ID == 1)
  Dominant_TSS_List$Percent_TSS_Reads <- logical(0)
  Dominant_TSS_List$Dist_from_Kallisto_TSS <- logical(0)
  Dominant_TSS_List$obsTSS_per_Gene.ID <- logical(0)
  
  Nondom_List <- subset(Unique_ID_Info, Unique_ID_Info$TSS.ID == 1)
  Nondom_List$Percent_TSS_Reads  <- logical(0)
  Nondom_List$obsTSS_per_Gene.ID <- logical(0)
  Nondom_List$Dist_from_Kallisto_TSS <- logical(0) 
  
  #Columns to Add to the Unique_ID_Info_Dataframe 
  Unique_ID_Info$Percent_TSS_Reads <- NA
  Unique_ID_Info$Dist_from_Kallisto_TSS <- NA
  Unique_ID_Info$obsTSS_per_Gene.ID <- NA
  
  while (i <= length(GeneID)) {
    
    #Pull All obsTSS annotations with the same gene ID
    obsTSS_perID <- subset(Unique_ID_Info, Unique_ID_Info$Unique_GeneID == GeneID[i])
    Total_Reads <- sum(obsTSS_perID$Reads)
    obsTSS_perID$Percent_TSS_Reads <- (obsTSS_perID$Reads / Total_Reads) * 100
    
    #If there are multiple obsTSSs per a given Gene.ID, the dominant obsTSS is...
    if (isTRUE(nrow(obsTSS_perID) > 1)){
      obsTSS_perID <- obsTSS_perID[order(-obsTSS_perID$Reads) , ] #Descending order by Read count, TSS with max counts will be row 1
      obsTSS_perID$obsTSS_per_Gene.ID = nrow(obsTSS_perID)  # For record purposes later 
      
      #If the dominant TSS, doesn't have a clear majority over the other obsTSS within the same gene model. We will pick the TSS Closest to the kalisto TSS
      Kallisto_Info <- subset(Kallisto, Kallisto$GeneID == obsTSS_perID$Unique_GeneID[1])
      
        if (isTRUE(obsTSS_perID$Percent_TSS_Reads[1] < 60 & nrow(Kallisto_Info) > 0 )){
          #Calculate the distance from the Kallisto Dominant TSS
          obsTSS_perID$Dist_from_Kallisto_TSS <- abs(obsTSS_perID$Position - Kallisto_Info$TSS)
          
          #If there is a direct Tie in read counts in the first two columns. Call the dominant TSS as the one that's closest to the Kalisto TSS
          if (isTRUE(obsTSS_perID$Reads[1] == obsTSS_perID$Reads[2])){
             obsTSS_perID <- obsTSS_perID[order(-obsTSS_perID$Reads, obsTSS_perID$Dist_from_Kallisto_TSS) , ]
             
             Dominant <- obsTSS_perID[1,]
             Nondominant <- obsTSS_perID[2:nrow(obsTSS_perID),]
             Dominant_TSS_List <- rbind(Dominant_TSS_List, Dominant) #export dominant TSS
             Nondom_List <- rbind(Nondom_List, Nondominant) #export flagged non-dominant TSSs
             
          } else if (isTRUE(obsTSS_perID$Percent_TSS_Reads[1] >= 40 &  obsTSS_perID$Percent_TSS_Reads[2] >= 40)){
            #Subset just the first two annotations. Pick the one thats closest to the the kallisto TSS
            temp <- subset(obsTSS_perID, obsTSS_perID$Percent_TSS_Reads >= 40)
            other <- subset(obsTSS_perID, obsTSS_perID$Percent_TSS_Reads < 40)
            temp <- temp[order(temp$Dist_from_Kallisto_TSS) , ]
            
            Dominant <- temp[1,]
            Nondominant <- rbind(temp[2,],other) 
            Dominant_TSS_List <- rbind(Dominant_TSS_List, Dominant) #export dominant TSS
            Nondom_List <- rbind(Nondom_List, Nondominant) #export flagged non-dominant TSSs
            
          } else {
            #If all the other conditiosn aren't true, call the dominant TSS as the obsTSS with the most reads.
            obsTSS_perID <- obsTSS_perID[order(-obsTSS_perID$Reads) , ] #Descending order by Read count, TSS with max counts will be row 1
            Dominant <- obsTSS_perID[1,]
            Nondominant <- obsTSS_perID[2:nrow(obsTSS_perID),]
            
            Dominant_TSS_List <- rbind(Dominant_TSS_List, Dominant) #export dominant TSS
            Nondom_List <- rbind(Nondom_List, Nondominant) #export flagged non-dominant TSSs
          }
          
          rm(Kallisto_Info)
        
      #If the obsTSS has more than 60% of the Percent TSS Coverage. Call the dominant TSS as the one with the most reads 
      } else {
        obsTSS_perID <- obsTSS_perID[order(-obsTSS_perID$Reads) , ] #Descending order by Read count, TSS with max counts will be row 1
        Dominant <- obsTSS_perID[1,]
        Nondominant <- obsTSS_perID[2:nrow(obsTSS_perID),]
        
        Dominant_TSS_List <- rbind(Dominant_TSS_List, Dominant) #export dominant TSS
        Nondom_List <- rbind(Nondom_List, Nondominant) #export flagged non-dominant TSSs
      }
      
    #If there is only 1 obsTSS for the Gene Model 
    } else {
      obsTSS_perID$obsTSS_per_Gene.ID = 1
      Dominant_TSS_List <- rbind(Dominant_TSS_List, obsTSS_perID) #export dominant TSS
    }
  print(i)
  i = i + 1
  
  }
  
  write.table(Dominant_TSS_List, paste0(args[2],"_Temp_Dominant_obsTSS_perGeneID_List.txt"), sep = "\t", row.names = F, col.names = T, quote = F)
  write.table(Nondom_List, paste0(args[2],"_Flagged_Nondominant_obsTSS_perGeneID_List.txt"), sep = "\t", row.names = F, col.names = T, quote = F)
  
  
  #Remove Position Duplicates caused by gene models with multiple unique Gene.IDs 
  Duplicates <- subset(Dominant_TSS_List, Dominant_TSS_List$Conflict_GeneID_Yes1 == 1 )
  Duplicates$Alternative_Gene.ID <- NA
  
  #Generate a unique list of Gene.IDs
  TSSIDs <- unique(Duplicates$TSS.ID)
  
  #Condition for while loop and files to hold dominant and nondominant TSS lists
  i = 1
  
  while (i <= length(TSSIDs)){
    
    #Pull All obsTSS annotations with the same gene ID
    obsTSS_perID <- subset(Duplicates, Duplicates$TSS.ID == TSSIDs[i])
    
    #If there are multiple obsTSSs per a given Gene.ID, the dominant obsTSS is the most with the most reads
    keep <- obsTSS_perID[1,]
    if (nrow(obsTSS_perID) == 2 ) {
      keep$Alternative_Gene.ID = obsTSS_perID$Gene.ID[2]
    } else if (nrow(obsTSS_perID) > 2) {
      keep$Alternative_Gene.ID = paste(obsTSS_perID$Gene.ID[2:nrow(obsTSS_perID)], collapse = ";")
    } else {
    }
    
    if (i == 1){
      Updated <- keep 
    } else {
      Updated <- rbind(Updated, keep)
    }
    
    i = i + 1
    
  }
  
  #Export Final Dominant TSS List 
  Nondup <- subset(Dominant_TSS_List, Dominant_TSS_List$Conflict_GeneID_Yes1 == 0 )
  Nondup$Alternative_Gene.ID <- NA
  Final <- rbind(Nondup, Updated)
  
  write.table(Final, paste0(args[2],"_FINAL_Dominant_obsTSS_perGeneID_List.txt"), sep = "\t", row.names = F, col.names = T, quote = F)
  
}

#obsTSS_Annotated
obsTSS_Annotated <- function(prefix){
  
  #Upload Files 
  obsTSS <- read.delim(paste0(prefix,"_TSScall_Output_obsTSS_List.txt"), as.is = T)                     #output from rvmTSS, This lists all active obsTSS (includes both dominant and non dominant TSSs)
  Dominant_TSS <- read.delim(paste0(prefix,"_FINAL_Dominant_obsTSS_perGeneID_List.txt"), as.is = T)     #Output of Dominant_TSScall
  Kallisto <- read.delim(paste0(prefix,"_avg_tpms_all_tx.txt"), as.is = T)      
  #Kallisto <- read.delim(paste0(prefix,"_avg_tpms_kallisto-filt-tx.txt"), as.is = T)      

  #Stable <- subset(Kallisto, Kallisto$AvgTPM > 0 )
  Stable <- subset(Kallisto, Kallisto$AvgTPM > as.numeric(args[5]) )
  
  #Simplify the obsTSS Table 
  obsTSS_Working <- obsTSS[,1:22]
  
  #Cross reference obsTSS IDs found in obsTSS and Dominant TSS to identify the dominant TSS 
  obsTSS_Working$Dominant_obsTSS <- obsTSS_Working$TSS.ID %in% Dominant_TSS$TSS.ID
  
  #Determine if associated with stable transcript 
  obsTSS_Working$Stable <- NA
  
  #For each obsTSS
  for (i in 1:nrow(obsTSS_Working)){

    #deduplicate the transcripts column 
    IDs <- str_split(obsTSS_Working$Transcripts[i], ";", n = Inf, simplify = TRUE) #Break the single string into a matric whos size is dependent on the number of ENSUMGS 
    IDs <- t(IDs) #Transpose the matrix 
    Unique_ID <- unique(IDs)

    #Determine overlap with Kallisto
    if (sum(Unique_ID %in% Stable$aTxID) >= 1){
      obsTSS_Working$Stable[i] = TRUE 
    } else {
      obsTSS_Working$Stable[i] = FALSE 
    }
    
    rm(IDs)
    rm(Unique_ID)
    print(i)
  }
  
  #Final output 
  write.table(obsTSS_Working, paste0(prefix,"_Annotated_Dominant_and_Nondominant_obsTSS_fordREG.txt"), sep = "\t", row.names = F, col.names = T, quote = F)
  
}



##### running functions

#read in TSS call file, remove nuTSSs
TSScall_obs <- rvm_nuTSS(read.delim(args[1], as.is = T))

#get unique ID info
Final <- Unique_ID_Info(TSScall_obs)

#call dominant TSSs
Dominant_TSScall(Final)

#Annotate active obsTSSs as Dominant vs Nondominant, Stable vs Nonstable 
obsTSS_Annotated(args[2])
