getphylumTableMetaPhlAn2 <- function(table) {
  table$clade_name[which(table$clade_name %like% "UNCLASSIFIED")] <- sapply(strsplit(table$clade_name[which(table$clade_name %like% "UNCLASSIFIED")], "\\|c__"), "[", 1)
   taxaTable <- table[(table$clade_name %like% "p__"),]
  taxaTable <- taxaTable[!(taxaTable$clade_name %like% "c__"),]
  taxaTable <- taxaTable[!(taxaTable$clade_name %like% "o__"),]
  taxaTable <- taxaTable[!(taxaTable$clade_name %like% "f__"),]
  taxaTable <- taxaTable[!(taxaTable$clade_name %like% "g__"),]
  taxaTable <- taxaTable[!(taxaTable$clade_name %like% "s__"),]
  names(taxaTable)[names(taxaTable) == "clade_name"] <- "phylum"
  return(taxaTable)
}

getclassTableMetaPhlAn2 <- function(table) {
   table$clade_name[which(table$clade_name %like% "UNCLASSIFIED")] <- sapply(strsplit(table$clade_name[which(table$clade_name %like% "UNCLASSIFIED")], "\\|o__"), "[", 1)
   taxaTable <- table[(table$clade_name %like% "p__"),]
  taxaTable <- taxaTable[(taxaTable$clade_name %like% "c__"),]
  taxaTable <- taxaTable[!(taxaTable$clade_name %like% "o__"),]
  taxaTable <- taxaTable[!(taxaTable$clade_name %like% "f__"),]
  taxaTable <- taxaTable[!(taxaTable$clade_name %like% "g__"),]
  taxaTable <- taxaTable[!(taxaTable$clade_name %like% "s__"),]
  names(taxaTable)[names(taxaTable) == "clade_name"] <- "class"
  return(taxaTable)
}

getorderTableMetaPhlAn2 <- function(table) {
   table$clade_name[which(table$clade_name %like% "UNCLASSIFIED")] <- sapply(strsplit(table$clade_name[which(table$clade_name %like% "UNCLASSIFIED")], "\\|f__"), "[", 1)
   
  taxaTable <- table[(table$clade_name %like% "p__"),]
  taxaTable <- taxaTable[(taxaTable$clade_name %like% "c__"),]
  taxaTable <- taxaTable[(taxaTable$clade_name %like% "o__"),]
  taxaTable <- taxaTable[!(taxaTable$clade_name %like% "f__"),]
  taxaTable <- taxaTable[!(taxaTable$clade_name %like% "g__"),]
  taxaTable <- taxaTable[!(taxaTable$clade_name %like% "s__"),]
  names(taxaTable)[names(taxaTable) == "clade_name"] <- "order"
  return(taxaTable)
}

getfamilyTableMetaPhlAn2 <- function(table) {
   table$clade_name[which(table$clade_name %like% "UNCLASSIFIED")] <- sapply(strsplit(table$clade_name[which(table$clade_name %like% "UNCLASSIFIED")], "\\|g__"), "[", 1)
   
  taxaTable <- table[(table$clade_name %like% "p__"),]
  taxaTable <- taxaTable[(taxaTable$clade_name %like% "c__"),]
  taxaTable <- taxaTable[(taxaTable$clade_name %like% "o__"),]
  taxaTable <- taxaTable[(taxaTable$clade_name %like% "f__"),]
  taxaTable <- taxaTable[!(taxaTable$clade_name %like% "g__"),]
  taxaTable <- taxaTable[!(taxaTable$clade_name %like% "s__"),]
  names(taxaTable)[names(taxaTable) == "clade_name"] <- "family"
  return(taxaTable)
}

getgenusTableMetaPhlAn2 <- function(table) {
   table$clade_name[which(table$clade_name %like% "UNCLASSIFIED")] <- sapply(strsplit(table$clade_name[which(table$clade_name %like% "UNCLASSIFIED")], "\\|s__"), "[", 1)
   
  taxaTable <- table[(table$clade_name %like% "p__"),]
  taxaTable <- taxaTable[(taxaTable$clade_name %like% "c__"),]
  taxaTable <- taxaTable[(taxaTable$clade_name %like% "o__"),]
  taxaTable <- taxaTable[(taxaTable$clade_name %like% "f__"),]
  taxaTable <- taxaTable[(taxaTable$clade_name %like% "g__"),]
  taxaTable <- taxaTable[!(taxaTable$clade_name %like% "s__"),]
  names(taxaTable)[names(taxaTable) == "clade_name"] <- "genus"
  return(taxaTable)
}

getspeciesTableMetaPhlAn2 <- function(table) {
   table$clade_name[which(table$clade_name %like% "UNCLASSIFIED")] <- sapply(strsplit(table$clade_name[which(table$clade_name %like% "UNCLASSIFIED")], "\\|t__"), "[", 1)
   
  taxaTable <- table[(table$clade_name %like% "p__"),]
  taxaTable <- taxaTable[(taxaTable$clade_name %like% "c__"),]
  taxaTable <- taxaTable[(taxaTable$clade_name %like% "o__"),]
  taxaTable <- taxaTable[(taxaTable$clade_name %like% "f__"),]
  taxaTable <- taxaTable[(taxaTable$clade_name %like% "g__"),]
  taxaTable <- taxaTable[(taxaTable$clade_name %like% "s__"),]
  taxaTable <- taxaTable[!(taxaTable$clade_name %like% "t__"),]
  names(taxaTable)[names(taxaTable) == "clade_name"] <- "species"
  return(taxaTable)
}

getstrainTableMetaPhlAn2 <- function(table) {
   # table$clade_name[which(table$clade_name %like% "UNCLASSIFIED")] <- sapply(strsplit(table$clade_name[which(table$clade_name %like% "UNCLASSIFIED")], "\\|t__"), "[", 1)
   
   taxaTable <- table[(table$clade_name %like% "p__"),]
   taxaTable <- taxaTable[(taxaTable$clade_name %like% "c__"),]
   taxaTable <- taxaTable[(taxaTable$clade_name %like% "o__"),]
   taxaTable <- taxaTable[(taxaTable$clade_name %like% "f__"),]
   taxaTable <- taxaTable[(taxaTable$clade_name %like% "g__"),]
   taxaTable <- taxaTable[(taxaTable$clade_name %like% "s__"),]
   taxaTable <- taxaTable[(taxaTable$clade_name %like% "t__"),]
   names(taxaTable)[names(taxaTable) == "clade_name"] <- "strain"
   return(taxaTable)
}

transposeTaxaTable <- function(table) {
  rownames(table) <- table[,1]
  table <- table[,-1]
  table <- t(table)
  return(table)
}

#Normalizes taxonomic tables
getNormTable <- function(table) {
  
  Normalized=table
  
  for (x in 1:nrow(table)){
    Normalized[x,1:ncol(table)] = 
      log10(((table[x,1:ncol(table)]) / (rowSums(table[,1:ncol(table)])[x]) * rowMeans(table[,1:ncol(table)])[x]) + 1)
  }
  
  return(Normalized)
  
}

getRelAbun <- function(table) {
  
  relabun <- table
  
  for (x in 1:nrow(table)){
    relabun[x,1:ncol(table)] = ((table[x,1:ncol(table)]) / (rowSums(table[,1:ncol(table)])[x]))
  }
  
  return(relabun)
  
}
##### This function rounds p values based on value #####
roundP <- function(pvalue) {
  
  if (is.nan(pvalue)) {
    message("TRUE")
  } else if (pvalue==0){
    p = "p < 2.2e-16"
  } else if (pvalue < 0.001) {
    p <- "p < 0.001"
  } else {
    p <- paste0("p = ", round(pvalue, digits = 3))
  }
  
}


##### This function assigns significance *s based on p value #####
sigStars <- function(pvalue) {
  
  if (pvalue < 0.05) {
    pStar <- "*"
  } else if (pvalue < 0.01){
    pStar <- "**"
  } else if (pvalue < 0.001){
    pStar <- "***"
  } else {
    pStar <- ""
  }
  
}


##### This function corrects the sign of p-values based on coefficient sign #####
getlog10p<-function(pval,coefficient){
  
  log10p<-sapply(1:length(pval), function(x){
    
    if (coefficient[x]>0) return(-log10(pval[x]))
    else return(log10(pval[x]))
  })
}


##### This function parses out concatenated taxa names. #####
  parseMetaPhlAn2Taxa<-function(table){
    
    lev <- colnames(table)[1]
    # table <- table[!(table[,1] %like% "k__Viruses"),]
    # table <- table[!(table[,1] %like% "k__Eukaryota"),]
    
    tempString <- table[,1]
    tempString <- gsub("[|]", "_", tempString)
    tempString <- str_split(tempString, pattern = "_[pcofgs]__")
    tempString=do.call(rbind,tempString)
    tempString <- as.data.frame(tempString)
    
    colnames(tempString)[colnames(tempString)=="V1"] <- "Domain"
    colnames(tempString)[colnames(tempString)=="V2"] <- "Phylum"
    colnames(tempString)[colnames(tempString)=="V3"] <- "Class"
    colnames(tempString)[colnames(tempString)=="V4"] <- "Order"
    colnames(tempString)[colnames(tempString)=="V5"] <- "Family"
    colnames(tempString)[colnames(tempString)=="V6"] <- "Genus"
    colnames(tempString)[colnames(tempString)=="V7"] <- "Species"
    
    table <- cbind(tempString, table)
    L <- which(colnames(table) == lev)-1
    rownames(table) <- table[,L]
    L1 <- L+1
    table <- table[,-(1:L1)]
    table <- t(table)
    table[is.na(table)] <- 0
    SampleID <- rownames(table)
    table <- cbind(SampleID, table)
    table <- as.data.frame(table)
    return(table)
    
  }


  ##### This function parses out concatenated taxa names. #####
  parseMetaPhlAn2Taxa<-function(tableInput){
    
    names(tableInput)[1] <- "V1"
    tempString <- tableInput[,1]
    tempString <- gsub("[|]", "_", tempString)
    tempString <- str_split(tempString, pattern = "_[pcofgs]__")
    tempString=do.call(rbind,tempString)
    tempString <- as.data.frame(tempString)
    
    colnames(tempString)[colnames(tempString)=="V1"] <- "Domain"
    colnames(tempString)[colnames(tempString)=="V2"] <- "Phylum"
    colnames(tempString)[colnames(tempString)=="V3"] <- "Class"
    colnames(tempString)[colnames(tempString)=="V4"] <- "Order"
    colnames(tempString)[colnames(tempString)=="V5"] <- "Family"
    colnames(tempString)[colnames(tempString)=="V6"] <- "Genus"
    colnames(tempString)[colnames(tempString)=="V7"] <- "Species"
    
    L <- ncol(tempString)
    lev <- colnames(tempString)[L]
    
    if (L == 7) {
      tempString$Strain <- sapply(strsplit(tempString$Species, "_t__"), "[", 2)
      tempString$Species <- sapply(strsplit(tempString$Species, "_t__"), "[", 1)
    }
    
    table <- cbind(tempString, tableInput)
    
    Unclassified <- which(table[,L] %like% "_unclassified")
    table[Unclassified, L] <- "Unclassified"
    
    NoName <- which(table[,L] %like% "_noname")
    table[NoName, L] <- "Unclassified"
    
    SampleStart <- which(colnames(table) == "V1") + 1
    
    # Loop through each column and convert the class
    for (col in SampleStart:ncol(table)) {
      table[, col] <- as.numeric(table[, col])
      table[is.na(table[, col]), col] <- 0
    }    

    agg_df <- aggregate(table[,SampleStart:ncol(table)], by=list(table[,L]), FUN = sum)
    rownames(agg_df) <- agg_df[,1]
    agg_df <- agg_df[,-1]
    # L1 <- L+1
    # table <- table[,-(1:L1)]
    table <- t(agg_df)
    # table[is.na(table)] <- 0
    SampleID <- rownames(table)
    table <- cbind(SampleID, table)
    table <- as.data.frame(table)
    
    # Loop through each column and convert the class
    for (col in 2:ncol(table)) {
      table[, col] <- as.numeric(table[, col])
    }    
    
    return(table)
    
  }

















getphylumTableKraken2 <- function(table) {
      table$clade_name[which(table$clade_name %like% "UNCLASSIFIED")] <- sapply(strsplit(table$clade_name[which(table$clade_name %like% "UNCLASSIFIED")], "\\|t__"), "[", 1)

   taxaTable <- table[(table$V1 %like% "p__"),]
   taxaTable <- taxaTable[!(taxaTable$V1 %like% "c__"),]
   taxaTable <- taxaTable[!(taxaTable$V1 %like% "o__"),]
   taxaTable <- taxaTable[!(taxaTable$V1 %like% "f__"),]
   taxaTable <- taxaTable[!(taxaTable$V1 %like% "g__"),]
   taxaTable <- taxaTable[!(taxaTable$V1 %like% "s__"),]
   names(taxaTable)[names(taxaTable) == "V2"] <- sampleID
   return(taxaTable)
}

getclassTableKraken2 <- function(table) {
   taxaTable <- table[(table$V1 %like% "p__"),]
   taxaTable <- taxaTable[(taxaTable$V1 %like% "c__"),]
   taxaTable <- taxaTable[!(taxaTable$V1 %like% "o__"),]
   taxaTable <- taxaTable[!(taxaTable$V1 %like% "f__"),]
   taxaTable <- taxaTable[!(taxaTable$V1 %like% "g__"),]
   taxaTable <- taxaTable[!(taxaTable$V1 %like% "s__"),]
   names(taxaTable)[names(taxaTable) == "V2"] <- sampleID
   return(taxaTable)
}

getorderTableKraken2 <- function(table) {
   taxaTable <- table[(table$V1 %like% "p__"),]
   taxaTable <- taxaTable[(taxaTable$V1 %like% "c__"),]
   taxaTable <- taxaTable[(taxaTable$V1 %like% "o__"),]
   taxaTable <- taxaTable[!(taxaTable$V1 %like% "f__"),]
   taxaTable <- taxaTable[!(taxaTable$V1 %like% "g__"),]
   taxaTable <- taxaTable[!(taxaTable$V1 %like% "s__"),]
   names(taxaTable)[names(taxaTable) == "V2"] <- sampleID
   return(taxaTable)
}

getfamilyTableKraken2 <- function(table) {
   taxaTable <- table[(table$V1 %like% "p__"),]
   taxaTable <- taxaTable[(taxaTable$V1 %like% "c__"),]
   taxaTable <- taxaTable[(taxaTable$V1 %like% "o__"),]
   taxaTable <- taxaTable[(taxaTable$V1 %like% "f__"),]
   taxaTable <- taxaTable[!(taxaTable$V1 %like% "g__"),]
   taxaTable <- taxaTable[!(taxaTable$V1 %like% "s__"),]
   names(taxaTable)[names(taxaTable) == "V2"] <- sampleID
   return(taxaTable)
}

getgenusTableKraken2 <- function(table) {
   taxaTable <- table[(table$V1 %like% "p__"),]
   taxaTable <- taxaTable[(taxaTable$V1 %like% "c__"),]
   taxaTable <- taxaTable[(taxaTable$V1 %like% "o__"),]
   taxaTable <- taxaTable[(taxaTable$V1 %like% "f__"),]
   taxaTable <- taxaTable[(taxaTable$V1 %like% "g__"),]
   taxaTable <- taxaTable[!(taxaTable$V1 %like% "s__"),]
   names(taxaTable)[names(taxaTable) == "V2"] <- sampleID
   return(taxaTable)
}

getspeciesTableKraken2 <- function(table) {
   taxaTable <- table[(table$V1 %like% "p__"),]
   taxaTable <- taxaTable[(taxaTable$V1 %like% "c__"),]
   taxaTable <- taxaTable[(taxaTable$V1 %like% "o__"),]
   taxaTable <- taxaTable[(taxaTable$V1 %like% "f__"),]
   taxaTable <- taxaTable[(taxaTable$V1 %like% "g__"),]
   taxaTable <- taxaTable[(taxaTable$V1 %like% "s__"),]
   taxaTable <- taxaTable[!(taxaTable$V1 %like% "t__"),]
   names(taxaTable)[names(taxaTable) == "V2"] <- sampleID
   return(taxaTable)
}


parseKraken2Taxa<-function(table){
   
   
   table <- table[!(table$V1 %like% "d__Viruses"),]
   table <- table[!(table$V1 %like% "d__Eukaryota"),]
   
   tempString <- table$V1
   tempString <- gsub("[|]", "_", tempString)
   tempString <- str_split(tempString, pattern = "_[pcofgs]__")
   tempString=do.call(rbind,tempString)
   tempString <- as.data.frame(tempString)
   
   colnames(tempString)[colnames(tempString)=="V1"] <- "Domain"
   colnames(tempString)[colnames(tempString)=="V2"] <- "Phylum"
   colnames(tempString)[colnames(tempString)=="V3"] <- "Class"
   colnames(tempString)[colnames(tempString)=="V4"] <- "Order"
   colnames(tempString)[colnames(tempString)=="V5"] <- "Family"
   colnames(tempString)[colnames(tempString)=="V6"] <- "Genus"
   colnames(tempString)[colnames(tempString)=="V7"] <- "Species"
   
   table <- cbind(tempString, table)
   L <- which(colnames(table) == "V1")-1
   rownames(table) <- table[,L]
   L1 <- L+1
   table <- table[,-(1:L1)]
   table <- t(table)
   table[is.na(table)] <- 0
   SampleID <- rownames(table)
   table <- cbind(SampleID, table)
   return(table)
   
}

parseMetaPhlAn2Taxa<-function(table){
   
   lev <- colnames(table)[1]
   # table <- table[!(table[,1] %like% "k__Viruses"),]
   # table <- table[!(table[,1] %like% "k__Eukaryota"),]
   
   tempString <- table[,1]
   tempString <- gsub("[|]", "_", tempString)
   tempString <- str_split(tempString, pattern = "_[pcofgs]__")
   tempString=do.call(rbind,tempString)
   tempString <- as.data.frame(tempString)
   
   colnames(tempString)[colnames(tempString)=="V1"] <- "Domain"
   colnames(tempString)[colnames(tempString)=="V2"] <- "Phylum"
   colnames(tempString)[colnames(tempString)=="V3"] <- "Class"
   colnames(tempString)[colnames(tempString)=="V4"] <- "Order"
   colnames(tempString)[colnames(tempString)=="V5"] <- "Family"
   colnames(tempString)[colnames(tempString)=="V6"] <- "Genus"
   colnames(tempString)[colnames(tempString)=="V7"] <- "Species"
   
   table <- cbind(tempString, table)
   L <- which(colnames(table) == lev)-1
   rownames(table) <- table[,L]
   L1 <- L+1
   table <- table[,-(1:L1)]
   table <- t(table)
   table[is.na(table)] <- 0
   SampleID <- rownames(table)
   table <- cbind(SampleID, table)
   table <- as.data.frame(table)
   return(table)
   
}

parseMetaPhlAn2Taxa<-function(tableInput){
   
   names(tableInput)[1] <- "V1"
   tempString <- tableInput[,1]
   tempString <- gsub("[|]", "_", tempString)
   tempString <- str_split(tempString, pattern = "_[pcofgs]__")
   tempString=do.call(rbind,tempString)
   tempString <- as.data.frame(tempString)
   
   colnames(tempString)[colnames(tempString)=="V1"] <- "Domain"
   colnames(tempString)[colnames(tempString)=="V2"] <- "Phylum"
   colnames(tempString)[colnames(tempString)=="V3"] <- "Class"
   colnames(tempString)[colnames(tempString)=="V4"] <- "Order"
   colnames(tempString)[colnames(tempString)=="V5"] <- "Family"
   colnames(tempString)[colnames(tempString)=="V6"] <- "Genus"
   colnames(tempString)[colnames(tempString)=="V7"] <- "Species"
   
   L <- ncol(tempString)
   lev <- colnames(tempString)[L]
   
   if (L == 7) {
      tempString$Strain <- sapply(strsplit(tempString$Species, "_t__"), "[", 2)
      tempString$Species <- sapply(strsplit(tempString$Species, "_t__"), "[", 1)
   }
   
   table <- cbind(tempString, tableInput)
   
   Unclassified <- which(table[,L] %like% "_unclassified")
   table[Unclassified, L] <- "Unclassified"
   
   NoName <- which(table[,L] %like% "_noname")
   table[NoName, L] <- "Unclassified"
   
   SampleStart <- which(colnames(table) == "V1") + 1
   
   # Loop through each column and convert the class
   for (col in SampleStart:ncol(table)) {
      table[, col] <- as.numeric(table[, col])
      table[is.na(table[, col]), col] <- 0
   }    
   
   agg_df <- aggregate(table[,SampleStart:ncol(table)], by=list(table[,L]), FUN = sum)
   rownames(agg_df) <- agg_df[,1]
   agg_df <- agg_df[,-1]
   # L1 <- L+1
   # table <- table[,-(1:L1)]
   table <- t(agg_df)
   # table[is.na(table)] <- 0
   SampleID <- rownames(table)
   table <- cbind(SampleID, table)
   table <- as.data.frame(table)
   
   # Loop through each column and convert the class
   for (col in 2:ncol(table)) {
      table[, col] <- as.numeric(table[, col])
   }    
   
   return(table)
   
}
