library(GEOquery)
library(readr)
setwd("~/work/GIT/RScripts") # has to be removed
setwd("..") #from [GIT]/Rscripts to main git directory
mainDr <- getwd()
source(paste0(mainDr, "/RScripts/outlRemoveGIT.R")) # has to be changed
#############################################################################
geoIDs <- c("GSE120795", "GSE112004")
dataDr = paste0(mainDr, "/Data")
if(!(paste0("./Data") %in% list.dirs(recursive = F))) dir.create(dataDr)
noTypeSpecification <- T  # has to be False when given as input from command line

buLS <- ls()
# geoID <- geoIDs[1]
for(geoID in geoIDs){
  cat("\n\n######################\n", geoID, "\n######################\n")
  
  rm(list = setdiff(ls(), c(buLS, "geoID", "buLS")))
  gc()
  
  #download data
  gse <- getGEO(geoID, GSEMatrix = TRUE)
  smplInfo <- pData(gse[[1]])
  dat <- as.data.frame(exprs(gse[[1]]))
  if(nrow(dat)==0){
    cat("\nDownloading data for ", geoID)
    setwd(dataDr)
    files <- rownames(getGEOSuppFiles(geoID))
    for(fl in files){
      cat("\rImporting data for file ", which(files==fl), " of ", length(files))
      tmp <- read.delim(fl, row.names=1)
      if(fl==files[1]){
        dat <- tmp
      }else{
        overlpngRownames <- intersect(rownames(dat), rownames(tmp))
        dat <- cbind(dat[overlpngRownames,], tmp[overlpngRownames,])
      }
      
    }
    rm(tmp)
    gc()
    #dat <- as.data.frame(t(dat))
  }
  setwd(dataDr)
  cat("\nSaving raw files for dataset ", geoID, " to data directory")
  write.table(dat, paste0(geoID, "DatRAW.txt"), sep="\t")
  write.table(smplInfo, paste0(geoID, "SmplInfoRAW.txt"), sep="\t")
  colnames(dat) <- gsub(".fastq.gz", "", colnames(dat))
   
  #specify type col
  if(any(colnames(smplInfo)=="type")) colnames(smplInfo)[which(colnames(smplInfo)=="type")] <- "oldTypeCol"
  if(noTypeSpecification){
    if(geoID=="GSE120795"){
      titleCol <- "description"
      typeCol <- "source_name_ch1" 
      timeCol <- ""
      colnames(smplInfo)[which(colnames(smplInfo)=="title")] <- "oldTitle"
      colnames(smplInfo)[which(colnames(smplInfo)==titleCol)] <- "title"
    }
    if(geoID=="GSE112004"){
      typeCol <- "characteristics_ch1.3"
      timeCol <- "time point:ch1"
      smplInfo$type <- paste0(smplInfo[,typeCol], smplInfo[,timeCol], sep = ";")
      typeCol <- "type"
    }
  }
  colnames(smplInfo)[which(colnames(smplInfo)==typeCol)] <- "type"
  
  
  # #################################
  # # sample titles
  # #################################
   if(!length(intersect(colnames(dat), smplInfo$title))==nrow(smplInfo)){
    if(length(intersect(colnames(dat), smplInfo$title))<ncol(dat)*0.8){
      if(all(sapply(colnames(dat), function(x) substr(x, 1,1))=="X")){
        smplInfo$title <- paste0("X", smplInfo$title)
      }else{
        if(length(which(sapply(colnames(dat), function(x) substr(x, 1,1))=="X"))>0.5*ncol(dat)){
          smplInfo$title <- sapply(smplInfo$title, function(x) ifelse(!is.na(as.numeric(substr(x, 1,1))), paste0("X", x), x))
        }
      }
      # if(length(intersect(colnames(dat), sapply(as.character(smplInfo$title), function(x) unlist(strsplit(x, split = " "))[1])))>ncol(dat)*0.8){
      #   smplInfo$title <- sapply(as.character(smplInfo$title), function(x) unlist(strsplit(x, split = " "))[1])
      # }
      if(length(intersect(colnames(dat), smplInfo$title))<ncol(dat)*0.8){
        if(length(intersect(colnames(dat), gsub("-",".",smplInfo$title)))>ncol(dat)*0.8){
          smplInfo$title <- gsub("-",".",smplInfo$title)
        }
      }
    }
   }
  
  ################################################################################
  # Match info and data files to same sample order
  smplInfo <- smplInfo[which(smplInfo$title %in% colnames(dat)),]
  dat <- dat[,match(smplInfo$title, colnames(dat))]
  smplInfo <- smplInfo[match(colnames(dat), smplInfo$title),]
  smplInfo$ID <- 1:nrow(smplInfo)
  gc()
  
  geneNames <- rownames(dat)
  
  if(all(colnames(dat)==as.character(smplInfo$title))){
    cat(paste0("\nWriting raw data to file: ",  paste0(geoID,".txt")))
    setwd(dataDr)
    write.table(data.frame(geneNames), paste0(geoID,"GenesRAW.txt"), quote = F, row.names = F,  sep="\t")
    write.table(dat, paste0(geoID,"DatRAW.txt"), quote = F, row.names = F,  sep="\t")
    write.table(smplInfo, paste0(geoID, "SampleInfoRAW.txt"), quote = F, row.names = F,  sep="\t")
  }
}
 