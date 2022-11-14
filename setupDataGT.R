library(GEOquery)
library(readr)
setwd("~/work/GIT/RScripts") # has to be removed
setwd("..") #from [GIT]/Rscripts to main git directory
mainDr <- getwd()
source(paste0(mainDr, "/RScripts/outlRemoveGIT.R")) # has to be changed
#############################################################################
geoIDs <- c("GSE120795", "GSE112004")
dataDr = paste0(mainDr, "/Data")
noTypeSpecification <- T  # has to be False when given as input from command line

buLS <- ls()
# geoID <- geoIDs[1]
for(geoID in geoIDs){
  cat("\n\n######################\n", geoID, "\n######################\n")
  
  rm(list = setdiff(ls(), c(buLS, "geoID", "buLS")))
  gc()
  
  setwd(dataDr)
  if("GSE120795DatRAW.txt" %in% list.files()){
    dat <- read_delim("GSE120795DatRAW.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    dat <- as.data.frame(dat)
    smplInfo <- read_delim("GSE120795SampleInfoRAW.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
    smplInfo <- as.data.frame(smplInfo)
  }else{
    cat("The specified file is not in the Data directory")
  }
  
  if(length(intersect(colnames(dat), smplInfo$title))>0.9*nrow(smplInfo)){
    #backup files
    BUSmplInfo <- smplInfo
    BUdat <- dat
    ################################################################################
    cat("\nOutlier detection...")
    outlOutp <- outlRemove(dat,smplInfo, F, T, T, T, T, T)
    
    #backup files
    smplInfo <- outlOutp[[2]]    
    xs_raw <- outlOutp[[1]]
    
    cat("\rOutlier detection done                    ")
    ################################################################################
    cat("\nSame number of samples per type...")
    sameNrOutp <- sameNrTypes(xs_raw, smplInfo, "type", NA, 8, sameNrSamples = T) 
    
    smplInfo <- sameNrOutp[[2]]    
    xs_raw <- sameNrOutp[[1]]
    
    cat("\rSame number of samples per type done")
    
    ################################################################################
    cat("\nSaving data to ", dataDr, "\n")
    if(all(colnames(xs_raw)==smplInfo$title)){
      prettyTables(dataDf = xs_raw, smplInfoSameSize = smplInfo, outDir = dataDr, tpColNam = "type", tpColNam2 = NA , name = geoID, reOrder = T, F)
    }else{
      stop("rownames xs_raw different than titles")
    }
  }
  cat("\nSaving nice setup files for dataset ", geoID, " to data directory")
  setwd(dataDr)
}
