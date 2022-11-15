outlRemove <- function(dataDf, smplInfo, plotChecks, sumsOutGPercentile, ratioOutG, sumsOutS, sumsOutSPercentile, ratiooutS, naOm=T){
  # plotChecks = F
  # sumsOutGPercentile = ratioOutG = sumsOutS = sumsOutSPercentile = ratiooutS = naOm=T
  cat("\nOutlier detection...")
  #sumsOutGPercentile <-  ratioOutG <-  sumsOutS <- sumsOutSPercentile <- ratiooutS <- T
  #plotChecks <- F
  if(ncol(dataDf)>nrow(dataDf)){
    buSampNames <- rownames(dataDf)
    buGeneNames <- colnames(dataDf)
    dataDf <- data.frame(t(dataDf))
    rownames(dataDf) <- buGeneNames
    colnames(dataDf) <- buSampNames
  }
  if(ncol(dataDf)>500){
    steps <- round(seq(0, ceiling(ncol(dataDf)), length.out =  ceiling(ncol(dataDf)/500)))
    steps2 <- round(seq(0, ceiling(nrow(dataDf)), length.out =  length(steps)))
  }else{
    steps <- c(0,ncol(dataDf))
    steps2 <- c(0,nrow(dataDf))
  }  
  ########################################################
  # Genes
  ########################################################
  # if(ratioOutG || sumsOutGPercentile){
  allRatiosG <- c()
  allSumsG <- c()
  # allSumsS <- c()
  cat("\nChecking genes\n")
  for(i in 2:length(steps2)){
    cat("\r", i-1, " of ", length(steps)-1)
    if(ratioOutG) ratiosGenes <- apply(dataDf[(steps2[i-1]+1):steps2[i],], 1, function(x) length(which(x==0))/length(x))
    sumsGenes <- apply(dataDf[(steps2[i-1]+1):steps2[i],], 1, sum)
    # }else{
    #   ratiosSampels <- apply(dataDf[,(steps[i-1]+1):steps[i]], 2, function(x) length(which(x<=0))/length(x))
    #   sumsSamples <- apply(dataDf[,(steps[i-1]+1):steps[i]], 2, sum)
    #   ratiosGenes <- apply(dataDf[(steps2[i-1]+1):steps2[i],], 2, function(x) length(which(x<=0))/length(x))
    #   sumsGenes <- apply(dataDf[(steps2[i-1]+1):steps2[i],], 1, sum)
    # }
    gc()
    if(ratioOutG) allRatiosG <- c(allRatiosG, ratiosGenes)
    allSumsG <- c(allSumsG, sumsGenes)
    # allSumsS <- c(allSumsS, sumsSamples)
  }
  ###################
  if(plotChecks){
    hist(allRatiosG, 1000)
    abline(v=c(median(allRatiosG, na.rm = T)+5*sd(allRatiosG, na.rm = T)), col="red")
    abline(v=c(median(allRatiosG, na.rm = T)+2*sd(allRatiosG, na.rm = T)), col="blue")
    abline(v=c(median(allRatiosG, na.rm = T)+sd(allRatiosG, na.rm = T)), col="orange")
    #b
    hist(allSumsG, 1000)
    abline(v=quantile(allSumsG, 0.995, na.rm = T), col="darkgreen")
    hist(allSumsG[intersect(which(allSumsG<quantile(allSumsG, 0.995, na.rm = T)), which(allSumsG!=0))], 1000)
    # abline(v=c((median(allSumsG[intersect(which(allSumsG<quantile(allSumsG, 0.995)), which(allSumsG!=0))])-5*sd(allSumsG[intersect(which(allSumsG<quantile(allSumsG, 0.995)), which(allSumsG!=0))])), median(allSumsG[intersect(which(allSumsG<quantile(allSumsG, 0.995)), which(allSumsG!=0))])+5*sd(allSumsG[intersect(which(allSumsG<quantile(allSumsG, 0.995)), which(allSumsG!=0))])), col="red")
    # abline(v=c((median(allSumsG[intersect(which(allSumsG<quantile(allSumsG, 0.995)), which(allSumsG!=0))])-2*sd(allSumsG[intersect(which(allSumsG<quantile(allSumsG, 0.995)), which(allSumsG!=0))])), median(allSumsG[intersect(which(allSumsG<quantile(allSumsG, 0.995)), which(allSumsG!=0))])+2*sd(allSumsG[intersect(which(allSumsG<quantile(allSumsG, 0.995)), which(allSumsG!=0))])), col="blue")
  }
  ###################
  outGenesZeroSum <- names(allSumsG[allSumsG==0])
  cat("\n", length(outGenesZeroSum), " genes out due to sum==0")
  if(sumsOutGPercentile){
    outGenesLrgSum <- c(names(allSumsG[allSumsG>quantile(allSumsG, 0.995, na.rm=T)]))
  }
  cat("\n", length(outGenesLrgSum), " genes out due to large sum")
  
  
  if(ratioOutG){
    #name <- paste0(name, "OUTLG")
    # if(ncol(dataDf) < 200){
    #   outGenes <- c(outGenes, names(allRatiosG[which(allRatiosG > (1-0.5*largeNr/(largeNr*length(unique(smplInfo[,tpCol])))))]))  # ?(median(allRatiosG, na.rm=T)+sd(allRatiosG, na.rm=T))
    # }else{
    outGenesLrgZeroRatio <- setdiff(names(allRatiosG[which(allRatiosG > 0.995)]), outGenesZeroSum)
    # outGenes <- c(outGenes, addOutGenes)  # ?(median(allRatiosG, na.rm=T)+sd(allRatiosG, na.rm=T))
    # }
  }else{
    outGenesLrgZeroRatio <- c()
  }  
  cat("\n", length(outGenesLrgZeroRatio), " genes out due to high zero-ratio")

  
  outGenes <- c(outGenesZeroSum, outGenesLrgSum, outGenesLrgZeroRatio)
  cat("\n", length(unique(outGenes)), " unique total genes out")
  
  # }
  
  # currentDir <- getwd()
  # setwd("//mi-buffer/ag_medbioinf/mrams/work/2017/IBOSS/IBOSS-data")
  # dir.create(orName)
  # setwd(orName)
  # save(list = c("outGenes", "allSumsG"), file = paste0(name, "_outGenes.RData"))
  # setwd(currentDir)
  
  cat("\nPercentage out genes:", 100*round(length(unique(outGenes))/nrow(dataDf), 3))
  gc()
  if(any(rownames(dataDf) %in% unique(outGenes))) dataDf <- dataDf[-which(rownames(dataDf) %in% unique(outGenes)),]
  gc()
  
  ########################################################
  # Samples
  ########################################################
  # if(ratiooutS || sumsOutS || sumsOutSPercentile){
    #name <- paste0(name, "OUTLS")
    if(ratiooutS) allRatiosS <- c()
    allSumsS <- c()
    allMediansS <- c()
    cat("\nChecking sampels\n")
    for(i in 2:length(steps)){
      cat("\r", i-1, " of ", length(steps)-1)
      # if(i<length(steps)){
      #check samples
      if(ratiooutS) ratiosSamples <- apply(dataDf[,(steps[i-1]+1):steps[i]], 2, function(x) length(which(x<=0))/length(x))
      sumsSamples <- apply(dataDf[,(steps[i-1]+1):steps[i]], 2, sum)
      mediansSamples <- apply(dataDf[,(steps[i-1]+1):steps[i]], 2, median)
      if(ratiooutS) allRatiosS <- c(allRatiosS, ratiosSamples)
      allSumsS <- c(allSumsS, sumsSamples)
      allMediansS <- c(allMediansS, mediansSamples)
    }
    
    ##################################
    # check 
    ##################################
    plotChecks <- F
    if(plotChecks){
      # set.seed(2)
      # shapiro.test(sample(sort(allRatiosS), size = 5000))$p.value
      # set.seed(2)
      # shapiro.test(sample(sort(allSumsS), size = 5000))$p.valueqqnorm(words1)
      # set.seed(2)
      # shapiro.test(sample(sort(allRatiosG), size = 5000))$p.value
      # set.seed(2)
      # shapiro.test(sample(sort(allSumsG), size = 5000))$p.value
      hist(allRatiosS, 1000)
      abline(v=c(median(allRatiosS)+5*sd(allRatiosS), median(allRatiosS)-5*sd(allRatiosS)), col="red")
      abline(v=c(median(allRatiosS)+3*sd(allRatiosS), median(allRatiosS)-3*sd(allRatiosS)), col="blue")
      abline(v=c(median(allRatiosS)+4*sd(allRatiosS), median(allRatiosS)-4*sd(allRatiosS)), col="grey")
      abline(v=c(median(allRatiosS)+sd(allRatiosS), median(allRatiosS)-sd(allRatiosS)), col="pink")
      
      #a.double check
      smallOut <- names(allRatiosS[which(allRatiosS<median(allRatiosS)-3*sd(allRatiosS))])
      sO <- smplInfo$organism_part[which(smplInfo[,1] %in% smallOut)]
      table(sO)
      largeOut <- names(allRatiosS[which(allRatiosS>median(allRatiosS)+3*sd(allRatiosS))])
      lO <- smplInfo$organism_part[which(smplInfo[,1] %in% largeOut)]
      table(lO)
      # #b
      hist(allSumsS, 1000)
      abline(v=c((median(allSumsS, na.rm = T)-5*sd(allSumsS, na.rm = T)), median(allSumsS, na.rm = T)+5*sd(allSumsS, na.rm = T)), col="red")
      abline(v=c((median(allSumsS, na.rm = T)-3*sd(allSumsS, na.rm = T)), median(allSumsS, na.rm = T)+3*sd(allSumsS, na.rm = T)), col="blue")
      abline(v=c((median(allSumsS, na.rm = T)-2*sd(allSumsS, na.rm = T)), median(allSumsS, na.rm = T)+2*sd(allSumsS, na.rm = T)), col="grey")
      abline(v=c((median(allSumsS, na.rm = T)-sd(allSumsS, na.rm = T)), median(allSumsS, na.rm = T)+sd(allSumsS, na.rm = T)), col="pink")
      #b.double check
      smallOut <- names(allSumsS[which(allSumsS<median(allSumsS, na.rm = T)-3*sd(allSumsS, na.rm = T))])
      sO <- smplInfo$organism_part[which(smplInfo[,1] %in% smallOut)]
      table(sO)
      largeOut <- names(allSumsS[which(allSumsS>median(allSumsS, na.rm = T)+3*sd(allSumsS, na.rm = T))])
      lO <- smplInfo$organism_part[which(smplInfo[,1] %in% largeOut)]
      table(lO)
      
      plot(allMediansS)
      abline(1.5*median(allMediansS),0)
      abline(0.5*median(allMediansS),0)
    }
    ##################################
    # remove 
    ##################################
    # if(ncol(dataDf) < 200){
    #   outSamples <- names(allRatiosS[which(allRatiosS < (median(allRatiosS, na.rm=T)-2*sd(allRatiosS, na.rm=T)))])
    #   outSamples <- c(outSamples, names(allRatiosS[which(allRatiosS > (median(allRatiosS, na.rm=T)+2*sd(allRatiosS, na.rm=T)))]))
    #   outSamples <- c(outSamples, names(allSumsS[which(allSumsS < (median(allSumsS, na.rm=T)-2*sd(allSumsS, na.rm=T)))]))
    #   outSamples <- c(outSamples, names(allSumsS[which(allSumsS > (median(allSumsS, na.rm=T)+2*sd(allSumsS, na.rm=T)))]))
    # }else{
    
    outSamplesZeroSum <- names(allSumsS[allSumsS==0])   
    outSamplesLrgSum <- names(allSumsS[which(allSumsS>=quantile(allSumsS, 0.995, na.rm=T))])   
    if(ratiooutS){
      #outSamples <- c(outSamples, names(allRatiosS[which(allRatiosS < (median(allRatiosS, na.rm=T)-3*sd(allRatiosS, na.rm=T)))]))
      outSamplesLrgZeroRatio <- setdiff(names(allRatiosS[which(allRatiosS > 0.995)]), outSamplesZeroSum) #(median(allRatiosS, na.rm=T)+3*sd(allRatiosS, na.rm=T)))]))
      #name <- paste0(name, "rat")
    }else{
      outSamplesLrgZeroRatio <- c()
    }
    
    if(sumsOutS){
      #name <- paste0(name, "sums")
      #outSamples <- c(outSamples, names(allSumsS[which(allSumsS < (median(allSumsS, na.rm=T)-3*sd(allSumsS, na.rm=T)))]))
      # if(sumsOutSPercentile){
        outSamplesLrgSum <- names(allSumsS[which(allSumsS > quantile(allSumsS, 0.995))])
       # name <- paste0(name, "Per")
      # }else{
      #   outSamples <- c(outSamples, names(allSumsS[which(allSumsS > (median(allSumsS, na.rm=T)+3*sd(allSumsS, na.rm=T)))])) 
      # }
    }
    
    #outSamples <- c(outSamples, names(allMediansS[c(which(allMediansS > 1.5*median(allMediansS)), which(allMediansS <0.5*median(allMediansS)))]))
    
    cat("\n", length(outSamplesZeroSum), " samples out due to sum==0")
    cat("\n", length(outSamplesLrgSum), " samples out due to large sum")
    cat("\n", length(outSamplesLrgZeroRatio), " samples out due to high zero-ratio")
    
    outSamples <- c(outSamplesZeroSum, outSamplesLrgSum, outSamplesLrgZeroRatio)
    
    cat("\nTotal unique samples:", length(unique(outSamples)))
    cat("\nPercentage out samples:", 100*round(length(unique(outSamples))/ncol(dataDf), 3))
    gc()
    
    ##################
    # setwd("//mi-buffer/ag_medbioinf/mrams/work/2017/IBOSS/IBOSS-data")
    # dir.create(name)
    # setwd(name)
    # save.image(paste0(name, "TablesNoOutl.RData"))
    ##################
    #load("//mi-buffer/ag_medbioinf/mrams/work/2017/IBOSS/IBOSS-data/ExprAOUTL/ExprAtlasTablesNoOutl.RData")
    if(any(colnames(dataDf) %in% unique(outSamples))) dataDf <- dataDf[,-which(colnames(dataDf) %in% unique(outSamples))]
#    smplInfo <- smplInfo[which(smplInfo[,1] %in% colnames(dataDf)),]
    gc()
    
    
    # for(i in 2:length(steps)){
    #   cat("\r", i-1, " of ", length(steps)-1)
    #   if(i<nrow(dataDf)){
    #     ratios <- apply(dataDf[,(steps[i-1]+1):steps[i]], 2, function(x) length(which(x<=0))/length(x))
    #   }else{
    #     ratios <- apply(dataDf[,(steps[i-1]+1):steps[i]], 2, function(x) length(which(x<=0))/length(x))
    #   }
    #   outSamples <- c(outSamples, names(which(ratios < (median(medians)-5*median(sds)))))
    # }
    
    ## }else{
    ##   sumsGenes <- apply(dataDf, 1, function(x) length(which(x<=0))/length(x))
    ##   sumsSamples <- apply(dataDf, 2, function(x) length(which(x<=0))/length(x))
    ##   ratios <- apply(dataDf, 2, function(x) length(which(x<=0))/length(x))
    ##   outSamples <- names(which(ratios < (median(ratios)-5*sd(ratios))))
    ## }
  
  cat("\ndim before na.omit: ", dim(dataDf))
  if(naOm){
    dataDf <- na.omit(dataDf)
    # dataDfNAo1 <- na.omit(dataDf)
    # dataDfNAo2 <- data.frame(t(na.omit(t(dataDf))))
    # dataDf <- ifelse(length(as.matrix(dataDfNAo1))>length(as.matrix(dataDfNAo2)), dataDfNAo1, dataDfNAo2)
  }
  smplInfo <- smplInfo[which(smplInfo$title %in% colnames(dataDf)),]
  cat("\ndim after na.omit: ", dim(dataDf), "\n")
  return(list(dataDf, smplInfo, allRatiosG))
}

sameNrTypes <- function(dataDf, smplInfo, tpColNam, tpColNam2, predefinedMin=NA, sameNrSamples){
  if(ncol(dataDf)>nrow(dataDf)){
    buSampNames <- rownames(dataDf)
    buGeneNames <- colnames(dataDf)
    dataDf <- data.frame(t(dataDf))
    rownames(dataDf) <- buGeneNames
    colnames(dataDf) <- buSampNames
  }
  # cat('\n', predefinedMin)
  cat('\n', dim(dataDf))
  #change organis_part2 <->  organis_part
  
  # if(detailed){
  #   tpColNam <-  colnames(smplInfo)[which(grepl("etail", colnames(smplInfo)))]
  #   saveName <- paste0(name, "detail")
  # }else{
  # saveName <- name
  # }
  ###################################################################################
  # same number samples
  ###################################################################################
  if(sameNrSamples){
    if(!is.na(tpColNam2)){
      tpColID <- which(colnames(smplInfo)==tpColNam) 
      smplInfo$bothTypeCols <- paste0(smplInfo[,tpColNam], " - ", smplInfo[,tpColNam2])
      tpTable <- table(smplInfo$bothTypeCols)
      if(is.na(predefinedMin)){
        largeNr <- min(tpTable)
        if(length(unique(tpTable))>2){
          if(length(unique(tpTable))>5){
            n <- min(tpTable[which(tpTable>7)])
          }else{
            n <- min(setdiff(tpTable, 1))
          }
        }else{
          n=largeNr
        }
      }else{
        n <- predefinedMin
      }
      largeNrSamples <- names(tpTable[which(tpTable>=n)])
      smplInfoL <- smplInfo[which(smplInfo$bothTypeCols %in% largeNrSamples),]
      tpTable <- data.frame(table(smplInfoL$bothTypeCols))
      
      # if(n!=1 && n!=predefinedMin){
      cat("\nNew n = ", n)
      smplInfoSameSize <- smplInfoL[0,]
      for(tp in tpTable[,1]){
        set.seed(which(unique(smplInfoL$bothTypeCols)==tp))
        alreadyIn <- which(smplInfoL$title %in% smplInfoSameSize$title)
        candidateSmpls <- which(smplInfoL$bothTypeCols == tp)
        newCand <- setdiff(candidateSmpls, alreadyIn)
        if(length(newCand)>1){
          newSmpls <- sample(x = newCand, size = n)
        }else{
          newSmpls <- newCand
        }
        smplInfoSameSize <- rbind(smplInfoSameSize, smplInfoL[newSmpls,])
      }
      
      #each type as n and each subtype has ~equally many
      smplInfoSameSize <- smplInfoL[0,]
      for(tp in sort(unique(smplInfoL[,tpColID]))){
        cat("\n ", tp)
        #alreadyIn <- which(smplInfoL$title %in% smplInfoSameSize$title)
        candidateSmpls <- which(smplInfoL[,tpColID] == tp)
        type2TableforType1 <- table(smplInfoL[candidateSmpls,tpColNam2])
        subSmplI <- smplInfoL[0,]
        if(all(type2TableforType1>ceiling(n/(length(type2TableforType1))))){
          newN <- ceiling(n/(length(type2TableforType1)))
          for(tp2 in names(type2TableforType1)){
            set.seed(which(unique(smplInfoL[,tpColID])==tp))
            newSmpls <- sample(x = which(smplInfoL[,tpColNam2]==tp2), size = newN)
            subSmplI <- rbind(subSmplI, smplInfoL[newSmpls,])
          }
          if(!(length(unique(table(subSmplI$SMTSD)))==1 && nrow(subSmplI)==n)){
            newTp2Table <- table(subSmplI$SMTSD)
            tooLargeTypes <- names(newTp2Table)[which(newTp2Table>newN-1)]
            counter <- 1
            while(nrow(subSmplI)>n){
              cat("\n in WHILE")
              tpRowIDs <- which(subSmplI$SMTSD==tooLargeTypes[counter])
              set.seed(which(unique(smplInfoL[,tpColID])==tp))
              subSmplI <- subSmplI[-sample(x = tpRowIDs, size = 1),]
              counter <- counter +1
            }
            smplInfoSameSize <- rbind(smplInfoSameSize, subSmplI)
          }else{
            smplInfoSameSize <- rbind(smplInfoSameSize, subSmplI)
          }
        }else{
          newTp2Table <- table(subSmplI$SMTSD)
          tooLargeTypes <- names(newTp2Table)[which(newTp2Table>newN-1)]
          counter <- 1
          while(nrow(subSmplI)>n){
            cat("\n in WHILE")
            tpRowIDs <- which(subSmplI$SMTSD==tooLargeTypes[counter])
            set.seed(which(unique(smplInfoL[,tpColID])==tp))
            smplInfoL <- smplInfoL[-sample(x = tpRowIDs, size = 1),]
            counter <- counter +1
          }
          smplInfoSameSize <- rbind(smplInfoSameSize, subSmplI)
        }
        # smplInfoSameSize <- rbind(smplInfoSameSize, subSmplI)
      }
      # }
    
    }else{
      # if only same nr of types without time or second sutype
      tpColID <- which(colnames(smplInfo)==tpColNam) 
      tpTable <- table(smplInfo[,tpColNam])
      if(is.na(predefinedMin)){
        largeNr <- min(tpTable)
        if(length(unique(tpTable))>2){
          if(length(unique(tpTable))>5){
            largeNr <- min(tpTable[which(tpTable>3)])
          }else{
            largeNr <- min(setdiff(tpTable, 1))
          }
        }
      }else{
        largeNr <- predefinedMin
      }
      largeNrSamples <- names(tpTable[which(tpTable>=largeNr)])
      smplInfoL <- smplInfo[which(smplInfo[,tpColNam] %in% largeNrSamples),]
      tpTable <- data.frame(tpTable)
      n <- min(tpTable$Freq[which(tpTable$Freq>=largeNr)])
      cat("\nNew n = ", n)
      smplInfoSameSize <- smplInfoL[0,]
      for(tp in tpTable$Var1[which(tpTable$Freq>=n)]){
        alreadyIn <- which(smplInfoL$title %in% smplInfoSameSize$title)
        set.seed(1)
        newSmpls <- sample(setdiff(which(smplInfoL[,tpColNam] == tp), alreadyIn), n)
        smplInfoSameSize <- rbind(smplInfoSameSize, smplInfoL[newSmpls,])
      }
    }
   # }else{
   #   smplInfoSameSize <- smplInfo
   #   tpTable <- table(smplInfo[,tpColNam])
   #   n <- min(tpTable)
   # }
  }
  subx <- dataDf[,match(as.character(smplInfoSameSize$title), colnames(dataDf))]
  ################################################################################
  # nice 
  ################################################################################
  cat('\n', dim(subx))
  cat('\nNormalising Zero Var genes')
  if(ncol(subx)>200){
    cat("\npart1...")
    varInGenes <- apply(subx[1:(nrow(subx)/3),], 1, function(subx) var(subx))
    outGenes <- c(names(varInGenes)[which(varInGenes==0)])
    cat("\npart2...")
    varInGenes <- apply(subx[(nrow(subx)/3 +1):(2*ncol(subx)/3),], 1, function(subx) var(subx))
    outGenes <- c(outGenes, c(names(varInGenes)[which(varInGenes==0)]))
    cat("\npart3...")
    varInGenes <- apply(subx[(2*ncol(subx)/3 + 1 ):nrow(subx),], 1, function(subx) var(subx))
    outGenes <- c(outGenes, c(names(varInGenes)[which(varInGenes==0)]))
  }else{
    varInGenes <- apply(subx, 1, function(subx) var(subx))
    outGenes <- c(names(varInGenes)[which(varInGenes==0)])
  }
  if(length(outGenes)>0) subx <- subx[-which(rownames(subx) %in% outGenes),]
  gc()
  cat('\n', dim(subx))
  
  return(list(subx, smplInfoSameSize))
}


prettyTables <- function(dataDf, smplInfoSameSize, outDir, tpColNam, tpColNam2, name, reOrder, sameNrBothTypes){#}, SUM1){
  
  if(ncol(dataDf)>nrow(dataDf)){
    buSampNames <- rownames(dataDf)
    buGeneNames <- colnames(dataDf)
    dataDf <- data.frame(t(dataDf))
    rownames(dataDf) <- buGeneNames
    colnames(dataDf) <- buSampNames
  }

  # if(reOrder){
  #   if("organism_partDetail" %in% colnames(smplInfoSameSize)) smplInfoSameSize <- smplInfoSameSize[order(as.character(smplInfoSameSize$organism_partDetail)),]
  #   dataDf <- dataDf[,match(colnames(dataDf), as.character(smplInfoSameSize$title))]
  #   
  #   smplInfoSameSize <- smplInfoSameSize[order(as.character(smplInfoSameSize[,tpColID])),]
  #   if(any(colnames(smplInfo)=="detailed")) smplInfoSameSize <- smplInfoSameSize[order(smplInfoSameSize$detailed),]
  #   if(grepl("GSE67310", name)) smplInfoSameSize <- smplInfoSameSize[order(match(smplInfoSameSize$type, c("iN2", "iN1", "iN3", "iN6", "iN7", "iN5", "iN4"))),]
  #   if(grepl("GSE67310", name)){
  #     smplInfoSameSize <- smplInfoSameSize[order(smplInfoSameSize$maintype),]
  #     smplInfoSameSize <- smplInfoSameSize[order(smplInfoSameSize[,tpColNam]),]
  #   }
  #   if(grepl("GSE90553", name)){
  #     smplInfoSameSize <- smplInfoSameSize[order(as.numeric(smplInfoSameSize$time)),]
  #     smplInfoSameSize <- smplInfoSameSize[order(smplInfoSameSize$type),]
  #   }
  #   if(grepl("GSE112004", name)){
  #     smplInfoSameSize <- smplInfoSameSize[order(smplInfoSameSize$generalType),]
  #     smplInfoSameSize <- smplInfoSameSize[order(as.numeric(smplInfoSameSize$time)),]
  #     smplInfoSameSize <- smplInfoSameSize[order(smplInfoSameSize$generalType),]
  #   }
  #   
  #   if(!is.na(tpColNam2)) smplInfoSameSize <- smplInfoSameSize[order(smplInfoSameSize[,tpColNam2]),]
  #   smplInfoSameSize <- smplInfoSameSize[order(smplInfoSameSize[,tpColNam]),]
  #   
  # }
  # 
  if(reOrder){
    if(!is.na(tpColNam2)){
      smplInfoSameSize <- smplInfoSameSize[order(smplInfoSameSize[,tpColNam2]),]
      smplInfoSameSize <- smplInfoSameSize[order(smplInfoSameSize[,tpColNam]),]
    }else{
      smplInfoSameSize <- smplInfoSameSize[order(smplInfoSameSize[,tpColNam]),]
    }
    dataDf <- dataDf[,match(colnames(dataDf), as.character(smplInfoSameSize$title))]
  }
  
  if(is.na(tpColNam2) || !sameNrBothTypes){
    tpTable <- data.frame(table(smplInfoSameSize[,tpColNam]))
    n <- min(tpTable$Freq)
    m <- nrow(tpTable)
  }else{
    bothTypeCols <- paste0(smplInfoSameSize[,tpColNam], " - ", smplInfoSameSize[,tpColNam2])
    tpTable <- data.frame(table(bothTypeCols))
    n <- min(tpTable[,2])
    m <- nrow(tpTable)
  }
  
  #smplInfoSameSize <- smplInfoSameSize[order(as.character(smplInfoSameSize[,tpColNam])),]
  dataDf <- dataDf[,match(as.character(smplInfoSameSize$title), colnames(dataDf))]
  geneNames <- rownames(dataDf)
  
  colClass <- sapply(1:ncol(smplInfoSameSize),function(x) class(smplInfoSameSize[,x]))
  if(any(colClass=="logical")){
    for(col in which(colClass=="logical")){
      smplInfoSameSize[,col] <- ifelse(smplInfoSameSize[,col]==TRUE, 1, 0)
    }
  }
  ###################################################################################
  # save
  ###################################################################################
  if(all(colnames(dataDf)==as.character(smplInfoSameSize$title)) && all(dim(na.omit(dataDf))==dim(dataDf))){
    mainNam <- paste0(n, "x", m, name)
    cat(paste0("\nWriting to file: ",  paste0(mainNam,".txt")))
    cat(paste0("\nFile dimension: ",  ncol(dataDf), " x ", nrow(dataDf)))
    setwd(outDir)
    write.table(data.frame(geneNames), paste0(mainNam,"GeneNamesPREPPED.txt"), quote = F, row.names = F,  sep="\t")
    write.table(t(dataDf), paste0(mainNam,"DatPREPPED.txt"), quote = F, row.names = F,  sep="\t")
    write.table(smplInfoSameSize, paste0(mainNam, "SampleInfoPREPPED.txt"), quote = F, row.names = F,  sep="\t")
    write.table(c(tpColNam, tpColNam2), paste0(mainNam, "TypeColsPREPPED.txt"), quote = F, row.names = F,  sep="\t")
  }
}

