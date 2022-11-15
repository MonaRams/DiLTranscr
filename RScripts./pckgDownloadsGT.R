pckgs <- c("BiocManager", "readr")
for(pckg in pckgs){
  if(!require(pckg, quietly = TRUE)){
    cat(paste0("\nInstalling package ",pckg))
    install.packages(pckg)
  }
}

#################################################################
bioCondPckgs <-  "GEOquery"
for(bioPckg in bioCondPckgs){
  BiocManager::install(c(bioPckg))
}
