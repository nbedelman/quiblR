#' make a matrix which scores the probability of introgression for each non-species tree pairing for each triplet
#'
#' @param quiblOut quibl output data
#' @param speciesTree phylo object of the species tree
#' @param summaryType either "max" or "mean"; describes how to deal with multiple instances of the same pair of taxa
#' @export
#'

getSummaryMatrix <- function(quiblOut,speciesTree,summaryType="mean"){
correlationDict <- hash()
for (row in 1:nrow(quiblOut)){
  numTrees <- sum(quiblOut[which(quiblOut$triplet==unique(quiblOut$triplet)[1]),]$count)
  taxa <- setdiff(unlist(strsplit(as.character(quiblOut[row,][,"triplet"]),"_")),quiblOut[row,][,"outgroup"])
  key <- paste(sort(taxa)[1], sort(taxa)[2],sep = "_")
  if(! isSpeciesTree(quiblOut[row,],speciesTree)){
    if (key %in% keys(correlationDict)){
      correlationDict[[key]] <- append(correlationDict[[key]],quiblOut[row,][,"mixprop2"]*quiblOut[row,][,"count"]/numTrees)
    } else {
      correlationDict[[key]] <- c(quiblOut[row,][,"mixprop2"]*quiblOut[row,][,"count"]/numTrees)
    }
  }
}

# mat <- matrix(NA, nrow=length(unique(quiblOut[,"outgroup"])),
#               ncol=length(unique(quiblOut[,"outgroup"])),
#               dimnames=list(unique(quiblOut[,"outgroup"]),unique(quiblOut[,"outgroup"])))

inp <- data.frame(tax1=character(),tax2=character(), value=numeric())
for(k in keys(correlationDict)){
  inp <- rbind(inp, data.frame(tax1=unlist(strsplit(k,"_"))[1],
                               tax2=unlist(strsplit(k,"_"))[2],
                               value=mean(correlationDict[[k]])))
}

# inp.mtx <- as.matrix(inp)
# mat[inp.mtx[,1:2] ] <- as.numeric(inp.mtx[,3])

return(inp)
}

