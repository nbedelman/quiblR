#' make a matrix which scores the probability of introgression for each non-species tree pairing for each triplet
#'
#' @param quiblOut quibl output data
#' @param speciesTree phylo object of the species tree
#' @param summaryType either "max" or "mean"; describes how to deal with multiple instances of the same pair of taxa
#' @export
#'

getIntrogressionSummary <- function(quiblOut,speciesTree,summaryType="mean"){
  species <- unique(quiblOut[,"outgroup"])
  orderedSpecies <- factor(species, levels=speciesTree$tip.label)
  correlationDict <- makeEmptyHash(orderedSpecies)

for (row in 1:nrow(quiblOut)){
  numTrees <- sum(quiblOut[which(quiblOut$triplet==unique(quiblOut$triplet)[1]),]$count)
  taxa <- setdiff(unlist(strsplit(as.character(quiblOut[row,][,"triplet"]),"_")),quiblOut[row,][,"outgroup"])
  key1 <- paste(sort(taxa)[1], sort(taxa)[2],sep = "_")
  key2 <- paste(sort(taxa)[2], sort(taxa)[1],sep = "_")
  if(! isSpeciesTree(quiblOut[row,],speciesTree)){
    if (key1 %in% keys(correlationDict)){
      correlationDict[[key1]] <- append(correlationDict[[key1]],quiblOut[row,][,"mixprop2"]*quiblOut[row,][,"count"]/numTrees)
    } else {
      correlationDict[[key2]] <- append(correlationDict[[key2]],quiblOut[row,][,"mixprop2"]*quiblOut[row,][,"count"]/numTrees)
    }
  }
}

inp <- data.frame(tax1=character(),tax2=character(), value=numeric())
for(k in keys(correlationDict)){
  taxa <- unlist(strsplit(k,"_"))[1:2]
  inp <- rbind(inp, data.frame(tax1=taxa,
                               tax2=rev(taxa),
                               value=rep(mean(correlationDict[[k]]),2)))
}
inp[,"tax1"] <- factor(inp[,"tax1"], levels = speciesTree$tip.label)
inp[,"tax2"] <- factor(inp[,"tax2"], levels = speciesTree$tip.label)
inp[,"value"][is.na(inp[,"value"])] <- 0
return(inp)
}

makeEmptyHash <- function(labels){
  out <- hash()
  for (lab1 in 1:(length(labels))){
    for (lab2 in (lab1):length(labels)){
      thisLab <- paste(labels[lab1],labels[lab2],sep="_")
      out[[thisLab]] <- numeric()
    }
  }
  return(out)
}



