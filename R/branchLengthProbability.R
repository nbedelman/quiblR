#' compute the branch length probability distribution from quibl data
#'
#' @param x branch length
#' @param C internodal distance
#' @param lmbd population size parameter (theta/2)
#' @param mixprop mixture proportion
#' @export

branchLengthProbability <-  function(x,C,lmbd,mixprop){
  if (x>0 & x>=C*lmbd){
    y=0.5/lmbd*exp(-x/lmbd)*(1+exp(C))
  }
  else if (x>0 & C*lmbd>x){
    y=sinh(x/lmbd)/(lmbd*(-1+exp(C)))
  }
  else{
    y=0
  }
  return (y*mixprop)
}

#' compute the branch length probability for an individual data point
#'
#' @param x branch length
#' @param quiblRow row from quibl output
#' @param distribution one of "ILSOnly", "ILSMix", "nonILSMix"
#' @export

getSingleProb <- function(x,quiblRow,distribution){
  if (distribution == "ILSOnly"){
    return (branchLengthProbability(x,C=0,lmbd=as.numeric(quiblRow[8]),mixprop=1))
  } else if (distribution == "ILSMix"){
    return (branchLengthProbability(x,C=0,lmbd=as.numeric(quiblRow[7]),mixprop=as.numeric(quiblRow[5])))
  } else if (distribution == "nonILSMix"){
    return (branchLengthProbability(x,C=as.numeric(quiblRow[4]),lmbd=as.numeric(quiblRow[7]),mixprop=as.numeric(quiblRow[6])))
  } else {print("please specify one of ILSOnly, ILSMix,  or nonILSMix as the distribution")}
    }

#' @export
getILSOnlyDist <- function(xmin,xmax,quiblRow){
  x <- seq(xmin,xmax,length.out=1000)
  y <- unlist(lapply(x,getSingleProb, quiblRow=quiblRow,distribution="ILSOnly"))
  return(data.frame(x=x,y=y))
}

#' @export
getILSMixtureDist <- function(xmin,xmax,quiblRow){
  x <- seq(xmin,xmax,length.out=1000)
  y <- unlist(lapply(x,getSingleProb, quiblRow=quiblRow,distribution="ILSMix"))
  return(data.frame(x=x,y=y))
}

#' @export
getNonILSMixtureDist <- function(xmin,xmax,quiblRow){
  x <- seq(xmin,xmax,length.out=1000)
  y <- unlist(lapply(x,getSingleProb, quiblRow=quiblRow,distribution="nonILSMix"))
  return(data.frame(x=x,y=y))
}

#' @export
getTotalMixtureDist <- function(xmin,xmax,quiblRow){
  ILS <- getILSMixtureDist(xmin,xmax,quiblRow)
  nonILS <- getNonILSMixtureDist(xmin,xmax,quiblRow)
  return(data.frame(x=ILS$x, y=ILS$y+nonILS$y))
}

#' @export
findOutgroup <- function(tripletTree){
  for (tip in tripletTree$tip.label){
    if (ape::is.monophyletic(tripletTree,setdiff(tip,tripletTree$tip.label))){
      return (tip)
    }
  }
}

#' @export
locusStats <- function(tree,tripFrame){
  testTree <- extractTripletTree(tree,unique(tripFrame$triplet))
  outgroup <- findOutgroup(testTree)
  branchLength <- ape::dist.nodes(testTree)[4,5]
  introProb <- getSingleProb(outgroup, tripFrame[which(tripFrame$outgroup==outgroup),],"nonILSMix")
  ILSProb <- getSingleProb(outgroup, tripFrame[which(tripFrame$outgroup==outgroup),],"ILSMix")
  return(data.frame(tree=tree,out=outgroup,branchLength=branchLength,introProb=introProb/(introProb+ILSProb)))
}

#' @export
getPerLocusStats <- function(quiblOutput,triplet,treeList){
  #first, extract the appropriate rows from the full quibl output
  thisTriplet <- subset(quiblOutput,triplet==triplet)
  perLocusOut <- sapply(treeList,locusStats,tripFrame=thisTriplet)
  return(perLocusOut)
}
