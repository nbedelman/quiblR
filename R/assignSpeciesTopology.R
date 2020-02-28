#' mark the topologies that are consistent with an input species tree, and add as a column to a quibl output data frame
#'
#' @param quiblOutput QuIBL output file
#' @param speciesTree species tree in phylo format
#' @export
#' @examples
#' assignSpeciesTopology(quiblOutput,speciesTree)

assignSpeciesTopology <- function(quiblOutput,speciesTree){
  return(mutate(quiblOutput, isSpeciesTopology=apply(quiblOutput,1,isSpeciesTree, sTree=speciesTree)))
}

#' takes as input one row of a quibl output data frame, as well as a species tree
#' returns TRUE if the row is concordant with the species tree; FALSE otherwise.
#'
#' @param quiblOutput QuIBL output file
#' @param speciesTree species tree in phylo format
#' @export
#' @examples
#' isSpeciesTree(quiblRow,speciesTree)
#'

isSpeciesTree <- function(quiblRow,sTree){
  triplet <- unlist(strsplit(as.character(quiblRow[1]),"_"))
  testTree <- extractTripletTree(sTree,triplet)
  return(ape::is.monophyletic(testTree,setdiff(triplet,as.character(quiblRow[2]))))
}

#' @export
#'

isDiscordant <- function(quiblRow,sTree){
  triplet <- unlist(strsplit(as.character(quiblRow[1]),"_"))
  testTree <- extractTripletTree(sTree,as.character(quiblRow[1]))
  return(! ape::is.monophyletic(testTree,setdiff(triplet,as.character(quiblRow[2]))))
}

#' @export
#'
#'
extractTripletTree <- function(fullTree,triplet){
  return(ape::drop.tip(fullTree,setdiff(fullTree$tip.label,triplet)))
}
