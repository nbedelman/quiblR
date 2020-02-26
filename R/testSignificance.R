#' test for significance at a given delta BIC threshold
#'
#' very basic function that just returns true if the difference between BIC1Dist and BIC2Dist is larger than some threshold
#' @param quiblRow
#' @param threshold
#' @export
#' @examples
#' testSignificance(quiblRow,threshold)

testSignificance <- function(quiblRow,threshold){
  BICdiff <- as.numeric(quiblRow[10]) - as.numeric(quiblRow[9])
  return(BICdiff > threshold)
}
