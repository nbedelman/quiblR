#' Read output file from QuIBL
#'
#' This function reads in the basic output generated by a QuIBL run
#' @param filepath QuIBL output file
#' @export
#' @examples
#' read.csv.quibl(quiblFile.csv)

read_csv_quibl <- function(filepath){
  utils::read.csv(opt$inFile,header=T,
           colClasses = c(rep("factor",2),rep("numeric",9)))
}


