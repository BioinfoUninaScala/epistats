#' numCG
#'
#' @param matrix epialleles matrix
#' @return integer indicating how many CpGs are in a given region
#' @export

numCG <- function(matrix){
  num_cg=ncol(matrix)
  return(num_cg)
}
