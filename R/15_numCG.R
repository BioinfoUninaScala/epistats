#' Counting the number of CpG sites contained in one interval
#'
#' @param matrix Epialleles binary matrix.
#' @return An integer indicating the number of CpGs displayed for one genomic interval.
#' @export

numCG <- function(matrix){
  num_cg=ncol(matrix)
  return(num_cg)
}
