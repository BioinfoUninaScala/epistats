#' Counting the number of reads mapping on one genomic interval.
#'
#' @param matrix Epialleles binary matrix.
#' @return An integer indicating the coverage observed in the genomic region.
#' @export

num_reads= function(matrix){
  num_reads=nrow(matrix)
  return(num_reads)
}
