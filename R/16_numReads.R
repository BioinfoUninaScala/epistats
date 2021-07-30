#' numReads
#'
#' @param matrix epialleles matrix
#' @return integer indicating the depth of a given region
#' @export

num_reads= function(matrix){
  num_reads=nrow(matrix)
  return(num_reads)
}
