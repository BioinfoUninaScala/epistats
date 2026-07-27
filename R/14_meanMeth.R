#' Calculating the average DNA methylation for one interval
#'
#' @param matrix Epialleles binary matrix.
#' @return A numeric indicating the average DNA methylation of one genomic region
#' @export

meanMeth <- function(matrix){
  mean_met=round(sum(as.matrix(matrix),na.rm=T)/(dim(matrix)[1]*dim(matrix)[2]),2)
  return(mean_met)
}
