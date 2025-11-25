#' Calculating the average DNA methylation for one interval
#'
#' @param matrix Epialleles binary matrix.
#' @return A numeric indicating the average DNA methylation of one genomic region
#' @export

meanMeth <- function(matrix){
  mean_met=round(sum(as.matrix(matrix),na.rm=T)/(dim(matrix)[1]*dim(matrix)[2]),2)
  return(mean_met)
}


#' Calculating the average DNA methylation for one interval for ONT data
#'
#' @param matrix Epialleles matrix.
#' @param methylated_code A numeric value indicating the methylation status in the epialelles matrix.
#' @return A numeric indicating the average DNA methylation of one genomic region
#' @export

meanMethONT <- function(matrix, methylated_code = 1L){
  mean_met=round(sum(as.matrix(matrix == methylated_code),na.rm=T)/(dim(matrix)[1]*dim(matrix)[2]),2)
  return(mean_met)
}
