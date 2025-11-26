#' Calculating the average DNA methylation for one interval
#'
#' @param matrix Epialleles binary matrix.
#' @return A numeric indicating the average DNA methylation of one genomic region
#' @export

meanMeth <- function(matrix){
  mean_met=round(sum(as.matrix(matrix),na.rm=T)/(dim(matrix)[1]*dim(matrix)[2]),2)
  return(mean_met)
}


#' Calculating the average DNA methylation for a not binary epialleles matrix.
#'
#' @param matrix Epialleles matrix.
#' @param param A list of two numeric values indicating the codes for methylated and unmethylated status in the epialleles matrix.
#' @return A numeric indicating the average DNA methylation of one genomic region
#' @export

meanMethONT <- function(matrix, param = list(meth_code = 1L, unmeth_code = 0L)){
  meth_code = param[[1]]
  unmeth_code = param[[2]]
  mean_met=round(sum(as.matrix(matrix == meth_code),na.rm=T)/sum(as.matrix(matrix %in% c(meth_code, unmeth_code)),na.rm=T),2)
  return(mean_met)
}

