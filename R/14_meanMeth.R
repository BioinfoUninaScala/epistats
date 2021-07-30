#' meanMeth
#'
#' @param matrix epialleles matrix
#' @return numeric indicating the mean meth of the entire region
#' @export

meanMeth <- function(matrix){
  mean_met=round(sum(as.matrix(matrix),na.rm=T)/(dim(matrix)[1]*dim(matrix)[2]),2)
  return(mean_met)
}
