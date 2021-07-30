#' cdist
#'
#' @param matrix epialleles matrix
#' @return integer indicating mean dist between CpG in the interval
#' @export

cdist <- function(matrix){
  Cpos=as.numeric(names(matrix))
  dist=round(mean(Cpos[-1]-Cpos[-length(Cpos)])/length(Cpos),2)
  return(dist)
}
