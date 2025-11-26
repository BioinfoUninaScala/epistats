#' Calcutating the mean distance among the CpG sites contained in the analysed region
#'
#' @param matrix Epiallele binary matrix for one genomic region
#' @return Integer indicating the mean distance between CpG in the interval
#' @export

cdist <- function(matrix){
  Cpos=as.numeric(names(matrix))
  dist=round(mean(Cpos[-1]-Cpos[-length(Cpos)]),2)
  return(dist)
}


#' Calcutating the mean distance among the CpG sites contained in the analysed region for ONT.
#'
#' @param matrix Epiallele matrix for one genomic region
#' @return Integer indicating the mean distance between CpG in the interval
#' @export

cdist_ont <- function(matrix){
  Cpos=as.numeric(colnames(matrix))
  dist=round(mean(Cpos[-1]-Cpos[-length(Cpos)]),2)
  return(dist)
}
