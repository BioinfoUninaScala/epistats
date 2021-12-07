#' Calculating the number of different epialleles observed for one genomic region
#'
#' @param matrix Epiallele binary matrix.
#' @return A numeric indicating the number of observed epialleles in one interval.
#' @export

epi <- function(matrix){
  matrix$epi= apply(matrix, 1, function(x) paste(x, collapse = ""))
  matrix= as.data.frame(table(matrix$epi))
  epi=nrow(matrix)
  return(epi)
}
