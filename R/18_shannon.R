#' Calcuting the Shannon Entropy in one genomic region.
#'
#' @param matrix Epialleles binary matrix.
#' @return A numeric indicating Shannon Entropy value relative to one genomic interval.
#' @export

shannon <- function(matrix){
  matrix$epi= apply(matrix, 1, function(x) paste(x, collapse = ""))
  matrix= as.data.frame(table(matrix$epi))
  prop= matrix$Freq/sum(matrix$Freq)
  shannon=round(-sum(prop * log2(prop), na.rm=T),2)
  return(shannon)
}
