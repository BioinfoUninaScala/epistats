#' Reporting the epiallele observed at the highest frequency in one genomic interval
#'
#' @param matrix Epialleles binary matrix.
#' @return A string reporting the epiallele observed at the highest frequency in the genomic interval.
#' @export

maxfreq <- function(matrix){
  matrix$epi= apply(matrix, 1, function(x) paste(x, collapse = ""))
  matrix= as.data.frame(table(matrix$epi))
  maxfreq=paste(matrix[matrix$Freq==max(matrix$Freq),1],sep=",")
  return(maxfreq)
}
