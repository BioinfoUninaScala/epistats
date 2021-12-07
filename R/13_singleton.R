#' Calculating the number of epialleles with just one observation
#'
#' @param matrix Epiallele binary matrix
#' @return An integer indicating the number of epiallele with a unique observation in a given interval
#' @export

singleton <- function(matrix){
  matrix$epi= apply(matrix, 1, function(x) paste(x, collapse = ""))
  matrix= as.data.frame(table(matrix$epi))
  singleton= length(matrix$Freq[matrix$Freq==1])
  return(singleton)
}
