#' singleton
#'
#' @param matrix epialleles matrix
#' @return integer indicating how many species are represented alone
#' @export

singleton <- function(matrix){
  matrix$epi= apply(matrix, 1, function(x) paste(x, collapse = ""))
  matrix= as.data.frame(table(matrix$epi))
  singleton= length(matrix$Freq[matrix$Freq==1])
  return(singleton)
}
