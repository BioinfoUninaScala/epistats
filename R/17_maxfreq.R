#' maxfreq
#'
#' @param matrix epialleles matrix
#' @return character
#' @export

maxfreq <- function(matrix){
  matrix$epi= apply(matrix, 1, function(x) paste(x, collapse = ""))
  matrix= as.data.frame(table(matrix$epi))
  maxfreq=paste(matrix[matrix$Freq==max(matrix$Freq),1],sep=",")
  return(maxfreq)
}
