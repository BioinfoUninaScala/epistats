#' shannon
#'
#' @param matrix epialleles matrix
#' @return numeric indicating shannon entropy value
#' @export

shannon <- function(matrix){
  matrix$epi= apply(matrix, 1, function(x) paste(x, collapse = ""))
  matrix= as.data.frame(table(matrix$epi))
  prop= matrix$Freq/sum(matrix$Freq)
  shannon=round(-sum(prop * log2(prop), na.rm=T),2)
  return(shannon)
}
