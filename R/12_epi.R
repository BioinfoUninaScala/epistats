#' epi
#'
#' @param matrix epialleles matrix
#' @return character indicating the max frequent epiallelic specie
#' @export

epi <- function(matrix){
  matrix$epi= apply(matrix, 1, function(x) paste(x, collapse = ""))
  matrix= as.data.frame(table(matrix$epi))
  epi=nrow(matrix)
  return(epi)
}
