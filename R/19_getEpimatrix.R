#' Getting the epialleles matrix for one genomic region
#'
#' @description
#' The function gives as output the epiallele composition matrix relative to one genomic interval.
#' @importFrom purrr map reduce
#' @importFrom dplyr select full_join filter
#' @param samples A list object containing the epiallele composition matrices from all the samples of the dataset.
#' @param region A character indicating the regionID wanted to perform the analysis.
#' @return Epiallele community data matrix.
#' @export

getEpimatrix <- function(samples, region){
  filt <- purrr::map(samples, ~ dplyr::filter(., id == region))
  filt <- filt %>% purrr::map(select, 1:2)
  counts <- filt %>% purrr::reduce(dplyr::full_join, by = "Var1")
  epinames = counts$Var1
  counts$Var1 = NULL
  colnames(counts) = names(samples)
  data = as.data.frame(t(counts))
  colnames(data) = epinames
  data[is.na(data)] = 0
  data = data %>% dplyr::filter(!rowSums(.) == 0)
  return(data)
}
