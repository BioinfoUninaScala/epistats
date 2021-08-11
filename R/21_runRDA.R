#' runRDA
#'
#' @importFrom purrr map reduce
#' @importFrom dplyr filter select full_join
#' @importFrom vegan rda
#' @param samples list of epimatrices
#' @param region one region of interest
#' @return PCA plot
#' @export

runRDA <- function(samples, region, printData = FALSE){
  data <- getEpimatrix(samples, region)
  if(length(colnames(data)) <= 1){
    print("Plotting is not possible with just one epiallele specie")
  } else {
    myrda <- vegan::rda(data)
    p <- biplot(myrda, scaling = "symmetric")
    if(printData == TRUE){
      print(data)
    }
    return(p)
  }
}

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
