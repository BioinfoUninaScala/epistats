#' epidiv
#'
#' @importFrom purrr map reduce
#' @importFrom dplyr filter select full_join
#' @importFrom vegan simper
#' @param samples list of epimatrices
#' @param region one region of interest 
#' @param metadata samples metadata 
#' @return list 
#' @export


epidiv <- function(samples, region, metadata){
  data <- mkmatr(samples, region)
  metadata = metadata %>% 
    dplyr::filter(., Samples == rownames(data))
  sim <- with(metadata, vegan::simper(data, Group, permutations = 999))
  return(sim)
}

mkmatr <- function(samples, region){
  filt = purrr::map(samples, ~ dplyr::filter(., id == region))
  filt = purrr::map(filt, ~ dplyr::select(., 1:2)) 
  counts = filt %>% purrr::reduce(dplyr::full_join, by = "Var1")
  rownames(counts)= counts$Var1
  counts$Var1=NULL
  colnames(counts) = names(samples)
  data = as.data.frame(t(counts))
  data[is.na(data)] = 0
  data = data %>% dplyr::filter(., !rowSums(.) == 0)
  return(data)
}


