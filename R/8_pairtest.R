#' pairtest
#'
#' @importFrom purrr map reduce
#' @importFrom dplyr filter select full_join
#' @importFrom pairwiseAdonis pairwise.adonis2
#' @param samples list of epimatrices
#' @param region one region of interest
#' @param metadata samples metadata
#' @return dataframe
#' @export


pairtest <- function(samples, region, metadata){
  data <- getEpimatrix(samples, region)
  metadata = metadata %>%
    dplyr::filter(Samples %in% rownames(data))
  data$Group = metadata$Group
  if (length(metadata$Group) > 2){
    p <- suppressMessages(pairwiseAdonis::pairwise.adonis2(data[,-ncol(data)] ~ Group,
                                          data = data,
                                          nperm = 999))
    o = purrr::map(2:length(p), ~ p[[.]]$`Pr(>F)`[1])
    o = as.data.frame(o)
    names(o) = names(p)[-1]
    rownames(o) = region
  }
  return(o)
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
