#' runCCA
#'
#' @importFrom purrr map reduce
#' @importFrom dplyr filter select full_join
#' @importFrom vegan cca ordiellipse
#' @param samples list of epimatrices
#' @param region one region of interest
#' @param metadata samples metadata
#' @return list
#' @export


runCCA <- function(samples, region, metadata){
  data <- getEpimatrix(samples, region)
  metadata = metadata %>%
    dplyr::filter(Samples %in% rownames(data))
  mod <- vegan::cca(data ~ Group, metadata)
  plot(mod, type="n", display = "sites")
  with(metadata, text(mod, display="sites", labels = as.character(Group),
                      col=as.numeric(Group)))
  pl <- with(metadata, vegan::ordiellipse(mod, Group, kind = "se", conf = 0.95,
                                          lwd = 2, draw = "polygon", col = 1:4,
                                          border = 1:4, alpha = 63))
  return(pl)
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
