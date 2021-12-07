#' Finding the epiallele species which are contributing most to the dissimilarity among groups
#'
#' @description
#' The function returns the information regarding which epialleles are the ones contributing for the dissimilarity of samples among groups.
#' @importFrom purrr map reduce
#' @importFrom dplyr filter select full_join
#' @importFrom vegan simper
#' @param samples A list object containing the epiallele composition matrices from all the samples of the dataset.
#' @param region A string containing the regionID wanted to perform the analysis.
#' @param metadata A dataframe object containing samples metadata. Dataframe should contain dedicated columns for samples IDs and the one indicating the group they belong to.
#' @return A list containing the dissimilarity results.
#' @export
#' @examples
#' samples_list <- list(Sample1_epiAnalysis.txt,
#'                      Sample2_epiAnalysis.txt,
#'                      Sample3_epiAnalysis.txt,
#'                      Sample4_epiAnalysis.txt)
#'
#' epidiv <- epidiv(samples = samples_list,
#'                  region = "chr1_34567876_34567923",
#'                  metadata = ann)
#'
#'
#'


epidiv <- function(samples, region, metadata){
  data <- getEpimatrix(samples, region)
  metadata = metadata %>%
    dplyr::filter(Samples %in% rownames(data))
  sim <- with(metadata, suppressMessages(vegan::simper(data, Group, permutations = 999)))
  return(sim)
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


