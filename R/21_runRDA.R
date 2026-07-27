#' Running principal component analysis for one genomic interval
#'
#' @importFrom purrr map reduce
#' @importFrom dplyr filter select full_join
#' @importFrom vegan rda
#' @param samples A list object containing the epiallele composition matrices from all the samples of the dataset.
#' @param region A character indicating the regionID wanted to perform the analysis.
#' @param printData Logical indicating whether the epiallele matrix should be printed in the standard output or not (Default = FALSE).
#' @param rmUnmeth Logical indicating if 0 methylated epialleles should be removed or not from the analysis.
#' @return Principal Component analysis results and plot.
#' @export
#' @examples
#' samples_list <- list(Sample1_epiAnalysis.txt,
#'                      Sample2_epiAnalysis.txt,
#'                      Sample3_epiAnalysis.txt,
#'                      Sample4_epiAnalysis.txt)
#'
#' runRDA(samples = samples_list,
#'        region = "chr1_34567876_34567923",
#'        printData = FALSE,
#'        rmUnmeth = FALSE)

runRDA <- function(samples, region, printData = FALSE, rmUnmeth = FALSE){
  data <- getEpimatrix(samples, region)
  if(rmUnmeth == TRUE){
    data <- data[grep("1", colnames(data))]
    data = data %>% dplyr::filter(!rowSums(.) == 0)
  }
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
