#' Visualizing the ordination plot for a selected region
#'
#' @description
#' Function used to display the ordination plot of samples relative to specific genomic regions.
#' @importFrom purrr map reduce
#' @importFrom dplyr filter select full_join
#' @importFrom vegan cca ordiellipse
#' @param samples A list object containing the epiallele composition matrices from all the samples of the dataset.
#' @param region A string containing the regionID wanted to perform the analysis.
#' @param metadata A dataframe object containing samples metadata. Dataframe should contain dedicated columns for samples IDs and the one indicating the group they belong to.
#' @param printData Logical indicating whether the epiallele matrix should be printed in the standard output or not.
#' @param rmUnmeth Logical indicating whether unmethylated epialleles should be discarded from the analysis.
#' @return A list object storing the ordination results and an ordination plot.
#' @export
#' @examples
#' samples_list <- list(Sample1_epiAnalysis.txt,
#'                      Sample2_epiAnalysis.txt,
#'                      Sample3_epiAnalysis.txt,
#'                      Sample4_epiAnalysis.txt)
#'
#' p <- runCCA(samples = samples_list,
#'             region = "chr1_34567876_34567923",
#'             metadata = ann,
#'             printData = FALSE,
#'             rmUnmeth = FALSE)
#' p

runCCA <- function(samples, region, metadata, printData = FALSE, rmUnmeth = FALSE){
  data <- getEpimatrix(samples, region)
  if(rmUnmeth == TRUE){
    data <- data[grep("1", colnames(data))]
    data = data %>% dplyr::filter(!rowSums(.) == 0)
  }
  if(length(colnames(data)) <= 1){
    print("Plotting is not possible with just one epiallele specie")
  } else {
    metadata = metadata %>%
      dplyr::filter(Samples %in% rownames(data))
    mod <- vegan::cca(data ~ Group, metadata)
    plot(mod, type="n", display = "sites")
    with(metadata, text(mod, display="sites", labels = as.character(Group),
                        col=as.numeric(Group)))
    pl <- with(metadata, vegan::ordiellipse(mod, Group, kind = "se", conf = 0.95,
                                            lwd = 2, draw = "polygon", col = 1:4,
                                            border = 1:4, alpha = 63))
    if(printData == TRUE){
      print(data)
    }
    return(pl)
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
