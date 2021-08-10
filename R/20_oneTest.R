#' oneTest
#'
#' @importFrom purrr map reduce
#' @importFrom dplyr select full_join filter group_by summarize
#' @param sample_list list of samples epiMatrices coming from epiAnalysis function
#' @param region genomic coordinate of your region of interest
#' @param metadata dataframe. Your samples metadata
#' @return Dataframe with statistics
#' @export

oneTest <- function(sample_list, region, metadata, rmUnmeth = FALSE, printData = FALSE){
  data <- getEpimatrix(sample_list, region)
  if(rmUnmeth == TRUE){
    data <- data[grep("1", colnames(data))]
    data = data %>% dplyr::filter(!rowSums(.) == 0)
  }
  dist = vegan::vegdist(data, method="bray")
  metadata <- metadata %>%
           dplyr::filter(Samples %in% rownames(data))
  check <- metadata %>%
              dplyr::group_by(Group) %>%
                 dplyr::summarize(count = length(Group))
  check$reps <- paste(check$Group, check$count, sep = ":")
  if(length(unique(metadata$Group)) > 1){
    ## problema formula
    p <- suppressMessages(adonis2(dist ~ Group, data = metadata))
    result= data.frame("Region"= region,
                       "F.statistics"= p$F[1],
                       "p.value" = p$`Pr(>F)`[1],
                       "num_EpiSpecies" = length(colnames(data)),
                       "num_Reps" = paste(check$reps, collapse = ","))
  } else {
    result = data.frame("Region"= region,
                        "F.statistics"= NA,
                        "p.value" = NA,
                        "num_EpiSpecies" = NA,
                        "num_Reps" = NA)
  }
  if(printData == TRUE){
    print(data)
  }
  return(result)
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
