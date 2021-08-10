#' oneTest
#'
#' @importFrom purrr map reduce map_dbl
#' @importFrom dplyr select full_join filter group_by summarize
#' @param sample_list list of samples epiMatrices coming from epiAnalysis function
#' @param region genomic coordinate of your region of interest
#' @param metadata dataframe. Your samples metadata
#' @return Dataframe with statistics
#' @export

oneTest <- function(sample_list, region, metadata, rmUnmeth = FALSE, minGroups = 2, minSampleSize = 2, printData = FALSE){
  data <- getEpimatrix(sample_list, region)
  if(rmUnmeth == TRUE){
    data <- data[grep("1", colnames(data))]
    data = data %>% dplyr::filter(!rowSums(.) == 0)
  }
  dist = vegan::vegdist(data, method="bray")
  groups = unique(metadata$Group)
  metadata <- metadata %>%
           dplyr::filter(Samples %in% rownames(data))
  check <- metadata %>%
              dplyr::group_by(Group) %>%
                 dplyr::summarize(count = length(Group))
  check <- check[check$count >= minSampleSize,]
  metadata <- metadata[metadata$Group %in% check$Group,]
  data <- data[rownames(data) %in% metadata$Samples,]
  if(length(unique(metadata$Group)) >= minGroups){
    a = groups %>% purrr::map_dbl(~ length(metadata$Samples[metadata$Group == .x]))
    a = t(as.data.frame(a, row.names = groups))
    rownames(a) = NULL
    ## problema formula
    p <- suppressMessages(adonis2(dist ~ Group, data = metadata))
    result = data.frame("Region"= region,
                       "F.statistics"= p$F[1],
                       "p.value" = p$`Pr(>F)`[1],
                       "num_EpiSpecies" = length(colnames(data)))
    result = cbind(result, a)
  } else {
    b = rep(NA, length(groups))
    b = t(as.data.frame(b, row.names = groups))
    rownames(b) = NULL
    result = data.frame("Region"= region,
                        "F.statistics"= NA,
                        "p.value" = NA,
                        "num_EpiSpecies" = NA)
    result = cbind(result, b)
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
