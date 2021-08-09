#' epiStat
#'
#' @importFrom purrr map reduce map_dfr
#' @importFrom dplyr select distinct mutate full_join filter
#' @importFrom future plan
#' @importFrom furrr future_map_dfr
#' @importFrom stats complete.cases p.adjust
#' @importFrom vegan vegdist adonis2
#' @param sample_list list of epiMatrix coming from epiAnalysis
#' @param metadata samples metadata
#' @param cores num of cores
#' @return dataframe w/ stats
#' @export


epiStat <- function(sample_list, metadata, rmUnmeth = FALSE, cores = 1){
  list = purrr::map(sample_list, function(x) x %>%
                      dplyr::select(id) %>%
                      dplyr::distinct() %>%
                      dplyr::mutate(present = TRUE))
  ## Names
  filt = list %>%
    purrr::reduce(dplyr::full_join, by = "id") %>%
    dplyr::mutate(count_na = rowSums(is.na(.))) %>%
    dplyr::filter(count_na < (length(sample_list)-2))
  ## Take regions IDs
  regions = filt$id
  ## Split regions
  regs <- split(regions, (seq(length(regions))-1) %/% (length(regions)/cores))
  ### Split samples
  filt_samples <- foreach(reg = regs) %do% {
    filt = purrr::map(sample_list, ~ dplyr::filter(., id %in% reg))
  }
  ## Plan parallel
  future::plan(multisession, workers = cores)
  ## Remove warning message
  options(future.rng.onMisuse = "ignore")
  ## Call allStat for all the blocks
  prova = furrr::future_map2_dfr(regs, filt_samples, allStat, metadata, rmUnmeth)
  ## Remove NAs
  prova = prova[stats::complete.cases(prova)==TRUE,]
  ## Adjust p-values
  prova$p.adjust = stats::p.adjust(prova$p.value, method = "fdr")
  return(prova)
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

oneStat <- function(sample_list, region, metadata, rmUnmeth){
  data <- getEpimatrix(sample_list, region)
  if(rmUnmeth == TRUE){
    data <- data[grep("1", colnames(data))]
  }
  dist = vegan::vegdist(data, method="bray")
  metadata = metadata %>%
    filter(Samples %in% rownames(data))
  if(length(unique(metadata$Group)) > 1){
    ## problema formula
    p <- suppressMessages(adonis2(dist ~ Group, data = metadata))
    result= data.frame("Region"= region,
                       "F.statistics"= p$F[1],
                       "p.value" = p$`Pr(>F)`[1])
  } else {
    result = data.frame("Region"= region,
                        "F.statistics"= NA,
                        "p.value" = NA)
  }
  return(result)
}


allStat <- function(regions, sample_list, metadata, rmUnmeth){
  out = regions %>% map_dfr(~ oneStat(sample_list, .x, metadata, rmUnmeth))
  return(out)
}
