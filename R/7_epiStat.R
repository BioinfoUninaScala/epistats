#' epiStat
#'
#' @importFrom purrr map reduce map_dfr map_dbl
#' @importFrom dplyr select distinct mutate full_join filter group_by summarize
#' @importFrom future plan
#' @importFrom furrr future_map2_dfr
#' @importFrom foreach foreach
#' @importFrom stats complete.cases p.adjust
#' @importFrom vegan vegdist adonis2
#' @param sample_list list of epiMatrix coming from epiAnalysis
#' @param metadata samples metadata
#' @param cores num of cores
#' @return dataframe w/ stats
#' @export


epiStat <- function(sample_list, metadata, rmUnmeth = FALSE, cores = 1, minGroups = 2, minSampleSize = 2){
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
  filt_samples <- foreach::foreach(reg = regs) %do% {
    filt = purrr::map(sample_list, ~ dplyr::filter(., id %in% reg))
  }
  ## Plan parallel
  future::plan(multisession, workers = cores)
  ## Remove warning message
  options(future.rng.onMisuse = "ignore")
  ## Call allStat for all the blocks
  prova = furrr::future_map2_dfr(regs, filt_samples, allStat, metadata, rmUnmeth, minGroups, minSampleSize)
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

oneStat <- function(sample_list, region, metadata, rmUnmeth, minGroups, minSampleSize){
  data <- getEpimatrix(sample_list, region)
  if(rmUnmeth == TRUE){
    data <- data[grep("1", colnames(data))]
    data <- data %>% dplyr::filter(!rowSums(.) == 0)
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
    if(length(unique(metadata$Group)) > minGroups){
      a = groups %>% purrr::map_dbl(~ length(metadata$Samples[metadata$Group == .x]))
      a = t(as.data.frame(a, row.names = groups))
      rownames(a) = NULL
      p <- suppressMessages(adonis2(dist ~ Group, data = metadata))
      result = data.frame("Region" = region,
                         "F.statistics" = p$F[1],
                         "p.value" = p$`Pr(>F)`[1],
                         "num_EpiSpecies" = length(colnames(data)))
      result = cbind(result, a)
    } else {
      b = rep(NA, length(groups))
      b = t(as.data.frame(b, row.names = groups))
      rownames(b) = NULL
      result = data.frame("Region" = region,
                          "F.statistics" = NA,
                          "p.value" = NA,
                          "num_EpiSpecies" = NA)
      result = cbind(result,b)
    }
  return(result)
}


allStat <- function(regions, sample_list, metadata, rmUnmeth, minGroups, minSampleSize){
  out = regions %>% map_dfr(~ oneStat(sample_list, .x, metadata, rmUnmeth, minGroups, minSampleSize))
  return(out)
}
