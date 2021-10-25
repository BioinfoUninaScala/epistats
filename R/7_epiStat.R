#' epiStat
#'
#' @importFrom purrr map reduce map_dfr map_dbl map2 map_df
#' @importFrom dplyr select mutate full_join filter group_by summarize group_by ungroup count n_distinct
#' @importFrom foreach foreach %dopar%
#' @importFrom stats complete.cases p.adjust
#' @importFrom vegan vegdist adonis2
#' @importFrom plyranges as_granges reduce_ranges group_by
#' @importFrom tidyr separate as_tibble
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @param sample_list list of epiMatrix coming from epiAnalysis function
#' @param metadata dataframe. Your samples metadata
#' @param colgroups character indicating the colname of Groups vector
#' @param colsamples character indicating the colname of Samples vector
#' @param rmUnmeth logical indicating if 0 methylated species should be removed of not from the analysis
#' @param cores num of cores
#' @param minGroups integer indicating the minimum number of unique Groups should be used for the dissimilarity analysis among samples
#' @param minSampleSize integer indicating the minimum number of samples per group needed to perform the statistical analysis
#' @param reduce logical indicating if intervals with the same statistics should be reduced or not
#' @return dataframe with rows indicating the genomic regions analysed and columns containing the statistics
#' @export


epiStat <- function(sample_list, metadata, colgroups, colsamples,
                    rmUnmeth = FALSE, cores = 1,
                    minGroups = 2, minSampleSize = 2, reduce = FALSE){
  ## Filtraggio regioni
  list = sample_list %>% map(select, id) %>%
    map(distinct, id) %>%
    purrr::map2(., names(.), ~ dplyr::mutate(.x, sample = .y)) %>%
    purrr::map_df(~ .) %>%
    dplyr::group_by(id) %>%
    dplyr::left_join(., metadata, by = c("sample" = colsamples)) %>%
    dplyr::group_by(id, .[colgroups]) %>% dplyr::count() %>%
    dplyr::filter(n >= minSampleSize) %>% dplyr::ungroup()
  ####
  list <- list %>% dplyr::group_by(id) %>% dplyr::mutate(ngroups = dplyr::n_distinct(get(colgroups))) %>% dplyr::filter(ngroups >= minGroups)

  ### Take filtered regions
  regs = unique(list$id)

  ### Split regions
  regs <- split(regs, (seq(length(regs))-1) %/% (length(regs)/cores))

  ### Split data
  filt_samples <- foreach::foreach(reg = regs) %do% {
    filt = purrr::map(sample_list, ~ dplyr::filter(., id %in% reg))
  }
  ### Do parallel
  cl <- parallel::makeCluster(cores, type = 'PSOCK')
  doParallel::registerDoParallel(cl)
  res = foreach::foreach(a = regs, b = filt_samples,
                         .combine = rbind,
                         .packages = c("magrittr", "dplyr", "purrr", "vegan"),
                         .export = c('oneStat', 'getEpimatrix')) %dopar% {
                          block = allStat(a, b, metadata, colgroups, colsamples, rmUnmeth, minGroups)
                         }
  parallel::stopCluster(cl)
  ## Filter results
  prova = res[stats::complete.cases(res)==TRUE,]
  ## Adjust p-values
  prova = prova %>% dplyr::mutate(p.adjust = p.adjust(p.value, method = "fdr"), .after = p.value)
  ## Separate columns
  prova = prova %>%
    tidyr::separate(id, c("seqnames", "start", "end"), "_",
                    remove = FALSE, convert = TRUE)

  if(reduce == TRUE){
    # Take groups names
    groups = colnames(prova[9:ncol(prova)])
    # Paste replicates number
    prova$replicates = do.call(paste, c(prova[9:ncol(prova)], sep = "_"))
    # Convert in granges
    gr = plyranges::as_granges(prova)
    # Group by all the statistcs and reduce ranges
    gr = gr %>% plyranges::group_by(F.statistics, p.value, num_EpiSpecies, replicates) %>%
      plyranges::reduce_ranges()
    # Convert again into a tibble
    prova = tidyr::as_tibble(gr)
    # Create the id column
    prova= prova %>% separate(replicates, groups, sep = "_") %>%
      mutate(id = paste(seqnames, start, end, sep = "_"), .before = seqnames) %>%
      mutate(p.adjust = p.adjust(p.value, method = "fdr"), .after = p.value)
  }
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

oneStat <- function(sample_list, region, metadata, colgroups, colsamples, rmUnmeth, minGroups){
  data <- getEpimatrix(sample_list, region)
  if(rmUnmeth == TRUE){
    data <- data[grep("1", colnames(data))]
    data <- data %>% dplyr::filter(!rowSums(.) == 0)
  }
    groups = unique(metadata[[colgroups]])
    metadata <- metadata %>%
      dplyr::filter(metadata[[colsamples]] %in% rownames(data))
    dist = vegan::vegdist(data, method="bray")
    if(length(unique(metadata[[colgroups]])) >= minGroups){
      a = groups %>% purrr::map_dbl(~ length(metadata[[colsamples]][metadata[[colgroups]] == .x]))
      a = t(as.data.frame(a, row.names = groups))
      rownames(a) = NULL
      p <- suppressMessages(vegan::adonis2(as.formula(paste("dist", colgroups, sep = "~")), data = metadata))
      result = data.frame("id" = region,
                         "F.statistics" = p$F[1],
                         "p.value" = p$`Pr(>F)`[1],
                         "num_EpiSpecies" = length(colnames(data)))
      result = cbind(result, a)
    } else {
      b = rep(NA, length(groups))
      b = t(as.data.frame(b, row.names = groups))
      rownames(b) = NULL
      result = data.frame("id" = region,
                          "F.statistics" = NA,
                          "p.value" = NA,
                          "num_EpiSpecies" = NA)
      result = cbind(result,b)
    }
  return(result)
}


allStat <- function(regions, sample_list, metadata, colgroups, colsamples, rmUnmeth, minGroups){
  out = regions %>% map_dfr(~ oneStat(sample_list, .x, metadata, colgroups, colsamples, rmUnmeth, minGroups))
  return(out)
}




