#' diffStat
#'
#' @importFrom purrr map reduce
#' @importFrom dplyr filter mutate select full_join group_by left_join ungroup nest_by n_distinct
#' @importFrom tidyr pivot_longer separate unnest
#' @param intervals_list list of samples intervals. It corresponds to the output 'epi' from the epiAnalysis function
#' @param metadata your samples metadata. The input file must contain the columns named "Group" and "Samples"
#' @param statistic parameter indicating the colname of the statistic the user wants to use to perform the test
#' @param groupcol character indicating the Group colname in the metadata used to compare the statistics
#' @param cores number of cores to use
#' @param reduce logical indicating whether adiacent overlapping intervals should be reduced or not
#' @return Table containing the p-values obtained comparing the given statistic between samples of different groups
#' @export

diffStat <- function(intervals_list,
                     metadata,
                     colsamples,
                     statistic,
                     groupcol,
                     cores=1,
                     reduce=FALSE){
  ## Create id for each dataframe and select only id and the statistic to compare
  lstdata <- intervals_list %>% purrr::map(dplyr::mutate, id = paste(seqnames, start, end, sep = "_")) %>%
                             purrr::map(dplyr::select, id, all_of(statistic))
  lstdata <- Map(function(x,y) mutate(x, sample = y), lstdata, names(lstdata))
  ## Take all the observed regions (also not in common)
  regs <- 1:length(lstdata) %>% map(function(x) lstdata[[x]]$id)
  ## Create an unique vector
  regs <- unique(c(unlist(regs)))
  ## Divide regions in blocks by number of cores
  regs <- split(regs, (seq(length(regs))-1) %/% (length(regs)/cores))
  ## Filter data in blocks
  datablocks <- foreach(reg = regs) %do% {
    map(lstdata, ~ dplyr::filter(., id %in% reg))
  }
  ## Parallel function
  cl <- parallel::makeCluster(cores, type = 'PSOCK')
  doParallel::registerDoParallel(cl)
  res = foreach::foreach(a = datablocks, .combine = rbind, .packages = c("magrittr", "dplyr")) %dopar% {
    out = onefun(datalst = a,
                 metadata = metadata,
                 colsamples = colsamples,
                 statistic = statistic,
                 groupcol = groupcol,
                 reduce = reduce)
  }
  parallel::stopCluster(cl)
  ### It adjusts p.values
  res$p.adjust <- p.adjust(res$p.value, method = "fdr")
  ###
  return(res)
}

## Table with all data
onefun <- function(datalst, metadata, colsamples, statistic, groupcol, reduce){
  datalst = datalst %>% purrr::map_df(~ .)
  ## Filter regions
  filter <- datalst %>% dplyr::left_join(., metadata, by = c("sample" = colsamples)) %>%
    dplyr::group_by(id, .[groupcol]) %>%
    dplyr::mutate(samples = dplyr::n_distinct(sample)) %>%
    dplyr::filter(samples >= 2) %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(groups = dplyr::n_distinct(all_of(across(groupcol)))) %>%
    dplyr::filter(groups >= 2 & samples >= 2) %>%
    dplyr::ungroup() %>%
    dplyr::select(1:4)
  ## Find Diff
  if(dim(filter)[1] == 0){
    res = filter
  } else {
    ### It performs the test for each region (wilcoxon if there are only two groups, otherwise the kruskal one)
    res <- filter %>% dplyr::nest_by(id) %>%
      dplyr::mutate(p.value = ifelse(dplyr::n_distinct(data[groupcol]) == 2,
                                     wilcox.test(as.formula(paste(statistic, groupcol, sep = "~")), data = data)$p.value,
                                     kruskal.test(as.formula(paste(statistic, groupcol, sep = "~")), data = data)$p.value),
                    type = ifelse(dplyr::n_distinct(data[groupcol]) == 2,
                                  "wilcox.test",
                                  "kruskal.test")) %>%
      tidyr::unnest(id) %>%
      dplyr::ungroup() %>%
      dplyr::select(1,4,3) %>%
      tidyr::separate(id, c("seqnames", "start", "end"), "_", remove = FALSE, convert = TRUE) %>%
      dplyr::filter(!is.nan(p.value) == TRUE)
    if(reduce == TRUE){
      res <- res %>% as_granges %>%
        group_by(type, p.value) %>%
        reduce_ranges() %>%
        as_tibble() %>%
        mutate(id = paste(seqnames, start, end, sep = "_"), .before = seqnames)
    }
  }
  return(res)
}
