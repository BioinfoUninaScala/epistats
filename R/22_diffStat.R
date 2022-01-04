#' Calculate regions that differ for a statistic
#'
#' Getting the regions that differ for one statistic (e.g. Shannon Entropy, Average DNA Methylation, ...) among distinct conditions.
#'
#' @importFrom purrr map reduce
#' @importFrom dplyr filter mutate select full_join group_by left_join ungroup nest_by n_distinct
#' @importFrom tidyr pivot_longer separate unnest
#' @param intervals_list A list object. It corresponds to a list containing the summary outputs obtained through epiAnalysis() function (the ones containing the summary statistics, such as Shannon Entropy, Mean CpGs distance, etc...).
#' @param metadata Samples metadata provided as table. This input data should contain specific columns for samples names and groups distinction.
#' @param statistic Character indicating the column name that contains the statistic vector the user wants to use to perform the test (e.g., "Shannon", "mean_met").
#' @param groupcol Character indicating the column name of the Groups in the metadata used to compare the statistics
#' @param min.per.group An integer indicating the minimum number of samples for group required to perform the test.
#' @param min.groups An integer indicating the minimum number of groups wanted to perform the test.
#' @param cores Integer indicating the number of cores used to perform the computation
#' @param reduce Logical indicating whether adjacent overlapping intervals should be reduced as a unique interval or not
#' @return A dataframe containing the statistical test results.
#'
#' id = Regions coordinates given as ID.
#'
#' p.value = significance of the difference found for each analysed region
#'
#' test = test used to perform the analysis (depending on the number of distinct groups used for the test)
#'
#' @export
#' @examples
#' data(epistats)
#' samples_list <- list(Sample1_intervals.bed,
#'                      Sample2_intervals.bed,
#'                      Sample3_intervals.bed,
#'                      Sample4_intervals.bed)
#'
#' diff <- diffStat(intervals_list = samples_list,
#'                  metadata = ann,
#'                  colsamples = "Samples",
#'                  statistic = "Shannon",
#'                  groupcol = "Group",
#'                  min.per.group = 2,
#'                  min.groups = 2,
#'                  cores = 40,
#'                  reduce = FALSE)

diffStat <- function(intervals_list,
                     metadata,
                     colsamples,
                     statistic,
                     groupcol,
                     min.per.group=2,
                     min.groups=2,
                     cores=1,
                     reduce=FALSE){
  ## Create id for each dataframe and select only id and the statistic to compare
  lstdata <- intervals_list %>% purrr::map(dplyr::mutate, id = paste(seqnames, start, end, sep = "_")) %>%
                                purrr::map(dplyr::filter, !duplicated(id)) %>%
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
  res = foreach::foreach(a = datablocks, .combine = rbind, .packages = c("magrittr", "dplyr"), .export = 'onefun') %dopar% {
    out = onefun(datalst = a,
                 metadata = metadata,
                 colsamples = colsamples,
                 statistic = statistic,
                 groupcol = groupcol,
                 min.per.group = min.per.group,
                 min.groups = min.groups,
                 reduce = reduce)
  }
  parallel::stopCluster(cl)
  ### It adjusts p.values
  res$p.adjust <- p.adjust(res$p.value, method = "fdr")
  ###
  return(res)
}

## Table with all data
onefun <- function(datalst, metadata, colsamples, statistic, groupcol, min.per.group, min.groups, reduce){
  datalst = datalst %>% purrr::map_df(~ .)
  ## Filter regions
  filter <- datalst %>% dplyr::left_join(., metadata, by = c("sample" = colsamples)) %>%
    dplyr::group_by(id, .[groupcol]) %>%
    dplyr::mutate(samples = dplyr::n_distinct(sample)) %>%
    dplyr::filter(samples >= min.per.group) %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(groups = dplyr::n_distinct(all_of(across(groupcol)))) %>%
    dplyr::filter(groups >= min.groups & samples >= min.per.group) %>%
    dplyr::ungroup() %>%
    dplyr::select(1:3, groupcol)
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
