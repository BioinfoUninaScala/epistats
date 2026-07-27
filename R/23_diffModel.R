#' Getting regions which differ for one statistic (such as Shannon Entropy, etc...) through linear model and possibly adjusting for a covariate
#' @description
#' This function is implemented to get different significant regions using a linear model adjusted for a covariate (different Points in time, different tissue, etc...).
#' @importFrom purrr map map2 map_df
#' @importFrom dplyr filter mutate select full_join group_by left_join ungroup nest_by summarise count as_tibble
#' @importFrom tidyr pivot_longer separate unnest
#' @importFrom broom tidy
#' @importFrom plyranges as_granges group_by reduce_ranges
#' @param intervals_list A list object. It corresponds to a list containing the summary outputs obtained through epiAnalysis() function (the ones containing the summary statistics, such as Shannon Entropy, Mean CpGs distance, etc...).
#' @param metadata Samples metadata provided as table. This input data should contain specific columns for samples names and groups distinction.
#' @param statistic Character indicating the column name that contains the statistic vector the user wants to use to perform the test (e.g., "Shannon", "mean_met").
#' @param groupcol Character indicating the column name of the Groups in the metadata used to compare the statistics.
#' @param covariate Character indicating the column name of the covariate to be used to adjust the linear model.
#' @param min.per.group An integer indicating the minimum number of samples for group required to perform the test.
#' @param min.groups An integer indicating the minimum number of groups wanted to perform the test.
#' @param cores An integer indicating the number of cores to be used to perform the computation.
#' @param reduce Logical indicating whether adjacent overlapping intervals should be reduced as a unique interval or not.
#' @return A dataframe containing the statistical test results.
#' @export
#' @examples
#' data(epistats)
#' samples_list <- list(Sample1_intervals.bed,
#'                      Sample2_intervals.bed,
#'                      Sample3_intervals.bed,
#'                      Sample4_intervals.bed)
#'
#' diffresult <- diffAdjust(intervals_list = samples_list,
#'                          metadata = ann,
#'                          samples = "Samples",
#'                          statistic = "Shannon",
#'                          groupcol = "Group",
#'                          covariate = "Time",
#'                          cores = 40,
#'                          reduce = FALSE)

diffModel <- function(intervals_list,
                      metadata,
                      samples,
                      statistic,
                      groupcol,
                      covariate,
                      min.per.group=2,
                      min.groups=2,
                      cores=1,
                      reduce=FALSE){
  # 1. Filter each dataframe of the samples list to keep only regions covered by at least 50 reads
  # 2. Add to each dataframe of the samples list two columns (shanNorm : normalized shannon, id : paste chr_start_end)
  # 3. Remove from each dataframe of the samples list the duplicated rows (you can have two epialleles with the highest frequency)
  # 4. Select from each dataframe only the shanNorm and the id columns
  # 5. Join by the id all the shanNorm columns from all the samples dataframes
  # 6. Take the rows that have no NAs (you have the ShanNorm value in all samples of your list)
  check <- intervals_list %>% purrr::map(mutate, id = paste(seqnames, start, end, sep = "_")) %>%
                              purrr::map(dplyr::filter, !duplicated(id)) %>%
                              purrr::map(dplyr::select, id, all_of(statistic)) %>%
                              purrr::map2(., names(.), ~ dplyr::mutate(.x, sample = .y)) %>%
                              purrr::map_df(~ .)
  ## Find all regions
  regs <- unique(check$id)
  ## Split regs
  regs <- split(regs, (seq(length(regs))-1) %/% (length(regs)/cores))
  ## Split dataframe
  dfs <- purrr::map(regs, ~ filter(check, id %in% .))
  ## Do parallel
  cl <- parallel::makeCluster(cores, type = 'PSOCK')
  doParallel::registerDoParallel(cl)
  res = foreach::foreach(a = dfs, .combine = rbind, .packages = c("magrittr", "dplyr")) %dopar% {
    out = onef(df = a,
               metadata = metadata,
               samples = samples,
               statistic = statistic,
               groupcol = groupcol,
               covariate = covariate,
               min.per.group = min.per.group,
               min.groups = min.groups,
               reduce = reduce)
  }
  parallel::stopCluster(cl)
  ## Adjust p.values
  res = res %>%
    dplyr::mutate(p.adjust_Group = p.adjust(p.value_Group, method = "fdr"), .after = p.value_Group) %>%
    dplyr::mutate(p.adjust_Covariate = p.adjust(p.value_Covariate, method = "fdr"), .after = p.value_Covariate) %>%
    dplyr::mutate(p.adjust_Interaction = p.adjust(p.value_Interaction, method = "fdr"), .after = p.value_Interaction)
  return(res)
}


onef <- function(df,
                 metadata,
                 samples,
                 statistic,
                 groupcol,
                 covariate,
                 min.per.group,
                 min.groups,
                 reduce){
  df <- df %>% dplyr::group_by(id) %>%
    dplyr::left_join(., metadata, by = c("sample" = samples)) %>%
    dplyr::ungroup()
  ### Check regions to test
  filt_step <- df %>% dplyr::group_by(id, .[groupcol]) %>% dplyr::count(.[covariate]) %>%
    dplyr::ungroup() %>% dplyr::filter(n >= min.per.group) %>%
    dplyr::group_by(id, .[groupcol]) %>% dplyr::count() %>%
    dplyr::ungroup() %>% dplyr::filter(n >= 2) %>%
    dplyr::group_by(id) %>% dplyr::count() %>%
    dplyr::ungroup() %>% dplyr::filter(n >= min.groups)
  ### Subset df regions that passed the filter step
  totest <- df %>% dplyr::filter(id %in% filt_step$id)
  ### Check if totest is an empty df
  if(dim(totest)[1] == 0){
    res = totest
  } else {
    ## Perform analysis
    res <- totest %>% dplyr::nest_by(id) %>%
      dplyr::mutate(Model = list(aov(as.formula(paste(statistic, paste(groupcol, covariate, sep = "*"), sep = "~")), data = data))) %>%
      dplyr::mutate(p.value_Group = broom::tidy(Model)[["p.value"]][1],
                    p.value_Covariate = broom::tidy(Model)[["p.value"]][2],
                    p.value_Interaction = broom::tidy(Model)[["p.value"]][3],
                    intercept = coefficients(Model)[[1]],
                    group = coefficients(Model)[[2]],
                    covariate = coefficients(Model)[[3]],
                    interaction = coefficients(Model)[[4]]) %>%
      tidyr::unnest(id) %>% dplyr::ungroup() %>%
      dplyr::select(1,4:10) %>%
      tidyr::separate(id, c("seqnames", "start", "end"), "_", remove = FALSE, convert = TRUE)
    ## Reduce
    if(reduce == TRUE){
      varcols <- colnames(res)[-c(1:4)]
      res <- res %>% plyranges::as_granges() %>%
        plyranges::group_by(!!!rlang::syms(varcols)) %>%
        plyranges::reduce_ranges() %>%
        dplyr::as_tibble() %>%
        dplyr::mutate(id = paste(seqnames, start, end, sep = "_"), .before = seqnames)
    }
  }
  return(res)
}
