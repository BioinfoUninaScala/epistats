#' diffAdjust
#'
#' @importFrom purrr map map2 map_df
#' @importFrom dplyr filter mutate select full_join group_by left_join ungroup nest_by summarise count as_tibble
#' @importFrom tidyr pivot_longer separate unnest
#' @importFrom broom tidy
#' @importFrom plyranges as_granges group_by reduce_ranges
#' @param intervals_list list of samples intervals. It corresponds to the output 'epi' from the epiAnalysis function
#' @param metadata your samples metadata. The input file must contain the columns named "Group" and "Samples"
#' @param statistic
#' @param groupcol
#' @param covariate
#' @param cores
#' @param reduce
#' @return Dataframe with statistics
#' @export

diffAdjust <- function(data,
                       metadata,
                       samples,
                       statistic,
                       groupcol,
                       covariate,
                       cores=1,
                       reduce=FALSE){
  # 1. Filter each dataframe of the samples list to keep only regions covered by at least 50 reads
  # 2. Add to each dataframe of the samples list two columns (shanNorm : normalized shannon, id : paste chr_start_end)
  # 3. Remove from each dataframe of the samples list the duplicated rows (you can have two epialleles with the highest frequency)
  # 4. Select from each dataframe only the shanNorm and the id columns
  # 5. Join by the id all the shanNorm columns from all the samples dataframes
  # 6. Take the rows that have no NAs (you have the ShanNorm value in all samples of your list)
  check <- data %>% purrr::map(mutate, id = paste(seqnames, start, end, sep = "_")) %>%
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


onef <- function(df, metadata, samples, statistic, groupcol, covariate, reduce){
  df <- df %>% dplyr::group_by(id) %>%
    dplyr::left_join(., metadata, by = c("sample" = samples)) %>%
    dplyr::ungroup()
  ### Check regions to test
  filt_step <- df %>% dplyr::group_by(id, .[groupcol]) %>% dplyr::count(.[covariate]) %>%
    dplyr::ungroup() %>% dplyr::filter(n >= 2) %>%
    dplyr::group_by(id, .[groupcol]) %>% dplyr::count() %>%
    dplyr::ungroup() %>% dplyr::filter(n >= 2) %>%
    dplyr::group_by(id) %>% dplyr::count() %>%
    dplyr::ungroup() %>% dplyr::filter(n >= 2)
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
                    p.value_Interaction = broom::tidy(Model)[["p.value"]][3]) %>%
      tidyr::unnest(id) %>% dplyr::ungroup() %>%
      dplyr::select(1,4:6) %>%
      tidyr::separate(id, c("seqnames", "start", "end"), "_", remove = FALSE, convert = TRUE)
    ## Reduce
    if(reduce == TRUE){
      res <- res %>% plyranges::as_granges() %>%
        plyranges::group_by(p.value_Group, p.value_Covariate, p.value_Interaction) %>%
        plyranges::reduce_ranges() %>%
        dplyr::as_tibble() %>%
        dplyr::mutate(id = paste(seqnames, start, end, sep = "_"), .before = seqnames)
    }
  }
  return(res)
}
