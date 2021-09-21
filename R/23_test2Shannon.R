#' test2Shannon
#'
#' @importFrom purrr map reduce
#' @importFrom dplyr filter mutate select full_join group_by left_join ungroup nest_by summarise
#' @importFrom tidyr pivot_longer separate unnest
#' @importFrom broom tidy
#' @param intervals_list list of samples intervals. It corresponds to the output 'epi' from the epiAnalysis function
#' @param metadata your samples metadata. The input file must contain the columns named "Group" and "Samples"
#' @return Dataframe with statistics
#' @export

test2Shannon <- function(intervals_list, metadata){
  # 1. Filter each dataframe of the samples list to keep only regions covered by at least 50 reads
  # 2. Add to each dataframe of the samples list two columns (shanNorm : normalized shannon, id : paste chr_start_end)
  # 3. Remove from each dataframe of the samples list the duplicated rows (you can have two epialleles with the highest frequency)
  # 4. Select from each dataframe only the shanNorm and the id columns
  # 5. Join by the id all the shanNorm columns from all the samples dataframes
  # 6. Take the rows that have no NAs (you have the ShanNorm value in all samples of your list)
  check <- intervals_list %>% purrr::map(dplyr::filter, num_reads >= 50) %>%
    purrr::map(dplyr::mutate, shanNorm = shannon / num_cg ,
               id = paste(seqnames, start, end, sep = "_")) %>%
    purrr::map(dplyr::filter, !duplicated(id)) %>%
    purrr::map(dplyr::select, id, shanNorm) %>%
    purrr::reduce(dplyr::full_join, by= "id")
  # Set the sample names as colnames
  colnames(check)[2:ncol(check)] <- names(intervals_list)
  # 1. It makes the dataframes longer
  # 2. It removes the rows with NAs (the samples with no shanNorm value)
  # 3. It groups the dataframe by id region
  # 4. It joins the Points value from metadata to the shan values
  totest <- check %>% tidyr::pivot_longer(cols = 2:ncol(check)) %>%
    dplyr::filter(!rowSums(is.na(.)) > 0) %>%
    dplyr::group_by(id) %>%
    dplyr::left_join(., metadata, by = c("name" = "Samples")) %>% ungroup()
  # Filter step to select only regions that have shan values for at least 2 points, each with 2 samples
  filt_step <- totest %>% group_by(id, Group) %>% count(Points) %>%
    ungroup() %>% filter(n >= 2) %>%
    group_by(id, Group) %>% count() %>%
    ungroup %>% filter(n >= 2) %>%
    group_by(id) %>% count() %>%
    ungroup() %>% filter(n >= 2)
  # It selects from totest only the regions that passed the filter
  totest <- totest %>% ungroup() %>% filter(id %in% filt_step$id)
  ##
  tests <- totest %>% dplyr::nest_by(id) %>%
    dplyr::mutate(Model = list(aov(value ~ Points*Group, data = data))) %>%
    dplyr::mutate(p.value_Points = broom::tidy(Model)[["p.value"]][1],
                  p.value_Group = broom::tidy(Model)[["p.value"]][2],
                  p.value_PointsGroup = broom::tidy(Model)[["p.value"]][3]) %>%
    tidyr::unnest(id) %>% dplyr::ungroup() %>%
    dplyr::select(1,4:6) %>%
    tidyr::separate(id, c("chr", "start", "end"), "_", remove = FALSE) %>%
    dplyr::mutate(p.adjust_Points = p.adjust(p.value_Points, method = "fdr"), .after = p.value_Points) %>%
    dplyr::mutate(p.adjust_Group = p.adjust(p.value_Group, method = "fdr"), .after = p.value_Group) %>%
    dplyr::mutate(p.adjust_PointsGroup = p.adjust(p.value_PointsGroup, method = "fdr"), .after = p.value_PointsGroup)
  ###
  return(tests)

}
