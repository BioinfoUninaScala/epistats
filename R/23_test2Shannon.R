#' test2Shannon
#'
#' @importFrom purrr map reduce
#' @importFrom dplyr filter mutate select full_join group_by left_join ungroup nest_by summarise
#' @importFrom tidyr pivot_longer separate unnest
#' @param intervals_list list of samples intervals. It corresponds to the output 'epi' from the epiAnalysis function
#' @param metadata your samples metadata. The input file must contain the columns named "Group" and "Samples"
#' @return Dataframe with statistics
#' @export

test2Shannon <- function(intervals_list, metadata){
  check <- intervals_list %>% purrr::map(dplyr::filter, num_reads >= 50) %>%
                              purrr::map(dplyr::mutate, shanNorm = shannon / num_cg ,
                                                 id = paste(seqnames, start, end, sep = "_")) %>%
                              purrr::map(dplyr::filter, !duplicated(id)) %>%
                              purrr::map(dplyr::select, id, shanNorm) %>%
                              purrr::reduce(dplyr::full_join, by= "id") %>%
                              dplyr::filter(!rowSums(is.na(.)) > 0)

  colnames(check)[2:ncol(check)] <- names(intervals_list)

  inter_test <- check %>% tidyr::pivot_longer(cols = 2:ncol(check)) %>%
                          dplyr::group_by(id) %>%
                          dplyr::left_join(., metadata, by = c("name" = "Samples")) %>%
                          dplyr::group_by(id, Group, Points) %>%
                          dplyr::summarise(shan = mean(value)) %>% dplyr::ungroup()

  tests <- inter_test %>% dplyr::nest_by(id) %>%
                          dplyr::mutate(Model = list(summary(lm(shan ~ Points, data = data)))) %>%
                          dplyr::mutate(grad = tidy(Model)[["estimate"]][2],
                                        p.value= tidy(Model)[["p.value"]][2]) %>%
                          tidyr::unnest(id) %>% dplyr::ungroup() %>%
                          dplyr::select(1,5) %>%
                          tidyr::separate(id, c("chr", "start", "end"), "_", remove = FALSE)

  tests$p.adjust <- p.adjust(tests$p.value, method = "fdr")

  return(tests)

}
