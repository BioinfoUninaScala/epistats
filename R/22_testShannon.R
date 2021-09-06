#' testShannon
#'
#' @importFrom purrr map reduce
#' @importFrom dplyr filter mutate select full_join group_by left_join ungroup nest_by n_distinct
#' @importFrom tidyr pivot_longer separate unnest
#' @param intervals_list list of samples intervals. It corresponds to the output 'epi' from the epiAnalysis function
#' @param metadata your samples metadata. The input file must contain the columns named "Group" and "Samples"
#' @return Dataframe with statistics and the test used to compare the Shannon between samples
#' @export

testShannon <- function(intervals_list, metadata){
  regs <- intervals_list %>% purrr::map(dplyr::filter, num_reads >= 50) %>%
                             purrr::map(dplyr::mutate, shanNorm = shannon / num_cg ,
                                         id = paste(seqnames, start, end, sep = "_")) %>%
                             purrr::map(dplyr::filter, !duplicated(id)) %>%
                             purrr::map(dplyr::select, id, shanNorm) %>%
                             purrr::reduce(dplyr::full_join, by= "id")

  colnames(regs)[2:ncol(regs)] <- names(intervals_list)

  filtered <- regs %>% tidyr::pivot_longer(cols = 2:ncol(.)) %>%
                        dplyr::filter(!rowSums(is.na(.)) > 0) %>%
                        dplyr::group_by(id) %>%
                        dplyr::left_join(., metadata, by = c("name" = "Samples")) %>%
                        dplyr::group_by(id, Group) %>%
                        dplyr::mutate(samples = dplyr::n_distinct(name)) %>%
                        dplyr::filter(samples >= 2) %>%
                        dplyr::group_by(id) %>%
                        dplyr::mutate(groups = dplyr::n_distinct(Group)) %>%
                        dplyr::filter(groups >= 2 & samples >= 2) %>%
                        dplyr::ungroup() %>%
                        dplyr::select(1:4)

  res <- filtered %>% dplyr::nest_by(id) %>%
                      dplyr::mutate(test = ifelse(dplyr::n_distinct(data$Group) == 2,
                                                  wilcox.test(value ~ Group, data = data)$p.value,
                                                  kruskal.test(value ~ Group, data = data)$p.value),
                                    type = ifelse(dplyr::n_distinct(data$Group) == 2,
                                                  "wilcox.test",
                                                  "kruskal.test")) %>%
                      tidyr::unnest(id) %>%
                      dplyr::ungroup() %>%
                      dplyr::select(1,3,4) %>%
                      tidyr::separate(id, c("chr", "start", "end"), "_", remove = FALSE) %>%
                      dplyr::filter(is.numeric(test) == TRUE)

  res$p.adjust <- p.adjust(res$test, method = "fdr")

  return(res)
}
