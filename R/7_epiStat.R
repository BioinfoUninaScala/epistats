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


epiStat <- function(sample_list, metadata, cores){
  list = purrr::map(sample_list, function(x) x %>%
                      dplyr::select(id) %>%
                      dplyr::distinct %>%
                      dplyr::mutate(present = TRUE))
  ## names
  filt = list %>%
    purrr::reduce(dplyr::full_join, by = "id") %>%
    dplyr::mutate(count_na = rowSums(is.na(.))) %>%
    dplyr::filter(count_na < (length(sample_list)-2))
  #### filter(!count_na == (length(sample_list)-1))
  regions = filt$id
  regs <- split(regions, (seq(length(regions))-1) %/% (length(regions)/cores))
  future::plan(multisession, workers = cores)
  prova = regs %>% furrr::future_map_dfr(~allStat(., sample_list, metadata))
  # prova = mclapply(1:length(regs), function(x) allStat(regs[[x]], sample_list, metadata),
  #                  mc.cores = cores)
  # prova = ldply(prova, data.frame)
  prova = prova[stats::complete.cases(prova)==TRUE,]
  prova$p.adjust = stats::p.adjust(prova$p.value, method = "fdr")
  return(prova)
}


mkmatr <- function(samples, region){
  filt = purrr::map(samples, ~ dplyr::filter(., id == region))
  filt = purrr::map(filt, ~ dplyr::select(., 1:2))
  counts = filt %>% purrr::reduce(dplyr::full_join, by = "Var1")
  rownames(counts)= counts$Var1
  counts$Var1=NULL
  colnames(counts) = names(samples)
  data = as.data.frame(t(counts))
  data[is.na(data)] = 0
  data = data %>% dplyr::filter(., !rowSums(.) == 0)
  return(data)
}


oneStat <- function(sample_list, region, metadata){
  data <- mkmatr(sample_list, region)
  #data = data[!rowSums(data)== 0,]
  dist = vegan::vegdist(data, method="bray")
  metadata = metadata %>%
    filter(., Samples %in% rownames(data))
  if(length(unique(metadata$Group)) > 1){
    ## problema formula
    p= vegan::adonis2(dist ~ Group, data = metadata)
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


allStat <- function(regions, sample_list, metadata){
  filt = purrr::map(sample_list, ~ dplyr::filter(., id %in% regions))
  out = regions %>% map_dfr(~ oneStat(filt, .x, metadata))
}
