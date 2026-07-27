#' regsToPathways
#'
#' @importFrom seq2pathway runseq2pathway
#' @param stat_out list of samples intervals. It corresponds to the output 'epi' from the epiAnalysis function
#' @param genome your samples metadata. The input file must contain the columns named "Group" and "Samples"
#' @return list of Dataframes
#' @export

regsToPathway <- function(stat_out, genome){
  out = seq2pathway::runseq2pathway(as.data.frame(stat_out), genome = genome)
  return(out)
}




