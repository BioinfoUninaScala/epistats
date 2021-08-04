#' filterByCoverage
#'
#' @importFrom GenomicAlignments coverage
#' @importFrom GenomicRanges GRanges width
#' @importFrom IRanges reduce
#' @param algn GAlignments object
#' @param threshold Integer indicating the minimum depth value to take genomic regions
#' @param minsize Interger indicating the minumum size of the genomic range
#' @return filtered GRanges object
#' @export


filterByCoverage <- function(algn,
                             threshold = 50,
                             minsize = 50)
{
  cat("filtering regions with minimum coverage", threshold, "\n")
  cvg <- GenomicAlignments::coverage(algn)                                      ### coverage per chr
  cvg <- GenomicRanges::GRanges(cvg)                                            ### estraggo i GRanges di queste regioni
  cvg <- cvg[cvg$score >= threshold]                                            ### seleziono i GRanges che hanno un n di reads(threshold) in particolare
  cvg <- cvg[cvg$score <= quantile(cvg$score, probs = seq(0, 1, 0.001))[1000]]
  cvg <- IRanges::reduce(cvg)                                                   ### unisco i GRanges contigui
  cvg = cvg[which(GenomicRanges::width(cvg)>= minsize)]                         ### elimino i ranges con width<1
  cat("done\n")
  return(cvg)
}
