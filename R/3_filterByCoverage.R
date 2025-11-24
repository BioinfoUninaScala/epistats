#' Selecting regions with a minimum user-defined coverage
#'
#' @description
#' The functions allows the selection of genomic regions that satisfy a minimum coverage threshold defined by the user as a parameter.
#' @importFrom GenomicAlignments coverage
#' @importFrom GenomicRanges GRanges width
#' @param algn GAlignments object to filter.
#' @param threshold Integer indicating the minimum coverage required to keep genomic regions
#' @param minsize Integer indicating the minimum size of the genomic range. This value is also depending from the method of choice that will be used by the dedicated functions makeBins() and makeWindows(). If you use makeBins() the minimum size will be equal to the length of n(CpG sites to be included in the target region)x2. Instead, if you use makeWindows, the minimum size of the genomic regions used to build target regions will be equal to the size of the sliding window that will be used next.
#' @return filtered GRanges object.
#' @export


filterByCoverage <- function(algn,
                             threshold = 50,
                             minsize = 50)
{
  cat("filtering regions with minimum coverage", threshold, "\n")
  cvg <- GenomicAlignments::coverage(algn)                                      ### coverage per chr
  cvg <- GenomicRanges::GRanges(cvg)                                            ### estraggo i GRanges di queste regioni
  cvg <- cvg[cvg$score >= threshold]                                            ### seleziono i GRanges che hanno un n di reads(threshold) in particolare
  cvg <- cvg[cvg$score <= stats::quantile(cvg$score, probs = seq(0, 1, 0.001))[1000]]
  cvg <- IRanges::reduce(cvg)                                                   ### unisco i GRanges contigui
  cvg = cvg[which(GenomicRanges::width(cvg)>= minsize)]                         ### elimino i ranges con width<1
  cat("done\n")
  return(cvg)
}
