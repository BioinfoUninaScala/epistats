#' filterBAM
#'
#' @importFrom GenomeInfoDb keepSeqlevels keepStandardChromosomes dropSeqlevels
#' @param algn GAlignments object
#' @param keepStChr Logical indicating if only standard chromosomes should be kept
#' @param keepM Logical indicating if the mitochondrial chromosome should be kept or not
#' @param retainChr Vector containing the strings of the chromosomes to be analysed
#' @return filtered GAlignments object
#' @export

filterBAM <- function(algn,
                      keepStChr=TRUE,
                      keepM=FALSE,
                      retainChr=NULL)
{
  cat("filtering reads in bam\n" )
  if(!is.null(retainChr)==TRUE){
    cat("filtering reads that align to user provided intervals\n")
    algn <- GenomeInfoDb::keepSeqlevels(algn, retainChr, pruning.mode = "coarse")
    cat("done\n")
  } else {
    if(keepStChr==TRUE){
      cat("filtering reads that align to standard chromosomes\n")
      algn <- GenomeInfoDb::keepStandardChromosomes(algn, pruning.mode = "coarse")
      cat("done\n")
    }
    if(keepM==FALSE){
      cat("filtering reads that align to chr M \n")
      algn <- GenomeInfoDb::dropSeqlevels(algn, c("chrM","MT"), pruning.mode = "coarse")
      cat("done")
    }
  }
  return(algn)
}
