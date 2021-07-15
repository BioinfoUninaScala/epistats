#' filterBAM
#'
#' @importFrom GenomeInfoDb keepSeqlevels keepStandardChromosomes dropSeqlevels
#' @param algn GAlignments object
#' @param keepStChr logical
#' @param keepM logical
#' @param retainChr logical
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
