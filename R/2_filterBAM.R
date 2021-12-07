#' Filtering function
#'
#' @description
#' This function is provided to filter only selected regions from bamfiles (specific chromosomes of interest, etc...).
#' @importFrom GenomeInfoDb keepSeqlevels keepStandardChromosomes dropSeqlevels
#' @param algn GAlignments object containing aligned reads
#' @param keepStChr Logical indicating whether only standard chromosomes should be kept
#' @param keepM Logical indicating whether the mitochondrial chromosome should be kept for the analysis or not
#' @param retainChr Vector containing the chromosomes to be analysed as characters. Note that the strings should correspond to the seqlevels contained in the GAlignment object.
#' @return A filtered GAlignments object.
#' @export
#' @examples
#' data <- loadInput(bamfile, genomefile)
#' algn <- data[[1]]
#'
#' ## Keeping only standard chromosomes
#' filtered <- filterBAM(algn,
#'                       keepStChr =TRUE,
#'                       keepM = FALSE,
#'                       retainChr = NULL)
#'
#' ## Keeping only selected chromosomes
#' chr <- c("chr1", "chr2, "chr3")
#' filtered <- filterBAM(algn,
#'                       keepStChr = FALSE,
#'                       keepM = FALSE,
#'                       retainChr = chr)

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
