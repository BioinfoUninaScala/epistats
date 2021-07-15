#' makeWindows
#'
#' @importFrom GenomicRanges slidingWindows start GRanges countOverlaps strand
#' @importFrom parallel mclapply
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom Biostrings matchPattern reverseComplement DNAString
#' @param gr GRanges
#' @param window length size of the sliding window
#' @param step step size to slide
#' @param Genome fasta Genome
#' @param mode pattern character
#' @param min.C min num of CpG in the sliding window
#' @param max.C max num of CpG in the sliding window
#' @param cores num of cores to use
#' @return GRanges
#' @export

makeWindows <- function(gr, window, step, Genome, mode, min.C, max.C, cores){
  ##split gr in sliding windows
  windows <- GenomicRanges::slidingWindows(gr, window, step)
  if (mode=="CG")
  {
    data=filterWindows(windows@unlistData, Genome, mode, min.C, max.C, cores = cores)
    windows_filt=data[["windows_filt"]]
    c_pos=data[["c_pos"]]
  }else{
    data_plus <- filterWindows(windows@unlistData, Genome, mode, min.C, max.C, cores = cores)
    data_plus=lapply(data_plus, function(x) GenomicRanges::strand(x)="+")
    data_minus=filterWindows(windows@unlistData, Genome, Biostrings::reverseComplement(Biostrings::DNAString(mode)), min.C, max.C, cores = cores)
    data_minus=lapply(data_minus, function(x) GenomicRanges::strand(x)="-")
    ####unisce plus e minus
    c_pos=c(data_plus[["c_pos"]], data_minus[["cpos"]])
    windows_filt=c(data_plus[["windows_filt"]], data_minus[["windows_filt"]])
  }
  return(list("windows"=windows_filt, "c_pos"=c_pos))
}


filterWindows = function(windows,
                         Genome,
                         mode,min.C,
                         max.C,
                         cores)
{
  ###find cmode coordinates in the Genome
  c_pos <- findPattern(Genome,cores=cores,mode=mode)
  #####filters windows with a number of cmode in the user supplied range
  windows_filt <- windows[GenomicRanges::countOverlaps(windows, c_pos)>=min.C]
  windows_filt <- windows_filt[GenomicRanges::countOverlaps(windows_filt, c_pos)<=max.C]
  #####select Cmode ranges that fall within the selected windows
  c_pos <- c_pos[GenomicRanges::countOverlaps(c_pos, windows_filt)>0]
  return(list("windows_filt"=windows_filt, "c_pos"=c_pos))
}


findPattern <- function(Genome, cores, mode="CG"){
  Pattern <- parallel::mclapply(GenomeInfoDb::seqlevels(Genome), function(x) GenomicRanges::start(Biostrings::matchPattern(mode, Genome[[x]])), mc.cores = cores)
  return(
    suppressWarnings(
      base::do.call(c, parallel::mclapply(1:length(GenomeInfoDb::seqlevels(Genome)), function(x) GenomicRanges::GRanges(base::names(Genome)[x],
                                                                                                                        IRanges::IRanges(Pattern[[x]], width = 2)
      ), mc.cores=cores))
    )
  )
}
