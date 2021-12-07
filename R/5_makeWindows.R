#' Designing sliding windows with a fixed length size
#'
#' @description
#' One of the two available methods to build target regions. This functions works by designing sliding windows with user-defined values of length and step size. Target regions will contain different numbers of CpG sites.
#' @importFrom GenomicRanges slidingWindows start GRanges countOverlaps strand
#' @importFrom parallel mclapply
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom Biostrings matchPattern reverseComplement DNAString
#' @param gr A GRanges object containing the filtered regions obtained to design target regions.
#' @param window An integer indicating the length size of the sliding windows
#' @param step An integer indicating the step size used to slide windows
#' @param Genome A DNAStringSet object containing the fasta Genome (obtained from the loadInput() function).
#' @param mode A character indicating the pattern to be used to analyse DNA methylation (one of "CG", "CC", "CA", "CT").
#' @param min.C Numeric indicating the minimum number of CpGs to be displayed in the sliding window (Default = 2).
#' @param max.C Numeric indicating the maximum number of CpGs to be displayed in the sliding window (Default = 50)
#' @param cores Numeric indicating the number of cores to be used for the computation
#' @return A GRanges object containing the newly designed sliding windows with user-defined features.
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
#' ## Selecting regions with a minimum coverage
#' covered <- filterByCoverage(algn, threshold = 50, minsize = 8)
#'
#' ## Designing customised sliding windows
#' windows <- makeWindows(gr = covered,
#'                        Genome = Genome,
#'                        window = 50,
#'                        step = 1,
#'                        mode = "CG",
#'                        min.C = 2,
#'                        max.C = 50,
#'                        cores = 40)

makeWindows <- function(gr,
                        Genome,
                        window = 50,
                        step = 1,
                        mode = "CG",
                        min.C = 1,
                        max.C = 50,
                        cores = 1){
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
                         mode,
                         min.C,
                         max.C,
                         cores)
{
  ###find cmode coordinates in the Genome
  c_pos <- findPattern(Genome,cores=cores,mode=mode)
  #####filters windows with a number of cmode in the user supplied range
  windows_filt <- windows[GenomicRanges::countOverlaps(windows, c_pos, minoverlap = 2)>=min.C]
  windows_filt <- windows_filt[GenomicRanges::countOverlaps(windows_filt, c_pos, minoverlap = 2)<=max.C]
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



