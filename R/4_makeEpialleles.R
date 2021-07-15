#' makeEpialleles
#'
#' @importFrom parallel mclapply
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom IRanges IRanges IRangesList
#' @importFrom Biostrings matchPattern
#' @importFrom GenomicRanges findOverlaps countOverlaps
#' @importFrom S4Vectors countQueryHits queryHits
#' @param gr GRanges
#' @param Genome fasta file
#' @param mode character pattern
#' @param n num CpG in the bin
#' @param min.binsize min length of a bin
#' @param max.binsize max length of a bin
#' @param cores number of cores to use
#' @return GRanges
#' @export

makeEpialleles <- function(gr,
                           Genome,
                           mode="CG",
                           n=4,
                           min.binsize=1,
                           max.binsize=200,
                           cores=10)
{
  if(mode=="CG"){
    data= findEpialleles(gr, Genome, mode, n, min.binsize, max.binsize, cores)
    bins= data[[1]]
    c_pos= data[[2]]
  } else {
    data_plus= findEpialleles(gr, Genome, mode, n, min.binsize, max.binsize, cores)
    data_minus= findEpialleles(gr, Genome, mode= Biostrings::reverseComplement(Biostrings::DNAString(mode)), n, min.binsize, max.binsize, cores)
    bins= c(data_plus[[1]], data_minus[[1]])
    c_pos= c(data_plus[[2]], data_minus[[2]])
  }
  return(list("bins"=bins,"c_pos"=c_pos))
}


findEpialleles <- function(gr,
                           Genome,
                           mode,
                           n,
                           min.binsize,
                           max.binsize,
                           cores)
{
  ###find Cmode in Genome
  Genome=Genome[names(Genome) %in% seqlevels(gr)]
  cat("retrieving", mode, "coordinates from genomic sequence")
  c_pos= findPattern(Genome, cores=cores, mode)
  ###find overlaps between Cmode and intervals
  overlaps <- GenomicRanges::findOverlaps(gr, c_pos)
  #####select intervals with number of Cmode >=n
  gr_filt = gr[S4Vectors::countQueryHits(overlaps) >= n]
  ###select only Cmode coord that overlap with the selected intervals
  c_filt = c_pos[GenomicRanges::countOverlaps(c_pos, gr_filt)>0]
  ####Find overlaps between filtered objects
  overlaps_filt <- GenomicRanges::findOverlaps(gr_filt, c_filt)
  ####add interval id to cpos. In this way cpg that overlap the same interval will share the same id
  c_filt$id <- S4Vectors::queryHits(overlaps_filt)
  ####split cpos in list by interval id
  byint <- split(c_filt, f=c_filt$id)
  ####apply find bins function to each element of the list
  cat("finding epialleles in Genome \n")
  bins <- do.call(c,lapply(byint, function(x) findBins(x,n)))
  cat("done\n")
  ####unlist  bins
  bins=unlist(IRanges::IRangesList(bins),use.names = FALSE)
  ####filter bins by size
  cat("filtering epialleles with max bin size", max.binsize, "and min bin size", min.binsize, "\n")
  bins=bins[bins@width<=max.binsize]
  bins=bins[bins@width>min.binsize]
  ###
  cat("done\n")
  bins=GenomicRanges::GRanges(seqnames = bins@NAMES, ranges= bins)
  c_pos<- c_filt[GenomicRanges::countOverlaps(c_filt, bins)>0]
  return(list("bins"= bins, "c_pos"=c_pos))
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


findBins <- function(gr,
                     n)
{
  cat("processing", as.character(gr@seqnames@values), "\n")
  return(IRanges::IRanges(start = ranges(gr)@start[c(1:(length(gr)-n+1))],end=(ranges(gr)@start[c(n:length(gr))])+1,
                          names = rep(gr@seqnames@values, times= length(gr)-n+1)))
}
