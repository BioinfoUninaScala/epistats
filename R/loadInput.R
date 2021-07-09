#' loadInput
#'
#' @importFrom Biostrings readDNAStringSet
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom GenomeInfoDb seqlevels
#' @param bamfile path to bamfile
#' @param genomefile path to genome/fasta
#' @return data
#' @export

loadInput=function(bamfile,
                   genomefile)
{
  algn=readBAM(bamfile)
  Genome=readGenome(genomefile)
  if (length(intersect(seqlevels(algn), names(Genome))) == length(Genome))
  {
    cat("chromosome names in input files match.\n Loading data done!\n")
    return(list("algn"=algn, "genome"=Genome))
  }else{
    cat("chromosome names in input files do not correspond.\n Matching chromosome names...\n")
    if(!length(names(Genome)) == length(seqlevels(algn)))
    {
      cat("Error: number of chromosomes in Genome not match number of chromosomes in Bam File.\n Match of chromosome names impeeded!\n")
      return(list("algn"=algn, "genome"=Genome))
    }else{
      Genome=matchGenomeWithBAM(Genome,algn)
      cat("match of chromosome names done.\n Loading data done!\n")
      return(list("algn"=algn, "genome"=Genome))
    }
  }
}

### Read Genome ----------------------------------------------------------------
readGenome <- function(GenomePath)
{
  cat("loading Genome from ", GenomePath, "...\n")
  Genome <- Biostrings::readDNAStringSet(GenomePath, format = "fasta")
  return(Genome)
}

### Read BAM -------------------------------------------------------------------
readBAM <- function(bamfile)
{
  cat("loading data from ", bamfile, "...\n")
  param= ScanBamParam(what = c("seq","strand"))
  algn <- GenomicAlignments::readGAlignments(bamfile,
                                             use.names=TRUE,
                                             param=param)
  return(algn)
}

### match Genome names with BAM names ------------------------------------------
matchGenomeWithBAM <- function(Genome,
                               algn)
{
  names(Genome)=GenomeInfoDb::seqlevels(algn)
  return(Genome)
}


