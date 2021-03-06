% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/6_epiAnalysis.R
\name{epiAnalysis}
\alias{epiAnalysis}
\title{Extracting epiallele information from customised target regions.}
\usage{
epiAnalysis(
  align,
  bin,
  aname,
  threshold,
  pathDir,
  genome,
  myfuns,
  bisu.Thresh = 0,
  stranded = FALSE,
  mode = "CG",
  remove.Amb = TRUE,
  cores = 1,
  retain.reads = FALSE,
  get.cPos = FALSE
)
}
\arguments{
\item{align}{A GAlignments object obtained trough the loading and filtering steps.}

\item{bin}{A GRanges object that is the output obtained from one of the functions makeWindows() and makeBins().}

\item{aname}{A character indicating the analysis name ("SampleName")}

\item{threshold}{An integer indicating the minimum coverage required to analyse the interval. This is specified here since it could be possible that DNA methylation status cannot be assessed for some reads, thus decreasing altering the coverage of the previously selected regions.}

\item{pathDir}{A string indicating the path of the directory desired to save the output.}

\item{genome}{A DNAStringSet object containing the fasta genome}

\item{myfuns}{A list object containing the functions to be applied to the binary epiallele matrix (to compute average DNA methylation, Shannon Entropy, etc...).}

\item{bisu.Thresh}{An integer indicating the bisulfite efficiency threshold to be used to discard reads or not (Default = 0).}

\item{stranded}{Logical indicating whether the analysis should be performed stranded or not (Default = FALSE).}

\item{mode}{A character indicating the pattern to be used to perform the methylation analysis (one of "CG", "CA", "CC", "CT").}

\item{remove.Amb}{Logical indicating whether ambiguous reads (containing not expected nucleotides) should be discarded or not (Default = TRUE).}

\item{cores}{An integer indicating the number of cores to be used for the computation.}

\item{retain.reads}{Logical indicating whether reads that have a low bisulfite conversion efficiency should be kept or not (Default = TRUE).}

\item{get.cPos}{Logical indicating whether the coordinates of the CpGs displayed through the analysis should be returned as a text output or not (Default = FALSE).}
}
\value{
A list object containing different dataframes.

intervals = Dataframe having all the analysed regions as rows and the different summary statistics as columns.

epi = Dataframe that is the compressed form of the binary matrix containing the epiallele composition for each analysed region.

log = Dataframe containing the coordinates of the genomic regions that have been discarded for the analysis.
}
\description{
This is the function which takes as input the user-defined regions of interest and it returns for each of them epiallele composition and summary statistics, such as average DNA methylation, Shannon Entropy, etc...
}
\examples{
data <- loadInput(bamfile, genomefile)
algn <- data[[1]]
Genome <- data [[2]]

## Keeping only standard chromosomes
filtered <- filterBAM(algn,
                      keepStChr =TRUE,
                      keepM = FALSE,
                      retainChr = NULL)

## Selecting regions with a minimum coverage
covered <- filterByCoverage(algn, threshold = 50, minsize = 8)

## Choose one of the two available methods to design target regions.
windows <- makeWindows(gr = covered,
                       Genome = Genome,
                       window = 50,
                       step = 1,
                       mode = "CG",
                       min.C = 2,
                       max.C = 50,
                       cores = 40)

## Extract the GRanges object containing the target regions
regions <- windows$windows

## Select the functions that work on the binary epiallele matrix
myfuns <- list("dist" = cdist,
               "epi" = epi,
               "singleton" = singleton,
               "maxfreq" = maxfreq,
               "shannon" = shannon,
               "mean_met" = meanMeth,
               "num_cg" = numCG,
               "num_reads" = num_reads)

## Extract epiallele information
aname = "Sample_1"

out <- epiAnalysis(align = algn,
                   bin = regions,
                   aname = aname,
                   threshold = 50,
                   bisu.Thresh = 0,
                   stranded = FALSE,
                   mode = "CG",
                   remove.Amb = TRUE,
                   genome = Genome,
                   pathDir = "/home/sarnataro/bam/out_epistats/",
                   cores = cores,
                   retain.reads = TRUE,
                   get.cPos = FALSE,
                   myfuns = myfuns)
}
