<p align="right">
 <img src="https://github.com/BioinfoUninaScala/epistats/blob/main/data-raw/logofinal.png" width="150" alt="EpiStatProfiler Logo">
</p>


## EpiStatProfiler
#### A new R package for the qualitative analysis of DNA methylation

### Introduction 

<p align = "justify"> 
 EpiStatProfiler is a new R package aimed to extract epialleles information from any type of bisulfite-sequencing data (from targeted to genome-wide experimental data). The epiallele-based analysis (EBA) relies on the characterization of the specific methylation patterns present on each sequenced read coming from bisulfite-sequencing experiments. 
 Given a genomic locus containing n Cs in the CpX context, all the possible combinations of the methylation states of these Cs observed at the single reads level are defined as epialleles. 
 Considering each sequenced read as coming from a single cell, this type of analysis can provide additional insights about the epigenetic cellular heterogeneity characterizing a sample. 
 EpiStatProfiler is implemented with the ability to extract the epialleles information from genomic units that can be easily designed by the user modifying several parameters. Beyond the profiling of the epialleles composition for each analysed sample, the package is then provided with dedicated statistical functions aimed to compare epialleles compositions among different biological conditions. Finally, epialleles extraction functions are designed to easily perform strand-specific and non-CG methylation analysis. 
</p>

---------

### Installation 
In R console, run 

```r
library(devtools)
install_github("BioinfoUninaScala/epistats", 
               build_vignettes=FALSE, 
               repos=BiocManager::repositories(),
               dependencies=TRUE, type="source")
```
----------

### Usage 

#### Main functions

* `loadInput`
* `filterByCoverage`
* `makeBins`
* `makeWindows`
* `epiAnalysis`
* `epiStat`
* `diffStat`

See [vignettes](/vignettes/my-vignette.Rmd) for details.
