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

### Workflow 

The functions that are provided by the tool can be grouped into 4 modules: 

**1. Loading functions**

Input data consist of the alignment file and a reference genome provided as BAM and FASTA files, rispectively. 
The paths of these files are to be specified as parameters in the loading function in order to be imported in the global environment. 
Input files are then loaded as GAlignment and DNAStringSet objects, rispectively.
You can get the single objects by subsetting them from the list returned by the `loadInput` function. 
Filtering functions allow to subset the aligned reads covering only desired regions. 

* `filterBAM` : it allows to subset from the GAlignments objects only reads mapped on the specified chromosomes.
* `filterByCoverage` : it allows to select those genomic regions which are covered by minimum number of reads.

**2. Target regions design**

EpiStatProfiler allows the customization of the genomic regions to be profiled. The user can design genomic regions using two different methods.
The figure below show the rationale that is behind these two methods. 

* `makeBins` : this function allows to build genomic regions having a different length size, but containing a fixed number of CpG sites (or any other pattern used to perform the analysis (e.g., "CA", "CT")).

* `makeWindows` : this function allows to build genomic regions by using sliding windows of fixed length and step size, containing a variable numer of CpX sites.


**3. Epialleles extraction**

**4. Statistics**



See [vignettes](/vignettes/my-vignette.Rmd) for details.
