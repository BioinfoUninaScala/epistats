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

`epiAnalysis` is the function that allows the extraction of the epiallele composition information for each genomic interval given as input (obtained trough the functions described above). A GRanges object containing all the intervals designed thorugh the previous steps is required. The user can adjust different parameters to perform a customized analysis. It can be decided to perform a stranded analysis or not. In this case, note that if you want to analyse non-CG methylation, only the stranded method is possible.
The check of bisulfite efficiency is also possible; the user can decide to exclude reads which show a low bisulfite efficiency from the analysis. The function also include an option to exclude reads with ambiguity.
The function returns a list containing three different objects:
* 1. a compressed binary matrix containing the epiallele composition for each analysed genomic region.

* 2. a summary dataframe containg the coordinates of the analysed regions and its relative metrics, such as mean CpG distance, average DNA methylation, etc...

* 3. a text file containing the coordinated of genomic intervals excluded from the analysis.


**4. Statistics**

EpiStatProfiler is provided with a set of functions that perform statistycal test on the output obtained from the previous functions.

* `epiStat` : this function takes as input the compressed matrices containing the epialleles counts of all the analysed samples. It perfoms the PERMANOVA test to identify those genomic loci which differ among groups for their epialleles composition.

* `diffStat` : this function takes as input the tables containing the summary metrics related to the analysed genomic loci in all samples. It implements non-parametric tests (Wilcoxon and Kruskal-Wallis tests) to identify significantly different genomic regions using the user-specified statistic (Shannon entropy, Average methylation, etc...)
</p>

