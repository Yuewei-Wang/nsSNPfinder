
<!-- README.md is generated from README.Rmd. Please edit that file -->

# nsSNPfinder

<!-- badges: start -->
<!-- badges: end -->

## Description

nsSNPfinder is an R package for integrating, analyzing, and visualizing
the SNP data from the genomic level to the protein level, finding
potential nsSNP-associated gene marker that might cause the protein
structure change. This package was developed using `R version 4.1.1` and
`x86_64-apple-darwin17.0 (64-bit)` platform runnning under
`macOS Monterey 12.0.1`.

## Installation

You can install the development version of nsSNPfinder from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
require("devtools")
devtools::install_github("Yuewei-Wang/nsSNPfinder", build_vignettes = TRUE)
library(nsSNPfinder)
```

To run the shinyApp:

``` r
runNSSNPfinder()
```

## Overview

``` r
ls("package:nsSNPfinder")
```

nsSNPfinder contains 4 functions to demonstrate the nsSNP distribution
in certain H. sapiens gene and address the protein structure of the
gene.

The *nsSNPCalculatebyRange* function calculates the percentage of
residues involving nsSNP event in the gene region.

The *SNPFreqPlot* function generates the plot to indicate the position
distribution of SNPs within the certain gene region.

![](./inst/extdata/overview.png)

The *nsSNPFreqPlot* function generates the plot similar as *SNPFreqPlot*
but the targets are all nsSNPs within the gene.

The *displayPDB* function displays the encoded-protein 3D structure
record in PBD or UniProt, and highlight the residue positions involved
nsSNPs.

``` r
browseVignettes("nsSNPfinder")
```

## Contributions

The author of the package is Yuewei Wang. The *nsSNPCalculatebyRange*
function utilized the `SNPlocs.Hsapiens.dbSNP144.GRCh38` and
`BSgenome.Hsapiens.UCSC.hg38` to obtain the initial information of SNP
locations and human chromosome sequences. The `biomaRt` package is used
for collecting the genomic sequences and type of variants.
`GENETIC_CODE` function in `Biostrings` package is used for compare the
codon, which contributes to find potential nsSNPs.

The *nsSNPFreqPlot* and *SNPFreqPlot* functions use map function from
`ggplot2` R package to generate frequency plots.

The *displayPDB* function use the function from `prort` and model
visualization function from `r3dmol` R packages to obtain the Uniprot
recorded protein sequence and pdb file of structure modelling..

## Reference

Durinck, S., Spellman, P., Birney, E.,& Huber, W. (2009). Mapping
identifiers for the integration ofgenomic datasets with the
R/Bioconductor package biomaRt. *Nature Protocols*, 4, 1184–1191.2.

Durinck, S., Moreau, Y., Kasprzyk, A., Davis, S., De Moor, B., Brazma,
A.,& Huber, W. (2005).BioMart and Bioconductor: a powerful link between
biological databases and microarray data analysis.*Bioinformatics*, 21,
3439–3440.

Pagès, H., Aboyoun, P., Gentleman, R., DebRoy, S. (2021). Biostrings:
Efficient manipulation of biological strings. R package version 2.62.0.

Pagès, H. (2017). SNPlocs.Hsapiens.dbSNP144.GRCh38: SNP locations for
Homo sapiens (dbSNPBuild 144). R package version 0.99.20.

Pagès, H. (2021). BSgenome: Software infrastructure for efficient
representation of full genomes and their SNPs. R package version 1.62.0.

Su, W (2020). r3dmol: Create Interactive 3D Visualizations of Molecular
Data. R package version 0.1.0. <https://github.com/swsoyee/r3dmol>

Team TBD (2021). BSgenome.Hsapiens.UCSC.hg38: Full genome sequences for
Homo sapiens (UCSC version hg38, based on GRCh38.p13). R package version
1.4.4.

Wickham, H. (2016). ggplot2: Elegant Graphics for Data Analysis.
Springer-Verlag New York. ISBN 978-3-319-24277-4,
<https://ggplot2.tidyverse.org>.

Xiao, N., Cao, DS., Zhu, MF., Xu, QS. (2015). protr/ProtrWeb: R package
and web server for generating various numerical representation schemes
of protein sequences. *Bioinformatics* 31 (11), 1857-1859.

## Acknowledgements:

This package was developed as part of an assessment for 2021 BCB410H:
Applied Bioinformatics, University of Toronto, Toronto, CANADA.

## Q&A

-   How to solve the *timeout error* in installation:

Due to the size of data packages (BSgenome.Hsapiens.UCSC.hg38 and
SNPlocs.Hsapiens.dbSNP144.GRCh38), sometimes the timeout error occurs in
the process of installing packages. You can increase the timeout setting
by:

``` r
library("utils")
# set to 300s (5mins)
options(timeout = 300)
```

to solve the error. Then running installation would work.
