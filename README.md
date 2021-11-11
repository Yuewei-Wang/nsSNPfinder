
<!-- README.md is generated from README.Rmd. Please edit that file -->

# nsSNPfinder

<!-- badges: start -->
<!-- badges: end -->

nsSNPfinder is an R package for integrating, analyzing, and visualizing
the SNP data from the genomic level to the protein level, finding
potential nsSNP-associated gene marker that might cause the protein
structure change.

## Installation

You can install the development version of nsSNPfinder from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
require("devtools")
devtools::install_github("Yuewei-Wang/nsSNPfinder"，build_vignettes = TRUE)
library(nsSNPfinder)
```

To run the shinyApp: Under construction

## Overview

``` r
ls("package:nsSNPfinder")
data(package = "nsSNPfinder") # optional
```

nsSNPfinder contains 4 functions to demonstrate the nsSNP distribution
in certain H. sapiens gene and address the protein structure of the
gene.

The *nsSNPCalculatebyRange* function calculates the percentage of
residues involving nsSNP event in the gene region.

The *nsSNPFreqPlot* function generates the plot to indicate the position
of nsSNP and frequency distribution within the certain gene region.

The *SNPFreqPlot* function generates the plot similar as *nsSNPFreqPlot*
but the targets are all SNPs within the gene.

The *3Dpdb* function displays the encoded-protein 3D structure record in
PBD or UniProt, and highlight the residue positions involved nsSNPs.

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

The *nsSNPFreqPlot* and *SNPFreqPlot* functions use of map function from
`ggplot2` R package to generate frequency plots. The `bio3d` R package
is used for generating 3D protein strucuture of the gene.

## Reference

Durinck, S., Spellman, P., Birney, E.,& Huber, W. (2009). Mapping
identifiers for the integration ofgenomic datasets with the
R/Bioconductor package biomaRt. *Nature Protocols*, 4, 1184–1191.2.

Durinck, S., Moreau, Y., Kasprzyk, A., Davis, S., De Moor, B., Brazma,
A.,& Huber, W. (2005).BioMart and Bioconductor: a powerful link between
biological databases and microarray data analysis.*Bioinformatics*, 21,
3439–3440.

Grant, B.J., Rodrigues, A.P.C., ElSawy, K.M., McCammon, J.A.,& Caves,
L.S.D. (2006). Bio3D: AnR package for the comparative analysis of
protein structures. *Bioinformatics*, 22, 2695–2696.

Pagès, H. (2017). SNPlocs.Hsapiens.dbSNP144.GRCh38: SNP locations for
Homo sapiens (dbSNPBuild 144). R package version 0.99.20.

## Acknowledgements:

This package was developed as part of an assessment for 2021 BCB410H:
Applied Bioinformatics, University of Toronto, Toronto, CANADA.
