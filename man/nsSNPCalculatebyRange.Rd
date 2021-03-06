% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsSNPCalculate.R
\name{nsSNPCalculatebyRange}
\alias{nsSNPCalculatebyRange}
\title{Calculates percentage of nsSNPs by gene range}
\usage{
nsSNPCalculatebyRange(chrName, startPosition, endPosition)
}
\arguments{
\item{chrName}{A integer value of class "numeric" indicating
the human chromosome 1-22, or a char either in 'X' or 'Y' indicating
human sex chromosome.}

\item{startPosition}{A positive integer indicating the starting coordinate
of gene range.}

\item{endPosition}{A positive integer indicating the end coordinate
of gene range.}
}
\value{
Returns a data frame object including the query and evaluated results
\itemize{
\item geneName - A value of class "character" indicating the gene name
\item lengths - A value of integer indicating the transcript length
\item nsSNPs - A value of class "numeric" indicating the nsSNP-involved
locations over the transcript
\item percent - A value of class "numeric" indicating the percentage of
nsSNP over the transcript
}
}
\description{
A function that integrate the information from bioMart and calculates the
percentage of positions could involve nsSNP of protein encoding range in
given chromosome, gene range (start coordinate, end coordinate).
}
\examples{
# The used dataset of SNPs is default from SNPlocs.Hsapiens.dbSNP144.GRCh38 package
# The used dataset of Human genome is default from BSgenome.Hsapiens.UCSC.hg38 package

# Calculates percentage of nsSNPs in chromosome 3 gene range 49359145 to 49411645
nsSNPCalculatebyRange(chrName = 3,
                     startPosition = 49359145,
                     endPosition = 49411645)

}
\references{
Durinck, S., Spellman, P., Birney, E.,& Huber, W. (2009). Mapping identifiers
for the integration ofgenomic datasets with the R/Bioconductor package
biomaRt. \emph{Nature Protocols}, 4, 1184–1191.2.

Durinck, S., Moreau, Y., Kasprzyk, A., Davis, S., De Moor, B., Brazma,
A.,& Huber, W. (2005).BioMart and Bioconductor: a powerful link between
biological databases and microarray data analysis.\emph{Bioinformatics},
21, 3439–3440.

Pagès, H., Aboyoun, P., Gentleman, R., DebRoy, S. (2021). Biostrings:
Efficient manipulation of biological strings. R package version 2.62.0,
https://bioconductor.org/packages/Biostrings.

Pagès, H. (2017). SNPlocs.Hsapiens.dbSNP144.GRCh38: SNP locations for Homo
sapiens (dbSNPBuild 144). R package version 0.99.20.
}
