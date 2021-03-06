% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/displayPDB.R
\name{displayPDB}
\alias{displayPDB}
\title{Display the 3D protein structure and nsSNP-involved residue}
\usage{
displayPDB(chrName, geneName, nsPos, choice = 0)
}
\arguments{
\item{chrName}{A integer value of class "numeric" indicating
the human chromosome 1-22, or a char either in 'X'or 'Y' indicating
human sex chromosome.}

\item{geneName}{A string that indicating the name of gene}

\item{nsPos}{A positive integer indicating the coordinate of interested
nsSNP coordinate within the sequence of geneName.}

\item{choice}{A integer value indicating the tag of pdb ID choice.
Defult value is 0, which means the function will automatically select
the first pdb ID. Users can make choice = 1 to specify that they would
like to choose the ID manually.}
}
\value{
Returns a 3D model indicating the potential nsSNP location within
the model
}
\description{
A function that shows the 3D protein structure and points the potential nsSNP
by query and integrate the information from bioMart, Uniprot and PDB. It
allows the user to enter the interested pdb ID from query list, then
provide the protein structure.
}
\examples{
# The used dataset of SNPs is default from SNPlocs.Hsapiens.dbSNP144.GRCh38 package
# The used dataset of Human genome is default from BSgenome.Hsapiens.UCSC.hg38 package

# Display the 3D protein structure for RHOA in chromosome 3 at coordinate 49395565

# Example 1: default choice of pdb ID
displayPDB(chrName = 3,
      geneName = 'RHOA',
      nsPos = 49395565,
      choice = 0)

\dontrun{
# Example 2: choose your own pdb ID
structure <- displayPDB(chrName = 3,
                   geneName = 'RHOA',
                   nsPos = 49395565,
                   choice = 1)
# [1] "The available pdb structrues associated with your genes are:"
# [1] "4D0N" "6BCA" "4XOI" "5JCP" "6BC0" "4XH9" "1XCG" "3KZ1" "3T06" "5JHG"
# "5JHH" "1X86" "2RGN" "1CC0"
# [15] "5ZHX" "5HPY" "1CXZ" "5C2K" "1OW3" "1TX4" "5M6X" "5M70" "6R3V" "3MSX"
# "1A2B" "1DPF" "1FTN" "1KMQ"
# [29] "1LB1" "1S1C" "3LW8" "3LWN" "3LXR" "4XSG" "4XSH" "5A0F" "5BWM" "5C4M"
# "5EZ6" "5FR1" "5FR2" "5IRC"
# [43] "6BCB" "6KX2" "6KX3"
# please enter one ID to display
# Here requires the user to input the interested pdb ID from the displayed list
# Note that the prompt is not case-sensitive, which means '6bca' also works
6BCA

structure
}

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

Su, W. (2020). r3dmol: Create Interactive 3D Visualizations of Molecular Data.
R package version 0.1.0. https://github.com/swsoyee/r3dmol

Xiao, N., Cao, DS., Zhu, MF., Xu, QS. (2015). protr/ProtrWeb: R package and
web server for generating various numerical representation schemes of protein
sequences. \emph{Bioinformatics} 31 (11), 1857-1859.
}
