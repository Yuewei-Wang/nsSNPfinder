#' Generate a plot of SNP locations versus number of nucleotide
#' variants by gene range
#'
#' A function that integrate the information from bioMart and plot the
#' SNP positions over the protein encoding range in
#' given chromosome, gene range (start coordinate, end coordinate).
#' For better resolution, the genome range is limited within 200-nt length.
#'
#' @param chrName A integer value of class "numeric" indicating
#'    the human chromosome 1-22, or a char either in 'X' or 'Y' indicating
#'    human sex chromosome.
#' @param startPosition A positive integer indicating the starting coordinate
#'    of gene range.
#' @param endPosition A positive integer indicating the end coordinate
#'    of gene range.
#'
#' @return Returns a plot indicating the overview situation of SNPs
#'
#' @examples
#' # The used dataset of SNPs is default from SNPlocs.Hsapiens.dbSNP144.GRCh38 package
#' # The used dataset of Human genome is default from BSgenome.Hsapiens.UCSC.hg38 package
#'
#' # Generate the SNP locations in chromosome 3 in region of 49395438 to 49395450
#' SNPFreqPlot(chrName = 3,
#'            startPosition = 49395438,
#'            endPosition = 49395566)
#'
#' @references
#' Durinck, S., Spellman, P., Birney, E.,& Huber, W. (2009). Mapping identifiers
#'  for the integration ofgenomic datasets with the R/Bioconductor package
#'  biomaRt. *Nature Protocols*, 4, 1184–1191.2.
#'
#' Durinck, S., Moreau, Y., Kasprzyk, A., Davis, S., De Moor, B., Brazma,
#' A.,& Huber, W. (2005).BioMart and Bioconductor: a powerful link between
#' biological databases and microarray data analysis.*Bioinformatics*,
#' 21, 3439–3440.
#'
#' Pagès, H., Aboyoun, P., Gentleman, R., DebRoy, S. (2021). Biostrings:
#' Efficient manipulation of biological strings. R package version 2.62.0,
#'  https://bioconductor.org/packages/Biostrings.
#'
#' Pagès, H. (2017). SNPlocs.Hsapiens.dbSNP144.GRCh38: SNP locations for Homo
#' sapiens (dbSNPBuild 144). R package version 0.99.20.
#'
#' Wickham, H. (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag
#' New York. ISBN 978-3-319-24277-4, https://ggplot2.tidyverse.org.
#'
#' @export
#' @import BSgenome
#' @import Biostrings
#' @import biomaRt
#' @import ggplot2
#' @import SNPlocs.Hsapiens.dbSNP144.GRCh38
#' @import BSgenome.Hsapiens.UCSC.hg38


SNPFreqPlot <- function(chrName, startPosition, endPosition){

  # Checking arguments and the input range validation
  numChroms <- c(1:22)
  chrChroms <- c('X','Y')
  if (typeof(chrName) == "character" && !chrName %in% chrChroms) {
    stop("Please use integer between 1 to 22 to express the human chromosome
         name, or character 'X' or 'Y' as human sex chromosome.
         Please re-enter a valid input.")
  }

  if (typeof(chrName) == "double" && !chrName %in% numChroms){
    stop("The numeric input for human chromosome is between 1-22.
         Please re-enter a valid input.")
  }

  if (typeof(startPosition) != "double" | typeof(endPosition) != "double" |
      startPosition < 0 | endPosition < 0) {
    stop("start and end coordinates should be in positive integer
         of class numeric.")
  }

  if (endPosition - startPosition > 200) {
    stop("The length of sequence input is over 200, please refine the range
         within 200-length for better resolution.")
  }

  # Query the information of transcripts and SNPs regarding the input range
  allGenesInfo <-findGeneInfo(chrName, startPosition, endPosition)
  allChrSNP<-BSgenome::snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh38,
                                     as.character(chrName))
  if (nrow(allGenesInfo) == 0){
    stop("no transcpits available for the input coordinates.")
  }

  # Filter out the sequence that not encoding with protein or exceeds the range
  snps <- c()
  for (i in 1:allChrSNP@elementMetadata@nrows){
    if (allChrSNP@ranges@pos[i] <= endPosition &&
        allChrSNP@ranges@pos[i] >= startPosition){
      snps <- append(snps, i)
    }
  }
  if (length(snps) == 0){
    stop("no SNPs involvoed in the input range.")
  }

  pc <- c()
  geneName <- c()
  lengths <- c()
  for (i in 1:nrow(allGenesInfo)){
    if (allGenesInfo$transcript_biotype[i] == 'protein_coding'){
      pc <- c(pc, i)
      geneName <- c(geneName, allGenesInfo$hgnc_symbol[i])
      len <- allGenesInfo$transcript_end[i] - allGenesInfo$transcript_start[i] + 1
      lengths <- c(lengths, len)
    }
  }
  if (length(pc) == 0){
    stop("no encoding protein involvoed in the input range, please re-enter another range.")
  }


  # Find the most wide width of encoded protein transcripts
  curMax <- 1
  for (i in 2:length(lengths)){
    if (lengths[curMax] < length(i)){
      curMax <- i
    }
  }


  # Obtain corresponding gene sequence
  # record all within SNP locations and number of single substitute nucleotides
  # within the transcripts
  #hsapiensSeq <- BSgenome.Hsapiens.UCSC.hg38
  chr <- paste('chr',as.character(chrName),sep="")
  tStart <- allGenesInfo$transcript_start[pc[curMax]]
  tEnd <- allGenesInfo$transcript_end[pc[curMax]]
  tSeq <-getSeq(x = hsapiensSeq, names = chr, start = tStart, end = tEnd)
  loc <- c()
  numVars <- c()
  for (i in 1:length(snps)){
    pos <- allChrSNP@ranges@pos[snps[i]]
    if (pos <= tEnd && pos >= tStart){
      loc <-c(loc, allChrSNP@ranges@pos[snps[i]])
      iupac <- unlist(strsplit(IUPAC_CODE_MAP[allChrSNP
                                              @elementMetadata
                                              $alleles_as_ambig[snps[i]]],""))
      numVars <- c(numVars, length(iupac))
    }
  }

  # Fill the information into the plot
  data <- data.frame(loc, numVars)
  plot<-ggplot(data = data,mapping = aes(x = loc, y = numVars),
                 color=numVars,fill=numVars) +
    geom_col(alpha=0.25, width = 2) +
    geom_text(aes(label = loc),angle = 60, vjust = 0.25,
              size = 2.5, colour = "black") +
    labs(x = "SNP location", y = "Count of different nts",
         title = "SNP locations VS Number of variant nucleotides",
         subtitle = paste('Gene Name',geneName[curMax],sep=": ")) +
    xlim(startPosition, endPosition) + ylim(c(0,4))


  return(plot)
}

# [END]
