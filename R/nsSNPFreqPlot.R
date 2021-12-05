#' Generate an overview plot for nsSNP locations versus number of nucleotide
#' variants by gene range
#'
#' A function that integrate the information from bioMart and plot the potential
#' nsSNP positions over the protein encoding range in
#' given chromosome, gene range (start coordinate, end coordinate).
#' For better resolution, the genome range is limited within 200-nt length.
#'
#' @param chrName A integer value of class "numeric" indicating
#'    the human chromosome 1-22, or a char either in 'X'or 'Y' indicating
#'    human sex chromosome.
#' @param startPosition A positive integer indicating the starting coordinate
#'    of gene range.
#' @param endPosition A positive integer indicating the end coordinate
#'    of gene range.
#'
#' @return Returns a plot indicating the overview potential nsSNPs
#'
#' @examples
#' # The used dataset of SNPs is default from SNPlocs.Hsapiens.dbSNP144.GRCh38 package
#' # The used dataset of Human genome is default from BSgenome.Hsapiens.UCSC.hg38 package
#'
#' # Generate the nsSNP locations in chromosome 3 in region of 49395439 to 49395566
#' nsSNPFreqPlot(chrName = 3,
#'               startPosition = 49395520,
#'               endPosition = 49395566)
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


nsSNPFreqPlot <- function(chrName, startPosition, endPosition){

  # Checking argument validation and valid range
  num_chroms <- c(1:22)
  chr_chroms <- c('X','Y')
  if (typeof(chrName) == "character" && !chrName %in% chr_chroms) {
    stop("Please use integer between 1 to 22 to express the human chromosome
         name, or character 'X' or 'Y' as human sex chromosome.
         Please re-enter a valid input.")
  }

  if (typeof(chrName) == "double" && !chrName %in% num_chroms){
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
  # Filter out the sequence that not encoding with protein
  snps <- c()
  for (i in 1:allChrSNP@elementMetadata@nrows){
    if (allChrSNP@ranges@pos[i] <= endPosition &&
        allChrSNP@ranges@pos[i] >= startPosition){
      snps <- append(snps, i)
    }
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
  snp <- c()
  numVars <- c()
  for (i in 1:length(snps)){
    pos <- allChrSNP@ranges@pos[snps[i]]
    if (pos <= tEnd && pos >= tStart){
      loc <-c(loc, allChrSNP@ranges@pos[snps[i]])
      localIndex <- pos-tStart+1
      if (localIndex %% 3 == 0){
        start <- localIndex - 2
        end <- localIndex
      } else{
        if (localIndex %% 3 == 1){
          start <- localIndex
          end <- localIndex + 2
        } else {
          start <- localIndex - 1
          end <- localIndex + 1
        }
      }
      if (end > length(tSeq)){
        break
      }
      iupac <- unlist(strsplit(IUPAC_CODE_MAP[allChrSNP
                                              @elementMetadata
                                              $alleles_as_ambig[snps[i]]],
                               split=""))
      oriCodon <- paste0(tSeq[start:end], collapse = "")
      acc <- 0
      for (k in 1:length(iupac)){
        replacedSeq <- replace(tSeq, localIndex, iupac[k])
        snpCodon <- paste0(replacedSeq[start:end], collapse = "")
        if (Biostrings::GENETIC_CODE[oriCodon] != Biostrings::GENETIC_CODE[snpCodon]){
          acc <- acc + 1
        }
      }
      if (acc > 0){
        numVars <- c(numVars, acc)
      } else {
        snp <- c(snp, i)
      }
    }
  }
  # Check whether nsSNP present, if no any nsSNP, gives message to user
  if (0 == length(numVars)){
    stop("no nsSNPs involvoed in the input range, you may re-define the range.")
  }


  # Fill the information into the plot
  loc <- loc[-snp]
  data <- data.frame(loc, numVars)
  plot<-ggplot(data = data, mapping = aes(x = loc, y = numVars),
               color=numVars,fill=numVars) +
    geom_col(alpha=0.25, width = 2) +
    geom_text(aes(label = loc),angle = 60, vjust = 0.25,
              size = 2.5, colour = "black") +
    labs(x = "nsSNP location", y = "Count of different nts",
         title = "nsSNP locations VS Number of variant nucleotides",
         subtitle = paste('Gene Name',geneName[curMax],sep=": ")) +
    xlim(startPosition, endPosition) + ylim(c(0,4))


  return(plot)
}

# [END]
