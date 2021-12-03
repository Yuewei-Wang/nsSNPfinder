#' Calculates percentage of nsSNPs by gene range
#'
#' A function that integrate the information from bioMart and calculates the
#' percentage of positions could involve nsSNP of protein encoding range in
#' given chromosome, gene range (start coordinate, end coordinate).
#'
#' @param chrName A integer value of class "numeric" indicating
#'    the human chromosome 1-22, or a char either in 'X' or 'Y' indicating
#'    human sex chromosome.
#' @param startPosition A positive integer indicating the starting coordinate
#'    of gene range.
#' @param endPosition A positive integer indicating the end coordinate
#'    of gene range.
#'
#' @return Returns a data frame object including the query and evaluated results
#' \itemize{
#'   \item geneName - A value of class "character" indicating the gene name
#'   \item lengths - A value of integer indicating the transcript length
#'   \item nsSNPs - A value of class "numeric" indicating the nsSNP-involved
#'    locations over the transcript
#'   \item percent - A value of class "numeric" indicating the percentage of
#'   nsSNP over the transcript
#' }
#'
#' @examples
#' # The used dataset of SNPs is default from SNPlocs.Hsapiens.dbSNP144.GRCh38 package
#' # The used dataset of Human genome is default from BSgenome.Hsapiens.UCSC.hg38 package
#'
#' # Calculates percentage of nsSNPs in chromosome 1 gene range 2321253 to 2391707
#' nsSNPCalculatebyRange(chrName = 1,
#'                      startPosition = 2321253,
#'                      endPosition = 2391707)
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
#' @export
#' @import BSgenome
#' @import Biostrings
#' @import biomaRt
#' @import SNPlocs.Hsapiens.dbSNP144.GRCh38
#' @import BSgenome.Hsapiens.UCSC.hg38


nsSNPCalculatebyRange <- function(chrName, startPosition, endPosition){

  # Checking argument validation
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

  # Query the information of transcripts and SNPs regarding the input range
  allGenesInfo <-findGeneInfo(chrName, startPosition, endPosition)
  allChrSNP<-BSgenome::snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh38,
                                     as.character(chrName))

  # Filter out the query results that are not within the input range
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
    if (allGenesInfo$transcript_biotype[i] == 'protein_coding' &&
        allGenesInfo$transcript_start[i] >= startPosition &&
        allGenesInfo$transcript_end[i] <= endPosition){
      pc <- c(pc, i)
      geneName <- c(geneName, allGenesInfo$hgnc_symbol[i])
      len <- allGenesInfo$transcript_end[i] - allGenesInfo$transcript_start[i] + 1
      lengths <- c(lengths, len)
    }
  }

  # No encoding protein, the query should terminate since no functional nsSNP at least
  if (length(pc) == 0){
    stop("no encoding protein involvoed in the input range, please re-enter another range.")
  }


  # Obtain corresponding gene sequence
  # Filter the potential nsSNPs positions within input coordinates
  #hsapiensSeq <- BSgenome.Hsapiens.UCSC.hg38
  chr <- paste('chr',as.character(chrName),sep="")
  nsSNPs <- c()
  for(i in 1:length(pc)) {
    tStart <- allGenesInfo$transcript_start[pc[i]]
    tEnd <- allGenesInfo$transcript_end[pc[i]]
    if(allGenesInfo$strand[pc[i]] == -1){
      strand <- "-"
    } else {
      strand <- "+"
    }
    tSeq <-getSeq(x = hsapiensSeq, names = chr,
                   start = tStart, end = tEnd, strand = strand)
    acc <- 0
    for (j in 1:length(snps)){
      snpPos <- allChrSNP@ranges@pos[snps[j]]
      if (snpPos <= tEnd && snpPos >= tStart){
        localIndex <- snpPos-tStart+1
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
        iupac <- unlist(strsplit(IUPAC_CODE_MAP[allChrSNP@elementMetadata$alleles_as_ambig[pc[i]]],split=""))
        oriCodon <- paste0(tSeq[start:end], collapse = "")
        for (k in 1:length(iupac)){
          replacedSeq <- replace(tSeq, snpPos-tStart, iupac[k])
          snpCodon <- paste0(replacedSeq[start:end], collapse = "")
          if (Biostrings::GENETIC_CODE[oriCodon] != Biostrings::GENETIC_CODE[snpCodon]){
            acc <- acc + 1
            break
          }
        }
        }
    }
    if (acc > 0){
      nsSNPs <- c(nsSNPs, acc)
    }
  }

  percent <- c()
  for (i in 1:length(pc)){
    percent <- append (percent, sprintf("%0.4f", nsSNPs[i]/lengths[i]))
  }

  # Combine all results and calculate the percentage
  result <- data.frame(geneName, lengths, nsSNPs, percent)

  return(result)
}
#' Helper
#'
#' load the retrieved information from the BioMart database according input
#' Note: the input should be valid since the preliminary check was done in the
#' main functions
#'
#' @param chrName the chromosome name of user input in the main function
#' @param startPosition the starting coordinate of user input in the main function
#' @param endPosition the end coordinate of user input in the main function
#' @import biomaRt
findGeneInfo <- function(chrName, startPosition, endPosition) {
  # for saving the run-time, here preload the example bm
  if (chrName == 1 && startPosition >= 2321253 && endPosition <= 2391707){
    return(demoAllGenesInfo)
  }
  bm <- biomaRt::getBM(attributes = c("hgnc_symbol","ensembl_gene_id",
                                "transcript_start","transcript_end",
                                "transcript_biotype","strand"),
                 filters = c('chromosome_name','start','end'),
                 values = list(chrName, startPosition, endPosition),
                 mart = hsapiens)
  return(bm)
}
# [END]
