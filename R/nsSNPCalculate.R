#' Calculates percentage of nsSNPs by gene range
#'
#' A function that integrate the information from bioMart and calculates the
#' percentage of positions could involve nsSNP of protein encoding range in
#' given chromosome, gene range (start coordinate, end coordinate).
#'
#' @param chrName A integer value of class "numeric" indicating
#'    the human chromosome 1-22, or a char either in 'X' or 'Y' indicating
#'    human sex chromosome.
#' @param start_position A positive integer indicating the starting coordinate
#'    of gene range.
#' @param end_position A positive integer indicating the end coordinate
#'    of gene range.
#'
#' @return Returns a data frame object including the query and evaluated results
#' \itemize{
#'   \item gene_Name - A value of class "character" indicating the gene name
#'   \item lengths - A value of integer indicating the transcript length
#'   \item all_nssnps - A value of class "numeric" indicating the nsSNP-involved
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
#' nsSNPResults <- nsSNPCalculatebyRange(
#'                           chrName = 1,
#'                           start_position = 2321253,
#'                           end_position = 2391707)
#' nsSNPResults
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


nsSNPCalculatebyRange <- function(chrName, start_position, end_position){

  # Checking argument validation
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

  if (typeof(start_position) != "double" | typeof(end_position) != "double" |
      start_position < 0 | end_position < 0) {
    stop("start and end coordinates should be in positive integer
         of class numeric.")
  }

  # Set up appropriate mart to extract information of transcripts
  hsapiens <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

  # Query the information of transcripts and SNPs regarding the input range
  all_genes_Info <-biomaRt::getBM(attributes = c("hgnc_symbol","ensembl_gene_id",
                                        "transcript_start","transcript_end",
                                        "transcript_biotype","strand"),
                         filters = c('chromosome_name','start','end'),
                         values = list(chrName, start_position, end_position),
                         mart = hsapiens)
  all_chr_snp<-BSgenome::snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh38, as.character(chrName))

  # Filter out the query results that are not within the input range
  snps <- c()
  for (i in 1:all_chr_snp@elementMetadata@nrows){
    if (all_chr_snp@ranges@pos[i] <= end_position &&
        all_chr_snp@ranges@pos[i] >= start_position){
      snps <- append(snps, i)
    }
  }
  pc <- c()
  gene_Name <- c()
  lengths <- c()
  for (i in 1:nrow(all_genes_Info)){
    if (all_genes_Info$transcript_biotype[i] == 'protein_coding' &&
        all_genes_Info$transcript_start[i] >= start_position &&
        all_genes_Info$transcript_end[i] <= end_position){
      pc <- c(pc, i)
      gene_Name <- c(gene_Name, all_genes_Info$hgnc_symbol[i])
      len <- all_genes_Info$transcript_end[i] - all_genes_Info$transcript_start[i] + 1
      lengths <- c(lengths, len)
    }
  }

  # No encoding protein, the query should terminate since no functional nsSNP at least
  if (length(pc) == 0){
    stop("no encoding protein involvoed in the input range, please re-enter another range.")
  }


  # Obtain corresponding gene sequence
  # Filter the potential nsSNPs positions within input coordinates
  hs_seq <- BSgenome.Hsapiens.UCSC.hg38
  chr <- paste('chr',as.character(chrName),sep="")
  all_nssnps <- c()
  for(i in 1:length(pc)) {
    t_start <- all_genes_Info$transcript_start[pc[i]]
    t_end <- all_genes_Info$transcript_end[pc[i]]
    if(all_genes_Info$strand[pc[i]] == -1){
      strand <- "-"
    } else {
      strand <- "+"
    }
    t_seq <-getSeq(x = hs_seq, names = chr,
                   start = t_start, end = t_end, strand = strand)
    acc <- 0
    for (j in 1:length(snps)){
      snp_pos <- all_chr_snp@ranges@pos[snps[j]]
      if (snp_pos <= t_end && snp_pos >= t_start){
        local_idx <- snp_pos-t_start+1
        if (local_idx %% 3 == 0){
          start <- local_idx - 2
          end <- local_idx
        } else{
          if (local_idx %% 3 == 1){
            start <- local_idx
            end <- local_idx + 2
          } else {
            start <- local_idx - 1
            end <- local_idx + 1
          }
        }
        if (end > length(t_seq)){
          break
        }
        iupac <- unlist(strsplit(IUPAC_CODE_MAP[all_chr_snp@elementMetadata$alleles_as_ambig[pc[i]]],split=""))
        ori_codon <- paste0(t_seq[start:end], collapse = "")
        for (k in 1:length(iupac)){
          replaced_seq <- replace(t_seq, snp_pos-t_start, iupac[k])
          snp_codon <- paste0(replaced_seq[start:end], collapse = "")
          if (Biostrings::GENETIC_CODE[ori_codon] != Biostrings::GENETIC_CODE[snp_codon]){
            acc <- acc + 1
            break
          }
        }
        }
    }
    if (acc > 0){
      all_nssnps <- c(all_nssnps, acc)
    }
  }

  percent <- c()
  for (i in 1:length(pc)){
    percent <- append (percent, sprintf("%0.4f", all_nssnps[i]/lengths[i]))
  }

  # Combine all results and calculate the percentage
  result <- data.frame(gene_Name, lengths, all_nssnps, percent)

  return(result)
}

# [END]
