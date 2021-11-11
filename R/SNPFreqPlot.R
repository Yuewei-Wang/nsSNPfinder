#' Generate an overview plot of SNP locations versus number of nucleotide
#' variants by gene range
#'
#' A function that integrate the information from bioMart and plot the
#' SNP positions over the protein encoding range in
#' given chromosome, gene range (start coordinate, end coordinate).
#'
#' @param chrName A integer value of class "numeric" indicating
#'    the human chromosome 1-22, or a char either in 'X'or 'Y' indicating
#'    human sex chromosome.
#' @param start_position A positive integer indicating the starting coordinate
#'    of gene range.
#' @param end_position A positive integer indicating the end coordinate
#'    of gene range.
#'
#' @return Returns a plot indicating the overview SNPs
#'
#' @examples
#' # Calculates percentage of nsSNPs in chromosome 1gene range 2321253 to 2391707
#' SNPplot <- SNPFreqPlot(chrName = 1,
#'                        start_position = 2321253,
#'                        end_position = 2391707)
#' SNPplot
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
#' @import ggplot2
#' @import SNPlocs.Hsapiens.dbSNP144.GRCh38
#' @import BSgenome.Hsapiens.UCSC.hg38


# 09 Nov 2021

SNPFreqPlot <- function(chrName, start_position, end_position){

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

  # Set uo appropriate mart
  hsapiens <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

  # Query the information of transcripts and SNPs regarding the input range
  all_genes_Info <-biomaRt::getBM(attributes = c("hgnc_symbol","ensembl_gene_id",
                                                 "transcript_start","transcript_end",
                                                 "transcript_biotype"),
                                  filters = c('chromosome_name','start','end'),
                                  values = list(chrName, start_position, end_position),
                                  mart = hsapiens)
  all_chr_snp<-BSgenome::snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh38, as.character(chrName))

  # Filter out the sequence that not encoding with protein
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
    if (all_genes_Info$transcript_biotype[i] == 'protein_coding'){
      pc <- c(pc, i)
      gene_Name <- c(gene_Name, all_genes_Info$hgnc_symbol[i])
      len <- all_genes_Info$transcript_end[i] - all_genes_Info$transcript_start[i] + 1
      lengths <- c(lengths, len)
    }
  }
  if (length(pc) == 0){
    stop("no encoding protein involvoed in the input range, please re-enter another range.")
  }


  # Find the most wide width of encoded protein transcripts
  cur_max <- 1
  for (i in 2:length(lengths)){
    if (lengths[cur_max] < length(i)){
      cur_max <- i
    }
  }


  # Obtain corresponding gene sequence
  # record all within SNP locations and number of single substitute nucleotides
  # within the transcripts
  hs_seq <- BSgenome.Hsapiens.UCSC.hg38
  chr <- paste('chr',as.character(chrName),sep="")
  t_start <- all_genes_Info$transcript_start[pc[cur_max]]
  t_end <- all_genes_Info$transcript_end[pc[cur_max]]
  t_seq <-getSeq(x = hs_seq, names = chr, start = t_start, end = t_end)
  loc <- c()
  num_of_var <- c()
  for (i in 1:length(snps)){
    pos <- all_chr_snp@ranges@pos[snps[i]]
    if (pos <= t_end && pos >= t_start){
      loc <-c(loc, all_chr_snp@ranges@pos[snps[i]])
      iupac <- unlist(strsplit(IUPAC_CODE_MAP[all_chr_snp
                                              @elementMetadata
                                              $alleles_as_ambig[snps[i]]],""))
      num_of_var <- c(num_of_var, length(iupac))
    }
  }

  # Fill the information into the plot
  data <- data.frame(loc, num_of_var)
  plot<-ggplot(data = data,mapping = aes(x = loc, y = num_of_var),
                 color=num_of_var,fill=num_of_var) +
    geom_col(alpha=0.25) +
    geom_text(aes(label = loc),angle = 60, vjust = 0.25,
              size = 2.5, colour = "black") +
    labs(x = "SNP location", y = "Count of different nts",
         title = "SNP locations VS Number of variant nucleotides",
         subtitle = paste('Gene Name',gene_Name[cur_max],sep=": ")) +
    xlim(start_position, end_position) + ylim(c(0,4))


  return(plot)
}
# [END]
