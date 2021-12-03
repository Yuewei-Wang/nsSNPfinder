#' Display the 3D protein structure and nsSNP-involved residue
#'
#' A function that shows the 3D protein structure and points the potential nsSNP
#' by query and integrate the information from bioMart, Uniprot and PDB. It
#' allows the user to enter the interested pdb ID from query list, then
#' provide the protein structure.
#'
#' @param chrName A integer value of class "numeric" indicating
#'    the human chromosome 1-22, or a char either in 'X'or 'Y' indicating
#'    human sex chromosome.
#' @param geneName A string that indicating the name of gene
#' @param nsPos A positive integer indicating the coordinate of interested
#'    nsSNP coordinate within the sequence of geneName.
#' @param choice A integer value indicating the tag of pdb ID choice.
#'    Defult value is 0, which means the function will automatically select
#'    the first pdb ID. Users can make choice = 1 to specify that they would
#'    like to choose the ID manually.
#'
#' @return Returns a 3D model indicating the potential nsSNP location within
#'    the model
#'
#' @examples
#' # The used dataset of SNPs is default from SNPlocs.Hsapiens.dbSNP144.GRCh38 package
#' # The used dataset of Human genome is default from BSgenome.Hsapiens.UCSC.hg38 package
#'
#' # Display the 3D protein structure for RHOA in chromosome 3 at coordinate 49395565
#'
#' # Example 1: default choice of pdb ID
#' displayPDB(chrName = 3,
#'       geneName = 'RHOA',
#'       nsPos = 49395565,
#'       choice = 0)
#'
#' \dontrun{
#' # Example 2: choose your own pdb ID
#' structure <- displayPDB(chrName = 3,
#'                    geneName = 'RHOA',
#'                    nsPos = 49395565,
#'                    choice = 1)
#' # [1] "The available pdb structrues associated with your genes are:"
#' # [1] "4D0N" "6BCA" "4XOI" "5JCP" "6BC0" "4XH9" "1XCG" "3KZ1" "3T06" "5JHG"
#' # "5JHH" "1X86" "2RGN" "1CC0"
#' # [15] "5ZHX" "5HPY" "1CXZ" "5C2K" "1OW3" "1TX4" "5M6X" "5M70" "6R3V" "3MSX"
#' # "1A2B" "1DPF" "1FTN" "1KMQ"
#' # [29] "1LB1" "1S1C" "3LW8" "3LWN" "3LXR" "4XSG" "4XSH" "5A0F" "5BWM" "5C4M"
#' # "5EZ6" "5FR1" "5FR2" "5IRC"
#' # [43] "6BCB" "6KX2" "6KX3"
#' # please enter one ID to display
#' # Here requires the user to input the interested pdb ID from the displayed list
#' # Note that the prompt is not case-sensitive, which means '6bca' also works
#' 6BCA
#'
#' structure
#' }
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
#' Su, W. (2020). r3dmol: Create Interactive 3D Visualizations of Molecular Data.
#'  R package version 0.1.0. https://github.com/swsoyee/r3dmol
#'
#' Xiao, N., Cao, DS., Zhu, MF., Xu, QS. (2015). protr/ProtrWeb: R package and
#' web server for generating various numerical representation schemes of protein
#' sequences. *Bioinformatics* 31 (11), 1857-1859.
#'
#' @export
#' @import BSgenome
#' @import Biostrings
#' @import biomaRt
#' @import ggplot2
#' @import r3dmol
#' @import protr
#' @import magrittr
#' @import SNPlocs.Hsapiens.dbSNP144.GRCh38
#' @import BSgenome.Hsapiens.UCSC.hg38


displayPDB <- function(chrName, geneName, nsPos, choice = 0){

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

  if (typeof(geneName) != "character") {
    stop("Gene name should be string with quote.")
  }

  if (typeof(nsPos) != "double" && nsPos <= 0) {
    stop("The coordinate of nsSNP should be positive integer.")
  }

  if (typeof(choice) != "double" | (choice != 1 && choice != 0)) {
    stop("The pdb ID choice should be integer in either
         0 (always the first pdb ID) or 1 (choose my own).")
  }


  # Set up appropriate mart for human genome from bioMart
  allGenesInfo <-biomaRt::getBM(attributes = c("hgnc_symbol","ensembl_gene_id",
                                                 "transcript_start","transcript_end",
                                                 "transcript_biotype",
                                                 "uniprot_gn_id","pdb","strand"),
                                  filters = c('chromosome_name','hgnc_symbol'),
                                  values = list(chrName, geneName),
                                  mart = hsapiens)

  # choose the source of pdb structure
  pc <- c()
  uniprotIDs <- c()
  pdbIDs <- c()
  for (i in 1:nrow(allGenesInfo)){
    if (allGenesInfo$transcript_biotype[i] == 'protein_coding' &&
        allGenesInfo$pdb[i] != "" &&
        allGenesInfo$transcript_start[i] <= nsPos &&
        allGenesInfo$transcript_end[i] >= nsPos){
      pc <- c(pc, i)
      uniprotIDs <- c(uniprotIDs, allGenesInfo$uniprot_gn_id[i])
      pdbIDs <- c(pdbIDs, allGenesInfo$pdb[i])
    }
  }
  if (length(pc) == 0){
    stop("no experimental structure of encoding protein in Uniprot or PDB
          associate with your input nsPos, please re-enter another gene.")
  }

  # List all available pdb IDs in PDB, for user to choose one structure to observe
  # Check the choice and determine the pdb ID
  print("The available pdb structrues associated with your genes are:")
  print(unique(pdbIDs))
  if (choice == 1){
    input <- readline("please enter one ID to display \n")
    input <- toupper(input)
    while (!(input %in% pdbIDs)){
      print(unique(pdbIDs))
      input <- readline("your input is not in the displayed list,
                      please re-enter the choice.")
      input<-toupper(input)
    }
  } else {
    input <- pdbIDs[1]
  }

  # Find the residue position of nsPos
  for (i in 1:nrow(allGenesInfo)){
    if (allGenesInfo$pdb[i] == input){
      idx <- i
      break
    }
  }

  localPos <- nsPos - allGenesInfo$transcript_start[idx]
  #hsapiensSeq <- BSgenome.Hsapiens.UCSC.hg38
  chr <- paste('chr',as.character(chrName),sep="")
  if(allGenesInfo$strand[idx] == -1){
    strand <- "-"
  } else {
    strand <- "+"
  }
  oriGeneSeq <- getSeq(x = hsapiensSeq, names = chr,
                         start = allGenesInfo$transcript_start[idx],
                         end = allGenesInfo$transcript_end[idx],
                         strand = strand)

  pdbSeq <- getUniProt(allGenesInfo$uniprot_gn_id[idx])[[1]]
  if (nchar(oriGeneSeq) %% 3 != 0){
    oriGeneSeq <- oriGeneSeq[1:(nchar(oriGeneSeq)%/%3 * 3)]
  }
  oriAASeq <- as.character(translate(oriGeneSeq))
  nchar(oriAASeq)
  nchar(pdbSeq)
  if (grepl(substr(pdbSeq,1,5), oriAASeq, fixed = TRUE)){
    firstAALoc <- unlist(gregexpr(substr(pdbSeq,1,5), oriAASeq))[1]
    locStart <- allGenesInfo$transcript_start[idx] + firstAALoc * 3
    locEnd <- locStart + nchar(pdbSeq) * 3
  }

  if (nsPos <= locEnd && nsPos >= locStart){
    nssnpLocPos <- (nsPos - locStart) %/% 3 + 1
  } else {
    stop('The input nsSNP position is not included in the selected PDB structure.')
  }
  # Draw the pdb as cartoon, entire structure overall,
  # Point the residue that involves in nsSNP in red color
  #stc1 <-
  r3dmol() %>%
    m_add_model(data = m_fetch_pdb(input)) %>%
    m_set_style(style = m_style_cartoon(color = "cyan")) %>%
    m_add_style(sel = m_sel(resi = nssnpLocPos),
                style = m_style_cartoon(color = "red")) %>%
    m_zoom_to() %>%
    m_spin()
}
