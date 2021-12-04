utils::globalVariables(c("BSgenome.Hsapiens.UCSC.hg38", "SNPlocs.Hsapiens.dbSNP144.GRCh38"))
# Set up global variables to extract information of transcripts for all functions
hsapiens <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl",host="https://www.ensembl.org")
hsapiensSeq <- BSgenome.Hsapiens.UCSC.hg38

# Set up global variable for saving run-time of examples
# for chromosome 1, coordinate (2321253, 2391707)
demoAllGenesInfo <-biomaRt::getBM(attributes = c("hgnc_symbol","ensembl_gene_id",
                                             "transcript_start","transcript_end",
                                             "transcript_biotype","strand"),
                              filters = c('chromosome_name','start','end'),
                              values = list(3, 49395401, 49395500),
                              mart = hsapiens)
