#' Launch Shiny App For Package nsSNPfinder
#'
#' A function that launches the shiny app for this package.
#' The app allows the user to preliminary check the information of nsSNP
#' locations within the gene of interest region.
#' The code has been placed in \code{./inst/shiny-scripts}.
#'
#' @return No return value but open up a shiny page.
#'
#' @examples
#' \dontrun{
#' nsSNPfinder::runNSSNPfinder()
#' }
#'
#' @references
#' RStudio Inc. (2013). Tabsets. Shiny Gallery. URL:https://shiny.rstudio.com/gallery/tabsets.html
#' Pagès, H. (2017). SNPlocs.Hsapiens.dbSNP144.GRCh38: SNP locations for Homo sapiens (dbSNPBuild 144). R package version 0.99.20.
#'
#' @export
#' @importFrom shiny runApp
runNSSNPfinder <- function() {
  appDir <- system.file("shiny-scripts",
                        package = "nsSNPfinder")
  shiny::runApp(appDir, display.mode = "normal")
  return()
}

# [END]
