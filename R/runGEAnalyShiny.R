#' Launch Shiny App for GEAnaly
#'
#' A function that launches the Shiny app for GEAnaly.
#' The purpose of this app is providing a shiny app page with
#' user-friendly UI that assists users without programming
#' experience to perform the analyses easily and
#' efficiently. The app includes the UI for all the analyses performed in
#' GEAnaly package. The code has been placed in \code{./inst/shiny-scripts}.
#'
#' @return No return value but open up a Shiny page.
#'
#' @examples
#' \dontrun{
#'
#' GEAnaly::runGEAnalyShiny()
#' }
#'
#' @references
#' Grolemund, G. (2015). Learn Shiny - Video Tutorials.
#' \href{https://shiny.rstudio.com/tutorial/}{Link}
#'
#' @export
#' @importFrom shiny runApp

runGEAnalyShiny <- function() {
  appDir <- system.file("shiny-scripts",
                        package = "GEAnaly")
  shiny::runApp(appDir, display.mode = "normal")
  return(invisible(NULL))
}

# [END]
