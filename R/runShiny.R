#' @title MoPAC pipeline
#'
#' @description Interactive use of the MoPAC pipeline.
#'
#' @author Oscar D Villarreal, \email{oscardvillarreal@gmail.com}
#' @keywords mopac pipeline shiny gui
#'
#' @usage run.MoPAC()
#'
#' @export
#'

run.MoPAC <- function() {
  appDir <- system.file("shiny", "pipeline", package = "MoPAC")
  shiny::runApp(appDir, display.mode = "normal")
}
