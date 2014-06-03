#' Run Poly Peak Parser
#' 
#' Runs the Poly Peak Parser shiny app in the system's default browser. Poly Peak Parser
#' is a web front end that reads, plots and parses double peaks from chromatogram files.
#' Instructions can be found on the webpage once it launches.
#' 
#' @return 
#' \code{\link{scf}} s4 object
#' 
#' @seealso
#' \code{\link{read.abif}}, \code{\link{readsangerseq}},  \code{\link{scf}}
#' 
#' @examples
#' PolyPeakParser()
#' 
#' @export

PolyPeakParser <- function() {
  runApp(system.file(package="sangerseqR", "PolyPeakParser"))
}