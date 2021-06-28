#' Uses reticulate package to plot RSS using python
#'
#' Plot RSS
#' @param cellAnnot cell annotations from pyscenic (one of the outputs of the `handle_pyscenic_results()` function)
#' @param rss_cellType regulon specificity scores by cell type / comparison feature (one of the outputs of the `handle_pyscenic_results()` function)
#' @keywords plot RSS
#' @export
#' @import reticulate
#' @examples
#' plot_RSS(cellAnnot, rss_cellType)

plot_RSS <- function(dir, cellAnnot, rss_cellType, title="RSS_plot.png") {
  plot_RSS_path <- system.file("python","plot_RSS.py", package="scRegNetwrapper")
  source_python(plot_RSS_path)
  plot_RSS_wrapper(dir, cellAnnot, rss_cellType, title)
}
