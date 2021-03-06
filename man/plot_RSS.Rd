% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_RSS.R
\name{plot_RSS}
\alias{plot_RSS}
\title{Uses reticulate package to plot RSS using python}
\usage{
plot_RSS(dir, cellAnnot, rss_cellType, title = "RSS_plot.png")
}
\arguments{
\item{cellAnnot}{cell annotations from pyscenic (one of the outputs of the \code{handle_pyscenic_results()} function)}

\item{rss_cellType}{regulon specificity scores by cell type / comparison feature (one of the outputs of the \code{handle_pyscenic_results()} function)}
}
\description{
Plot RSS
}
\examples{
plot_RSS(cellAnnot, rss_cellType)
}
\keyword{RSS}
\keyword{plot}
