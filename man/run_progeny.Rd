% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_progeny.R
\name{run_progeny}
\alias{run_progeny}
\title{Calculating pathway activity scores from PROGENy pathways}
\usage{
run_progeny(seurat_obj, num_genes, organism = "Human")
}
\arguments{
\item{seurat_obj}{seurat object}

\item{num_genes}{number of pathway genes to use for calculations}

\item{organism}{"Human" or "Mouse", defaults to "Human"}
}
\description{
Takes a seurat object and a number of pathway genes and returns the same Seurat object with a "progeny" assay.
}
\examples{
run_progeny(pbmc10k, 500, "Human")
}
\keyword{PROGENy}
\keyword{calculations}
