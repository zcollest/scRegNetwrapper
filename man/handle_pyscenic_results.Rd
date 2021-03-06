% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/handle_pySCENIC_results.R
\name{handle_pyscenic_results}
\alias{handle_pyscenic_results}
\title{Uses reticulate package to handle and process the results from pySCENIC}
\usage{
handle_pyscenic_results(
  dir,
  output_loom_path,
  anndata_path,
  comparison_feature = "cell_type",
  regulon_path
)
}
\arguments{
\item{dir}{directory path to output pySCENIC results}

\item{output_loom_path}{path to output loom file}

\item{anndata_path}{path to initial anndata object}

\item{comparison_feature}{which obs feature from the anndata should we look at? for RSS.}

\item{regulon_path}{path to the regulon output from CisTarget}

\item{regulon_output_fname}{file name for writing processed regulons (.csv)}

\item{auc_mtx_fname}{filename for writing AUC matrix (.csv)}

\item{cellannot_fname}{filename for writing cellAnnot matrix (.csv)}

\item{RSS_fname}{filename for writing regulon specificity scores matrix (.csv)}
}
\description{
handle pySCENIC results
}
\examples{
handle_pyscenic_results(dir, output_loom_path, anndata_path, comparison_feature="cell_type", regulon_path)
}
\keyword{pySCENIC}
