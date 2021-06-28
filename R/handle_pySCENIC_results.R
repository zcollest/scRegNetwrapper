#' Uses reticulate package to handle and process the results from pySCENIC
#'
#' handle pySCENIC results
#' @param dir directory path to output pySCENIC results
#' @param output_loom_path path to output loom file
#' @param anndata_path path to initial anndata object
#' @param comparison_feature which obs feature from the anndata should we look at? for RSS.
#' @param regulon_path path to the regulon output from CisTarget
#' @param regulon_output_fname file name for writing processed regulons (.csv)
#' @param auc_mtx_fname filename for writing AUC matrix (.csv)
#' @param cellannot_fname filename for writing cellAnnot matrix (.csv)
#' @param RSS_fname filename for writing regulon specificity scores matrix (.csv)
#' @keywords pySCENIC
#' @export
#' @import reticulate
#' @examples
#' handle_pyscenic_results(dir, output_loom_path, anndata_path, comparison_feature="cell_type", regulon_path)

handle_pyscenic_results <- function(dir, output_loom_path, anndata_path, comparison_feature="cell_type", regulon_path) {
  handle_pyscenic_results_path <- system.file("python","handle_SCENIC_results.py", package="scRegNetwrapper")
  source_python(handle_pyscenic_results_path)
  results <- scenic_results_wrapper(dir, output_loom_path, anndata_path, comparison_feature, regulon_path)
  results_named <- list("regulons" = results[[1]], "AUC_mtx" = results[[2]],"cellAnnot" = results[[3]],"RSS"=results[[4]])
  return(results_named)
}
