#' Uses reticulate package to run pySCENIC from R
#'
#' Run pySCENIC pipeline
#' @param outdir directory path to output pySCENIC results
#' @param anndata_path path to anndata object
#' @param loom_path path to write initial loom file
#' @param tfs_path path to list of TFs
#' @param rank_db_path path to ranking database
#' @param motif_path path to motif list
#' @param output_loom_path path to write output loom file
#' @keywords pySCENIC
#' @export
#' @import reticulate
#' @examples
#' run_pyscenic(outdir, anndata_path, loom_path, tfs_path, rank_db_path, motif_path, output_loom_path)

run_pyscenic <- function(outdir, anndata_path, loom_path, tfs_path, rank_db_path, motif_path, output_loom_path) {
  pyscenic_wrapper(outdir, anndata_path, loom_path, tfs_path, rank_db_path, motif_path, output_loom_path)
}

