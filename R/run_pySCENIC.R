#' Uses reticulate package to run pySCENIC from R
#'
#' Run pySCENIC pipeline
#' @param dir directory path to output pySCENIC results
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


run_pyscenic <- function(dir, anndata_path, loom_path, tfs_path, rank_db_path, motif_path, output_loom_path) {
  pyscenic_wrapper(outdir, anndata_path, loom_path, tfs_path, rank_db_path, motif_path, output_loom_path)
}


#' Uses reticulate package to set up loom file for pySCENIC pipeline
#'
#' Run GRN from pySCENIC
#' @param anndata_path path to anndata object
#' @param loom_path path to write initial loom file
#' @keywords pySCENIC loom set up
#' @export
#' @import reticulate
#' @examples
#' run_loom_setup(anndata_path, loom_path)


run_loom_setup <- function(anndata_path, loom_path){
  setup_loom_wrapper(anndata_path, loom_path)
}

#' Uses reticulate package to run GRN from R
#'
#' Run GRN from pySCENIC pipeline
#' @param dir directory path to output pySCENIC results
#' @param loom_path path to write initial loom file
#' @param tfs_path path to list of TFs
#' @keywords GRN
#' @export
#' @import reticulate
#' @examples
#' run_grn(dir, loom_path, tfs_path,)

run_grn <- function(dir, loom_path, tfs_path) {
  grn_wrapper(dir, loom_path, tfs_path)
}


#' Uses reticulate package to run CisTarget (regulon prediction) from R
#'
#' Run pySCENIC pipeline
#' @param dir directory path to output pySCENIC results
#' @param loom_path path to write initial loom file
#' @param rank_db_path path to ranking database
#' @param motif_path path to motif list
#' @keywords CisTarget
#' @export
#' @import reticulate
#' @examples
#' run_cistarget(dir, loom_path, rank_db_path, motif_path)


run_cistarget <- function(dir, loom_path, rank_db_path, motif_path) {
  cistarget_wrapper(dir, loom_path, rank_db_path, motif_path)
}

#' Uses reticulate package to run AUCell from R
#'
#' Run AUCell from pySCENIC pipeline
#' @param dir directory path to output pySCENIC results
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
#' run_aucell(dir, loom_path, output_loom_path)

run_aucell <- function(dir, loom_path, output_loom_path) {
  aucell_wrapper(dir, loom_path, output_loom_path)
}
