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
#' run_pyscenic(dir, anndata_path, loom_path, tfs_path, rank_db_path, motif_path, output_loom_path)


run_pyscenic <- function(dir, anndata_path, loom_path, tfs_path, rank_db_path, motif_path, output_loom_path) {
  run_pyscenic_path <- system.file("python","run_pySCENIC.py", package="scRegNetwrapper")
  source_python(run_pyscenic_path)
  pyscenic_wrapper(dir, anndata_path, loom_path, tfs_path, rank_db_path, motif_path, output_loom_path)
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
#' @param adjacencies_fname filename of adjacencies output
#' @keywords GRN
#' @export
#' @import reticulate
#' @examples
#' run_grn(dir, loom_path, tfs_path,adjacencies_fname)

run_grn <- function(dir, loom_path, tfs_path, adjacencies_fname) {
  grn_wrapper(dir, loom_path, tfs_path, adjacencies_fname)
}


#' Uses reticulate package to run CisTarget (regulon prediction) from R
#'
#' Run pySCENIC pipeline
#' @param dir directory path to output pySCENIC results
#' @param adjacencies_fname filename of adjacencies output
#' @param loom_path path to write initial loom file
#' @param rank_db_path path to ranking database
#' @param motif_path path to motif list
#' @param regulons_fname filename of regulons output
#' @keywords CisTarget
#' @export
#' @import reticulate
#' @examples
#' run_cistarget(dir, adjacencies_fname, loom_path, rank_db_path, motif_path, regulons_fname)


run_cistarget <- function(dir, adjacencies_fname, loom_path, rank_db_path, motif_path, regulons_fname) {
  cistarget_wrapper(dir, adjacencies_fname, loom_path, rank_db_path, motif_path, regulons_fname)
}

#' Uses reticulate package to run AUCell from R
#'
#' Run AUCell from pySCENIC pipeline
#' @param dir directory path to output pySCENIC results
#' @param regulons_fname filename of regulons output
#' @param loom_path path to write initial loom file
#' @param output_loom_path path to write output loom file
#' @keywords pySCENIC
#' @export
#' @import reticulate
#' @examples
#' run_aucell(dir, regulons_fname, loom_path, output_loom_path)

run_aucell <- function(dir, regulons_fname, loom_path, output_loom_path) {
  aucell_wrapper(dir, regulons_fname, loom_path, output_loom_path)
}
