#' Wrapper function to run dorothea and progeny and get basic downstream results
#'
#' Takes a seurat object an
#' @param seurat_obj seurat object
#' @param dorothea_conf_scores vector of DoRothEA confidence scores (default = A->D)
#' @param topTFs topTFs for downstream heatmap (default = 30)
#' @param cores number of cores to use for calculation (default = 16)
#' @param top_progeny_genes number of genes per progeny pathway to use for calculations (default = 500)
#' @param organism for progeny, use Human or Mouse pathways (default = "Human")
#' @param comparison_feature which metadata do you want to do comparisons for?
#' @keywords all-in-one wrapper function for dorothea and progeny
#' @export
#' @import Seurat
#' @import dorothea
#' @import progeny
#' @import dplyr
#' @import tibble
#' @import effsize
#' @examples
#' scRegNetwrapper(seurat_obj = pbmc10k, dorothea_conf_scores = c("A","B","C","D"), topTFs = 30, cores = 20, top_progeny_genes = 100, organism = "Human", comparison_feature = Idents(pbmc10k))

scRegNetwrapper <- function(seurat_obj,dorothea_conf_scores = c("A","B","C","D"), topTFs = 30, cores=16, top_progeny_genes=500, organism="Human", comparison_feature){
  seurat_obj <- run_dorothea(seurat_obj, dorothea_conf_scores,cores)
  seurat_obj <- run_progeny(seurat_obj, top_progeny_genes, organism)
  dorothea_data <- handle_dorothea_scores(seurat_obj, comparison_feature, topTFs)
  progeny_data <- handle_progeny_scores(seurat_obj,comparison_feature)

  if (length(levels(comparison_feature))==2){
    progeny_effectsize <- pathway_effsize_calc(progeny_data$proportionadjusted_scores_bycell)
    dorothea_effectsize <- tf_effsize_calc(dorothea_data$proportionadjusted_scores_bycell)
    obj_name <- seurat_obj@project.name
    results <- list(seurat_obj = seurat_obj, "dorothea_data" = dorothea_data, "progeny_data" = progeny_data,
                    "progeny_cohenD" = progeny_effectsize, "dorothea_cohenD" = dorothea_effectsize)
    return(results)
  }
  else{
    results <- list(seurat_obj = seurat_obj, "dorothea_data" = dorothea_data, "progeny_data" = progeny_data)
    return(results)
  }
}

