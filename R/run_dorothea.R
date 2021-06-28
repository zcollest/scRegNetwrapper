#' Quantifying TF activity from DoRothEA Regulons
#'
#' Takes a Seurat object and a vector of DoRothEA confidence levels and returns the same Seurat object with a "dorothea" assay of TF activity scores.
#' @param seurat_obj Seurat object
#' @param conf_scores vector of DoRothEA confidence scores (A,B,C,D,E)
#' @param cores number of cores to use for calculation (default = 16)
#' @keywords TF Activity with DoRothEA Regulons
#' @export
#' @import Seurat
#' @import dorothea
#' @import dplyr
#' @import tibble
#' @examples
#' run_dorothea(pbmc10k, c("A","B","C","D"),cores=20)


run_dorothea <- function(seurat_obj, conf_scores, cores=16){

  ## read Dorothea regulons for Human:
  dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))

  ## obtain the regulons based on interactions with confidence levels in conf_scores
  regulon <- dorothea_regulon_human %>%
    dplyr::filter(confidence %in% conf_scores)

  ## compute Viper Scores
  DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj <- run_viper(seurat_obj, regulon,
                          options = list(method = "scale", minsize = 4,
                                         eset.filter = FALSE, cores = cores,
                                         verbose = FALSE))
  ## scale the tf activity scores
  seurat_obj <- Seurat::ScaleData(seurat_obj, assay = "dorothea")
  return(seurat_obj)
}


#' Returning DoRothEA Regulons
#'
#' Returns DoRothEA regulons from a vector of DoRothEA confidence levels.
#' @param conf_scores vector of DoRothEA confidence scores (A,B,C,D,E)
#' @keywords DoRothEA Regulon
#' @export
#' @examples
#' get_dorothea_regulons(c("A","B","C","D"))


get_dorothea_regulons <- function(conf_scores){
  ## read Dorothea regulons for Human:
  dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))
  ## obtain the regulons based on interactions with confidence levels in conf_scores
  regulon <- dorothea_regulon_human %>%
    dplyr::filter(confidence %in% conf_scores)
  return(regulon)
}


#' Quantifying TF activity using a data frame of custom regulons and the VIPER algorithm
#'
#' This function uses VIPER to quantify TF activity given a data frame of custom regulons and returns a Seurat object assay with the results.
#' @param seurat_obj seurat object
#' @param regulons data frame of custom regulons (must be same format as DoRothEA regulons)
#' @param assay_name string to name the resulting Seurat object assay
#' @param cores number of cores to use for calculation (default is 16)
#' @keywords Quantifying TF activity with custom regulons and the VIPER algorithm
#' @export
#' @examples
#' custom_regulons_calc(data,regulons,"custom_TFactivity"))

custom_regulons_calc <- function(seurat_obj, regulons, assay_name, cores=16){
  ## compute Viper Scores
  if ("dorothea" %in% Assays(seurat_obj)){
    dorothea_data <- as.matrix(Seurat::GetAssayData(seurat_obj, assay = "dorothea",
                                                           slot = "data"))
    dorothea_data_scaled <- as.matrix(Seurat::GetAssayData(seurat_obj, assay = "dorothea",
                                                    slot = "scale.data"))
    seurat_obj[["dorothea"]] <- NULL
  }
  DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj <- run_viper(seurat_obj, regulons,
                          options = list(method = "scale", minsize = 4,
                                         eset.filter = FALSE, cores=cores,
                                         verbose = FALSE))
  seurat_obj[[assay_name]] <- seurat_obj[["dorothea"]]
  seurat_obj@assays[["dorothea"]]@data <- dorothea_data
  seurat_obj@assays[["dorothea"]]@scale.data <- dorothea_data_scaled

  ## Scale the tf activity scores
  seurat_obj <- Seurat::ScaleData(seurat_obj, assay = assay_name)
  return(seurat_obj)
}


