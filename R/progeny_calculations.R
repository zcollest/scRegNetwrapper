#' Calculating pathway activity scores from PROGENy pathways
#'
#' This function takes a seurat object and a number of pathway genes and returns the same Seurat object with a "progeny" assay.
#' @param seurat_obj seurat object
#' @param num_genes number of pathway genes to use for calculations
#' @param organism "Human" or "Mouse", defaults to "Human"
#' @keywords PROGENy calculations
#' @export
#' @import Seurat
#' @import progeny
#' @import dplyr
#' @import tibble
#' @examples
#' progeny_score_calc(seurat_obj, 500, "Human")

progeny_score_calc <- function(seurat_obj, num_genes, organism="Human"){
  seurat_obj <- progeny(seurat_obj, scale=FALSE, organism=organism, top=num_genes, perm=1,
                        return_assay = TRUE)

  # scale the pathway activity scores
  seurat_obj <- Seurat::ScaleData(seurat_obj, assay = "progeny")
}

#' Returning PROGENy pathways
#'
#' This function returns PROGENy pathways
#' @param full_model return the full model? Defaults to TRUE.
#' @param num_genes number of top genes to return for each pathway if full_model = FALSE.
#' @param organism "Human" or "Mouse", defaults to "Human".
#' @keywords PROGENy pathways
#' @export
#' @examples
#' get_progeny_pathways(full_model=FALSE, num_genes=100, organism = "Human")

get_progeny_pathways <- function(full_model = TRUE, num_genes = NULL, organism = "Human"){
  if (full_model == TRUE){
    return(get("model_human_full", envir = .GlobalEnv))
  }
  else {
    return(getModel(organism = organism, top = num_genes))
  }
}



#' Calculate pathway activity using a dataframe of custom pathways
#'
#' This function uses the PROGENy framework to quantify pathway activity given a dataframe of custom pathways and returns a seurat object assay with the results.
#' @param seurat_obj seurat object
#' @param pathways dataframe of custom pathways (must be same format as PROGENy pathways)
#' @param num_genes number of top genes to return for each pathway if full_model = FALSE.
#' @param organism "Human" or "Mouse", defaults to "Human".
#' @param assay_name string to name the resulting Seurat object assay
#' @keywords custom pathway activity calculation
#' @export
#' @examples
#' custom_pathway_calc(seurat_obj, pathways, 500, "Human", "custom_pathway_activity"))

custom_pathway_calc <- function(seurat_obj, pathways, num_genes, organism="Human", assay_name){
  # check if progeny assay exists before rewriting to that assay - if exists, save that data separately
  if ("progeny" %in% Assays(seurat_obj)){
    progeny_data <- as.matrix(Seurat::GetAssayData(seurat_obj, assay = "progeny", slot = "data"))
    progeny_data_scaled <- as.matrix(Seurat::GetAssayData(seurat_obj, assay = "progeny", slot = "scale.data"))
    seurat_obj[["progeny"]] <- NULL
  }
  DefaultAssay(seurat_obj) <- "RNA"
  if (organism == "Human"){
    model_human_full <- pathways
    seurat_obj <- progeny_score_calc(seurat_obj,num_genes,"Human")
  }
  else if (organism == "Mouse"){
    model_mouse_full <- pathways
    seurat_obj <- progeny_score_calc(seurat_obj,num_genes,"Mouse")
  }
  seurat_obj[[assay_name]] <- seurat_obj[["progeny"]]
  seurat_obj@assays[["progeny"]]@data <- progeny_data
  seurat_obj@assays[["progeny"]]@scale.data <- progeny_data_scaled

  ## Scale the pathway activity scores
  seurat_obj <- Seurat::ScaleData(seurat_obj, assay = assay_name)
  return(seurat_obj)
  # reload the normal model_human_full and model_mouse_full
  detach("package:progeny", unload=TRUE)
  library(progeny)
}




