#' Given a list of genes, find which TFs they are associated with and see how statistically significant those TFs are for a given comparison (e.g. healthy vs disease)
#'
#' This function takes a vector of genes and TF activity data (by group, not by cell...optional) and a dataframe of cohen's D scores (if the comparison is 2 groups...optional).
#' If only the gene vector is supplied, the function will return all TFs for which a given gene is part of the regulon for the entire dorothea collection of TFs.
#' @param gene_vector vector of genes
#' @param tf_data_bygroup TF data summarized by comparison group
#' @param effect_size_data dataframe of cohen's D scores for a given comparison (optional, only if comparison is 2 groups)
#' @param custom_regulons searches custom TF regulons for genes in gene vector, default is FALSE.
#' @keywords find associated TFs
#' @export
#' @import Seurat
#' @import dorothea
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @examples
#' find_associated_TFs(gene_vector, tf_data_bygroup, effect_size_data, custom_regulons=custom_regulons)

find_associated_TFs <- function(gene_vector, tf_data_bygroup=NULL, effect_size_data=NULL, custom_regulons=NULL){
  associated_TFs <- list()
  if (missing(custom_regulons)){
    dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))
    regulon <- dorothea_regulon_human %>% filter(confidence %in% c("A","B","C","D","E"))
  }
  else {
    regulon <- custom_regulons
  }
  if (missing(tf_data_bygroup)){
    for (gene in gene_vector) {
      associated_TFs[[gene]] <- regulon %>% filter(target %in% gene)
    }
  }
  else{
    rankings_by_var <- tf_data_bygroup %>% group_by(tf) %>% summarise(var = var(avg))
    rankings_by_var <- rankings_by_var[order(-rankings_by_var$var),]
    rankings_by_var$rank <- 1:nrow(rankings_by_var)
    for (gene in gene_vector) {
      associated_TFs[[gene]] <- regulon %>% filter(target %in% gene)
      tfs <- rankings_by_var %>% filter(tf %in% associated_TFs[[gene]]$tf)
      associated_TFs[[gene]] <- merge(associated_TFs[[gene]], tfs, by="tf")
      associated_TFs[[gene]] <- associated_TFs[[gene]][order(associated_TFs[[gene]]$rank),]
      row.names(associated_TFs[[gene]]) <- NULL
      if (missing(effect_size_data)){
        next
      }
      else{
        associated_TFs[[gene]] <- merge(associated_TFs[[gene]], effect_size_data, by="tf")
        associated_TFs[[gene]] <- associated_TFs[[gene]][order(associated_TFs[[gene]]$rank),]
        associated_TFs[[gene]]$var <- NULL
        row.names(associated_TFs[[gene]]) <- NULL
      }
    }
  }
  return(associated_TFs)
}
