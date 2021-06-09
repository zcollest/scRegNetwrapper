#' Given a list of genes, find which pathways they are associated with and see how statistically significant those pathways are for a given comparison (e.g. healthy vs disease)
#'
#' This function takes a vector of genes and pathway activity data (by group, not by cell...optional) and a data frame of cohen's D scores (if the comparison is 2 groups...optional).
#' If only the gene vector is supplied, the function will return all pathways for which a given gene is part of the pathway's gene set. 
#' @param gene_vector vector of genes
#' @param pathway_data_bygroup pathway data summarized by comparison group
#' @param effect_size_data dataframe of cohen's D scores for a given comparison (optional, only if comparison is 2 groups)
#' @param custom_pathways searches custom pathways for genes in gene vector, default is FALSE.
#' @param pathway_data_bygroup pathway data summarized by comparison group
#' @param effect_size_data dataframe of cohen's D scores for a given comparison (optional, only if comparison is 2 groups)
#' @param custom_pathways searches custom pathways for genes in gene vector, default is FALSE.
#' @keywords find associated TFs
#' @export
#' @import Seurat
#' @import progeny
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @examples
#' find_associated_pathways(gene_vector, pathway_data_bygroup, effect_size_data, custom_pathways=custom_pathways, progeny_pathways=FALSE)

find_associated_pathways <- function(gene_vector, pathway_data_bygroup=NULL, effect_size_data=NULL, custom_pathways=NULL, progeny_pathways=TRUE, num_genes=500, organism="Human"){
  pathway_associations <- list()
  if (missing(custom_pathways)){
    pathways <- get_progeny_pathways(num_genes = num_genes, organism = organism)
  }
  else {
   pathways  <- custom_pathways
  }
  if (missing(pathway_data_bygroup)){
    for (gene_name in gene_vector) {
      pathway_associations[[gene_name]] <- pathways %>% filter(gene %in% gene_name)
    }
  }
  else{
    rankings_by_var <- pathway_data_bygroup %>% group_by(pathway) %>% summarise(var = var(avg))
    rankings_by_var <- rankings_by_var[order(-rankings_by_var$var),]
    rankings_by_var$rank <- 1:nrow(rankings_by_var)
    for (gene_name in gene_vector) {
      pathway_associations[[gene_name]] <- pathways %>% filter(gene %in% gene_name)
      path <- rankings_by_var %>% filter(pathway %in% pathway_associations[[gene_name]]$pathway)
      pathway_associations[[gene_name]] <- merge(pathway_associations[[gene_name]], path, by="pathway")
      pathway_associations[[gene_name]] <- pathway_associations[[gene_name]][order(pathway_associations[[gene_name]]$rank),]
      row.names(pathway_associations[[gene_name]]) <- NULL
      if (missing(effect_size_data)){
        next
      }
      else{
        pathway_associations[[gene_name]] <- merge(pathway_associations[[gene_name]], effect_size_data, by="pathway")
        pathway_associations[[gene_name]] <- pathway_associations[[gene_name]][order(pathway_associations[[gene_name]]$rank),]
        pathway_associations[[gene_name]]$var <- NULL
        row.names(pathway_associations[[gene_name]]) <- NULL
      }
    }
  }
  return(pathway_associations)
}
