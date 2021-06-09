#' Reformat the progeny pathway data when you select only the top n genes
#'
#' When selecting progeny pathway data and choosing only a subset of the top genes from each pathway, the output is a sparse matrix and is not conducive to downstream analysis. This reformats to a long data frame. 
#' @param pathway_data progeny pathway data for top n genes (from `get_progeny_pathways()` function)
#' @keywords reformat progeny pathways
#' @export
#' @examples
#' reformat_progeny_pathways(pathway_data = top500genes_progenypathways)

reformat_progeny_pathways <- function(pathway_data){
  index <- data.frame(which(pathway_data!=0, arr.ind=TRUE))
  for (i in 1:nrow(index)){
    weight_index <- as.numeric(index[i,])
    weight <- pathway_data[weight_index[1],weight_index[2]]
    index$weight[i] <- weight
    index$pathway[i] <- colnames(pathway_data)[index$col[i]]
    }
  index$gene <- rownames(index)
  pathways_top_genes <- index[,c(5,4,3)]
  rownames(pathways_top_genes) <- NULL
  return(pathways_top_genes)
}