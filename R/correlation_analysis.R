#' Generating a correlation matrix of pathways and TFs
#'
#' Takes cell-wise data from the 'scores_handling' functions and generates a correlation matrix
#' @param tf_data transcription factor data by cell
#' @param pathway_data pathway data by cell
#' @param render_plot render plot? default = TRUE
#' @param return_corr_data return the correlation data? default = FALSE
#' @keywords correlation plot
#' @export
#' @import Seurat
#' @import dplyr
#' @import tidyr
#' @import pheatmap
#' @import corrplot
#' @examples
#' correlation_analysis(tf_data_bycell, pathway_data_bycell)

correlation_analysis <- function(tf_data, pathway_data,render_plot=TRUE,return_corr_data=FALSE){
  names(pathway_data) <- paste(names(pathway_data),"(P)", sep="")
  joint_df <- cbind(tf_data, pathway_data)
  res <- cor(joint_df, method="pearson")
  cols <- brewer.pal(8, "RdBu")
  corrplot(res, type = "upper", order = "hclust",
           tl.col = "black", tl.srt = 45, tl.cex=0.6, col=colorRampPalette(c("navy","grey","dark red"))(200))
}
