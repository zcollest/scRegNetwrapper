#' Plotting dorothea/progeny scores by comparison group on a heatmap
#'
#' This function takes data from the 'scores_handling' functions and plots a heatmap of the top TFs or pathways
#' @param data Data comes from the scores_handling functions, plot the top TFs or pathways
#' @param title title for the heatmap
#' @keywords downstream heatmap
#' @export
#' @import Seurat
#' @import dplyr
#' @import tidyr
#' @import pheatmap
#' @examples
#' downstream_heatmap(data, "top 30 TFs, by indication")

downstream_heatmap <- function(df, title){
  paletteLength = 100
  myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)

  Breaks = c(seq(min(df), 0,
                        length.out=ceiling(paletteLength/2) + 1),
                    seq(max(df)/paletteLength,
                        max(df),
                        length.out=floor(paletteLength/2)))
  hmap = pheatmap(t(df[,-1]),fontsize=14,
                          fontsize_row = 10,
                          color=myColor, breaks = Breaks,
                          main = title, angle_col = 45,
                          treeheight_col = 0,  border_color = NA)
}
