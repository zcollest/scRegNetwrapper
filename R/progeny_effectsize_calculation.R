#' Cohen's D effect size calculation for progeny data (requires 2 group comparisons...e.g. healthy vs disease)
#'
#' This function takes progeny (or more generally pathway activity) data by cell and produces a dataframe of cohen's D scores for each pathway.
#' @param pathway_data_bycell dataframe of pathway data by cell
#' @keywords progeny effect size calculation
#' @export
#' @import Seurat
#' @import dorothea
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @import effsize
#' @examples
#' pathway_effsize_calc(pathway_data_bycell)

pathway_effsize_calc <- function(pathway_data_bycell){
  stats_df <- data.frame()
  for (pathway_name in unique(pathway_data_bycell$Pathway)){
    df <- pathway_data_bycell %>% filter(Pathway %in% pathway_name)
    binvar <- df$group
    scalevar <- as.numeric(df$Activity)
    test <- cohen.d(scalevar~binvar)
    stats <- data.frame(cbind(pathway_name, test[["estimate"]]))
    if (length(stats_df)<1){
      stats_df <- stats
    }
    else{
      stats_df <- rbind(stats_df,stats)
    }
  }
  names(stats_df) <- c("pathway","cohen_d")
  stats_df <- arrange(stats_df,desc(stats_df$cohen_d))
  return(stats_df)
}
