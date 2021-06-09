#' Cohen's D effect size calculation for dorothea data (requires 2 group comparisons...e.g. healthy vs disease)
#'
#' This function takes a dorothea (or more generally tf activity) data by cell and produces a dataframe of cohen's D scores for each TF.
#' @param tf_data_bycell dataframe of TF data by cell
#' @keywords dorothea effect size calculation
#' @export
#' @import Seurat
#' @import dorothea
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @import effsize
#' @examples
#' tf_effsize_calc(tf_data_bycell=tf_data_bycell)

tf_effsize_calc <- function(tf_data_bycell){
  stats_df <- data.frame()
  for (TFname in unique(tf_data_bycell$tf)){
    df <- tf_data_bycell %>% filter(tf %in% TFname)
    binvar <- df$group
    scalevar <- as.numeric(df$activity)
    test <- cohen.d(scalevar~binvar)
    stats <- data.frame(cbind(TFname, test[["estimate"]]))
    if (length(stats_df)<1){
      stats_df <- stats
    }
    else{
      stats_df <- rbind(stats_df,stats)
    }
  }
  names(stats_df) <- c("tf","cohen_d")
  stats_df <- arrange(stats_df,desc(stats_df$cohen_d))
  return(stats_df)
}
