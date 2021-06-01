effsize_test <- function(pathway_data_bycell){
  stats_df <- data.frame()
  for (TFname in unique(tf_activity_data$tf)){
    df <- tf_activity_data %>% filter(tf %in% TFname)
    binvar <- df$group
    scalevar <- as.numeric(df$activity)
    test <- cohen.d(scalevar~binvar)
    stats <- data.frame(cbind(TFname, abs(test[["estimate"]])))
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