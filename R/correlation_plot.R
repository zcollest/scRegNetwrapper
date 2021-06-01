correlation_analysis <- function(tf_data, pathway_data){
  names(pathway_data) <- paste(names(pathway_data),"(P)", sep="")
  joint_df <- cbind(tf_data, pathway_data)
  res <- cor(joint_df, method="pearson")
  cols <- brewer.pal(8, "RdBu")
  corrplot(res, type = "upper", order = "hclust",
           tl.col = "black", tl.srt = 45, tl.cex=0.6, col=colorRampPalette(c("navy","grey","dark red"))(200))
}
