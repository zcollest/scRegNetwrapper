#' Handling the dorothea scores (proportion adjustments, summarizing by cell type, etc.)
#'
#' Takes a Seurat object and produces a list of data frames that are necessary for downstream analysis
#' @param seurat_obj Seurat object
#' @param comparison_feature which metadata do you want to do comparisons for?
#' @param topTFs number of top TFs to be us ed for heatmap downstream
#' @param assay_name if using custom TFs, what's the assay name?
#' @keywords dorothea score handling
#' @export
#' @import Seurat
#' @import dorothea
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @examples
#' handle_dorothea_scores(seurat_obj,seurat_obj$meta.data$cell_type, 30)

handle_dorothea_scores <- function(seurat_obj,comparison_feature, topTFs, assay_name=NULL){
  # extract scores
  if (missing(assay_name)){
    viper_scores_df <- GetAssayData(seurat_obj, slot = "scale.data",
                                    assay = "dorothea") %>%
      data.frame(check.names = F) %>%
      t()
  }
  else{
    viper_scores_df <- GetAssayData(seurat_obj, slot = "scale.data",
                                    assay = assay_name) %>%
      data.frame(check.names = F) %>%
      t()
  }

  # get a dataframe of cell cell types, frequency, and cell type proportions
  seurat_obj$cell_type <- comparison_feature
  cell_type_table <- data.frame(table(seurat_obj$cell_type))
  totalcells <- dim(seurat_obj)[2]
  # loop to add cell proportions to df and duplicate rows according to cell type counts
  i <- 1
  for (freq in cell_type_table$Freq){
    cell_type_table$prop[i] <- freq/totalcells
    cell_type_table <- rbind(cell_type_table, cell_type_table[rep(i, freq-1), ])
    i <- i + 1
  }
  cell_type_table <- cell_type_table[order(cell_type_table$Var1),]

  # adjust the viper scores by cell proportion
  cell_types <- data.frame(seurat_obj@meta.data[["cell_type"]])
  viper_scores_df_prop <- cbind(cell_types, viper_scores_df)
  viper_scores_df_prop <- viper_scores_df_prop[order(viper_scores_df_prop[,1]),]
  viper_scores_df_prop$prop <- cell_type_table$prop
  viper_scores_df_prop[,1] <- NULL

  # actually calculate a df with correct cell type proportions
  viper_scores_df_proportions <- viper_scores_df_prop[,1:ncol(viper_scores_df_prop)-1] * viper_scores_df_prop[,ncol(viper_scores_df_prop)]

  ## create a dataframe containing the cells and their clusters
  CellsClusters <- data.frame(cell = names(Idents(seurat_obj)),
                              group = seurat_obj@meta.data[["cell_type"]],
                              check.names = F)

  # create a dataframe with the Viper score per cell and its clusters (both normal and adjusted for cell proportions)
  viper_scores_clusters <- viper_scores_df  %>%
    data.frame() %>%
    rownames_to_column("cell") %>%
    gather(tf, activity, -cell) %>%
    inner_join(CellsClusters)

  viper_scores_clusters_proportions <- viper_scores_df_proportions  %>%
    data.frame() %>%
    rownames_to_column("cell") %>%
    gather(tf, activity, -cell) %>%
    inner_join(CellsClusters)


  # summarize the Viper scores by cell population (both normal and adjusted for cell proportions)
  summarized_viper_scores <- viper_scores_clusters %>%
    group_by(tf, group) %>%
    summarise(avg = mean(activity),
              std = sd(activity))

  summarized_viper_scores_proportions <- viper_scores_clusters_proportions %>%
    group_by(tf, group) %>%
    summarise(avg = mean(activity),
              std = sd(activity))

  # calculate pooled SD
  pooledSD_bytf <- summarized_viper_scores_proportions %>%
    group_by(tf) %>%
    summarise(poolSD = pooledSD(std))
  total_celltypes <- length(unique(seurat_obj$cell_type))
  pooledSD_bytf <- pooledSD_bytf[rep(seq_len(nrow(pooledSD_bytf)), each = total_celltypes), ]
  summarized_viper_scores_proportions <- cbind(pooledSD_bytf,summarized_viper_scores_proportions$group,summarized_viper_scores_proportions$avg)
  names(summarized_viper_scores_proportions) <- c("tf","poolSD","group","avg")
  summarized_viper_scores_proportions$adj_avg <- summarized_viper_scores_proportions$avg/summarized_viper_scores_proportions$poolSD

  # get the most variable TFs
  highly_variable_tfs_adj <- summarized_viper_scores_proportions %>%
    group_by(tf) %>%
    mutate(var = var(avg))  %>%
    ungroup() %>%
    top_n(topTFs*total_celltypes, var) %>%
    distinct(tf)

  highly_variable_tfs <- summarized_viper_scores %>%
    group_by(tf) %>%
    mutate(var = var(avg))  %>%
    ungroup() %>%
    top_n(topTFs*total_celltypes, var) %>%
    distinct(tf)

  topTFs_adjavg <- summarized_viper_scores_proportions %>%
    semi_join(highly_variable_tfs_adj, by = "tf") %>%
    dplyr::select(-c(poolSD,avg)) %>%
    spread(tf, adj_avg) %>%
    data.frame(row.names = 1, check.names = FALSE)

  topTFs_avg <- summarized_viper_scores %>%
    semi_join(highly_variable_tfs, by = "tf") %>%
    dplyr::select(-std) %>%
    spread(tf, avg) %>%
    data.frame(row.names = 1, check.names = FALSE)


  # return everything
  handled_scores <- list("scores_bycell" = data.frame(viper_scores_clusters),"proportionadjusted_scores_bycell" = data.frame(viper_scores_clusters_proportions),
                    "scores_bygroup" = data.frame(summarized_viper_scores), "proportionadjusted_scores_bygroup" = data.frame(summarized_viper_scores_proportions),
                    "topTFs" = topTFs_avg, "proportionadjusted_topTFs" = topTFs_adjavg)
  return(handled_scores)
}
