#' Handling the progeny scores (proportion adjustments, summarizing by cell type, etc.)
#'
#' Takes a Seurat object and produces a list of data frames that are necessary for downstream analysis
#' @param seurat_obj seurat object
#' @param comparison_feature which metadata feature do you want to do comparisons for?
#' @param assay_name if using custom TFs, what's the assay name?
#' @keywords progeny score handling
#' @export
#' @import Seurat
#' @import progeny
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @examples
#' handle_progeny_scores(pnmc10k,30)

handle_progeny_scores <- function(seurat_obj, comparison_feature, assay_name = NULL){
  if (missing(assay_name)){
    progeny_scores_df <- as.data.frame(t(GetAssayData(seurat_obj, slot = "scale.data",
                                                      assay = "progeny"))) %>% data.frame(check.names = F)
  }
  else{
    progeny_scores_df <- as.data.frame(t(GetAssayData(seurat_obj, slot = "scale.data",
                                                      assay = assay_name))) %>% data.frame(check.names = F)
  }
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


  ## We adjust the viper scores by cell proportion
  cell_types <- data.frame(seurat_obj@meta.data[["cell_type"]])
  progeny_scores_df_prop <- cbind(cell_types, progeny_scores_df)
  progeny_scores_df_prop <- progeny_scores_df_prop[order(progeny_scores_df_prop[,1]),]
  progeny_scores_df_prop$prop <- cell_type_table$prop
  progeny_scores_df_prop[,1] <- NULL
  # actually calculate a df with correct cell type proportions
  progeny_scores_df_proportions <- progeny_scores_df_prop[,1:ncol(progeny_scores_df_prop)-1] * progeny_scores_df_prop[,ncol(progeny_scores_df_prop)]


  ##create a data frame with the specification of the cells that belong toeach cluster to match with the Progeny scores
  CellsClusters <- data.frame(Cell = names(Idents(seurat_obj)),
                              group = as.character(seurat_obj@meta.data[["cell_type"]]),
                              stringsAsFactors = FALSE)

  ## We create a data frame with the Viper score per cell and its clusters (both normal and adjusted for cell proportions)
  progeny_scores_clusters <- progeny_scores_df  %>%
    data.frame() %>%
    rownames_to_column("Cell") %>%
    gather(Pathway, Activity, -Cell) %>%
    inner_join(CellsClusters)

  progeny_scores_clusters_proportions <- progeny_scores_df_proportions  %>%
    data.frame() %>%
    rownames_to_column("Cell") %>%
    gather(Pathway, Activity, -Cell) %>%
    inner_join(CellsClusters)

  ## We summarize the progeny scores by cell population (both normal and adjusted for cell proportions)
  summarized_progeny_scores <- progeny_scores_clusters %>%
    group_by(Pathway, group) %>%
    summarise(avg = mean(Activity),
              std = sd(Activity))

  summarized_progeny_scores_proportions <- progeny_scores_clusters_proportions %>%
    group_by(Pathway, group) %>%
    summarise(avg = mean(Activity),
              std = sd(Activity))

  pooledSD_bypathway <- summarized_progeny_scores_proportions %>%
    group_by(Pathway) %>%
    summarise(poolSD = pooledSD(std))
  total_celltypes <- length(unique(seurat_obj$cell_type))
  pooledSD_bypathway <- pooledSD_bypathway[rep(seq_len(nrow(pooledSD_bypathway)), each = total_celltypes), ]

  summarized_progeny_scores_proportions <- cbind(pooledSD_bypathway,summarized_progeny_scores_proportions$group,summarized_progeny_scores_proportions$avg)
  names(summarized_progeny_scores_proportions) <- c("pathway","poolSD","group","avg")
  summarized_progeny_scores_proportions$adj_avg <- summarized_progeny_scores_proportions$avg/summarized_progeny_scores_proportions$poolSD


  Pathways_adjavg <- summarized_progeny_scores_proportions %>%
    dplyr::select(-c(poolSD,avg)) %>%
    spread(pathway, adj_avg) %>%
    data.frame(row.names = 1, check.names = FALSE)

  Pathways_avg <- summarized_progeny_scores %>%
    dplyr::select(-c(std)) %>%
    spread(Pathway, avg) %>%
    data.frame(row.names = 1, check.names = FALSE)

  # Return everything
  progeny_handled <- list("scores_bycell" = data.frame(progeny_scores_clusters), "proportionadjusted_scores_bycell" = data.frame(progeny_scores_clusters_proportions),
                          "scores_bygroup" = data.frame(summarized_progeny_scores), "proportionadjusted_scores_bygroup" = data.frame(summarized_progeny_scores_proportions),
                          "pathways_bygroup" = Pathways_avg, "proportionadjusted_pathways_bygroup" = Pathways_adjavg)
  return(progeny_handled)
}
