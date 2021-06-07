# scRegNetwrapper
Wrapper functions for single cell regulatory network analysis (DoRothEA, PROGENy) as well as tools for downstream analysis.

## Installation
 ```R 
library(devtools)
install_github("zcollest/scRegNetwrapper")
library(scRegNetwrapper)
```

## Input 
The calculation of DoRothEA/PROGENy scores requires a Seurat object of scRNAseq gene expression data. Downstrean analyses utilize data that is outputted from the `handle_dorothea_scores()` and `handle_progeny_scores()` functions. 

## Tutorial

To obtain more information about any specific function, run `?"function_name"()` in your R session.

#### Quantifying transcription factor and pathway activity
The first step is to quantify the activity of transcription factors and/or pathways from the scRNAseq gene expression data. To better understand the DoRothEA confidence scores used in the `run_dorothea()` function, please refer to the DoRothEA package from the Saez Lab: https://github.com/saezlab/dorothea/. To better understand the inputs for the `run_progeny()` function, please refer to the PROGENy package from the Saez Lab: https://github.com/saezlab/progeny/. 

```R
# Running dorothea and progeny
pbmc <- run_dorothea(seurat_obj = pbmc, conf_scores = c("A","B","C","D"), cores = 16)
pbmc <- run_progeny(seurat_obj = pbmc, num_genes = 500, organism = "Human")

# Using custom transcription factor regulons and pathway information 
pbmc <- custom_pathways_calc(seurat_obj = pbmc, regulons = regulons, assay_name = "custom_regulon_scores")
pbmc <- custom_regulons_calc(seurat_obj = pbmc, pathways = pathways, num_genes = 100, organism = "Human", assay_name = "custom_pathway_scores") 
```
#### Handling scores for downstream analysis
These functions organize the cell-wise transcription factor and pathway scores into various data frames that can be used for downstream analyses. The functions require a `comparison_feature` argument, which is important for two reasons: <br>
1. It defines the output of the heatmap used for downstram analyses
2. It adjusts the cell-wise scores based on the proportion of cells in each category of the comparison feature (i.e. healthy vs disease).
```R
# Handling scores 
tf_scores <- handle_dorothea_scores(seurat_obj = pbmc, comparison_feature = pbmc@meta.data$indication, topTFs = 30)
pathway_scores <- handle_progeny_scores(seurat_obj = pbmc, num_genes = 500, organism = "Human")

```
