# scRegNetwrapper
Wrapper functions for single cell regulatory network analysis (DoRothEA, PROGENy) as well as tools for downstream analysis.

## Installation
 ```R 
library(devtools)
library(reticulate)
use_python("/usr/bin/python3") # specify whatever python path you'd prefer
library(progeny) # loading progeny loads global variables that are required for some functions
install_github("zcollest/scRegNetwrapper")
```
In order to use the Python functions that are contained in this package (pySCENIC analysis), it is necessary to ensure that the proper python dependencies are installed. At the root of this package is a `requirements.txt` file with the Python packages and versions required. Download that file and use the following shell script to download those packages -- note: make sure that these packages are being installed in the python environment that you specified above. 
```bash
pip install -r requirements.txt
```
Lastly, load the python functions into the R session with the following R code:

```R
run_pyscenic_path <- system.file("python","run_pySCENIC.py", package="scRegNetwrapper")
source_python(run_pyscenic_path)
```

## Input 
The calculation of DoRothEA/PROGENy scores requires a Seurat object of scRNAseq gene expression data. Downstrean analyses utilize data that is outputted from the `handle_dorothea_scores()` and `handle_progeny_scores()` functions. 

## Tutorial

To obtain more information about any specific function, run `?"function_name"()` in your R session.

### Quantifying transcription factor and pathway activity
The first step is to quantify the activity of transcription factors and/or pathways from the scRNAseq gene expression data. To better understand the DoRothEA confidence scores used in the `run_dorothea()` function, please refer to the DoRothEA package from the Saez Lab: https://github.com/saezlab/dorothea/. To better understand the inputs for the `run_progeny()` function, please refer to the PROGENy package from the Saez Lab: https://github.com/saezlab/progeny/. If interested in inserting custom regulons or pathways, the format must be the same as the regulons and pathways from DoRothEA and PROGENy, respectively. 

```R
# Running dorothea and progeny
pbmc <- run_dorothea(seurat_obj = pbmc, conf_scores = c("A","B","C","D"), cores = 16)
pbmc <- run_progeny(seurat_obj = pbmc, num_genes = 500, organism = "Human")

# Using custom transcription factor regulons and pathway information 
pbmc <- custom_pathways_calc(seurat_obj = pbmc, regulons = regulons, assay_name = "custom_regulon_scores")
pbmc <- custom_regulons_calc(seurat_obj = pbmc, pathways = pathways, num_genes = 100, organism = "Human", assay_name = "custom_pathway_scores") 
```
### Handling scores for downstream analysis
These functions organize the cell-wise transcription factor and pathway scores into various data frames that can be used for downstream analyses. The functions require a `comparison_feature` argument, which is important for two reasons: <br>
1. It defines the output of the heatmap used for downstram analyses
2. It adjusts the cell-wise scores based on the proportion of cells in each category of the comparison feature (i.e. healthy vs disease).
```R
# Handling scores 
tf_scores <- handle_dorothea_scores(seurat_obj = pbmc, comparison_feature = pbmc@meta.data$indication, topTFs = 30)
pathway_scores <- handle_progeny_scores(seurat_obj = pbmc, comparison_feature = pbmc@meta.data$indication)
```

### Visualizing the results with heatmaps 
There are multiple ways that you can visualize these scores. Perhaps one of the best ways is to look at a heatmap of transcription factor or pathway activity for a given comparison (cluster-wise comparisons, healthy vs disease comparisons, etc.). For example, this can be useful for visualizing cell-type heterogeneity from the perspective of transcription factors and pathways. The `downstream_heatmap` allows a user to quickly and easily output a basic heatmap, but the source code can easily be modified to produce heatmaps more aligned to the user's preferences.

```R
downstream_heatmap(data = dorothea_scores$proportionadjusted_tfs_bygroup, title = "progeny pathways, by indication (healthy vs disease)")
downstream_heatmap(data = progeny_scores$proportionadjusted_pathways_bygroup, title = "progeny pathways, by indication (healthy vs disease)")
```

### Effect size calculations (Cohen's D) 
To determine statistical significance for comparisons between two groups (i.e. healthy vs disease), the user can calculate effect sizes (Cohen's D scores) for each transcription factor and pathway. The output is a data frame with columns that include transcription factor / pathway name and its corresponding Cohen's D score. 

```R
cohend_dorothea <- tf_effsize_calc(data = dorothea_scores$proportionadjusted_scores_bycell)
cohend_progeny <- pathway_effsize_calc(data = progeny_scores$proportionadjusted_scores_bycell)
```

### Finding transcription factors associated with a vector of genes  
Given an input vector of gene names, this function searches the regulons of the transcriptions for those genes. Optional arguments include effect size data as well as transcription factor activity summarized by comparison group (one of the outputs from the `handle_dorothea_scores()` function.  It returns a list of data frames in which each data frame is a target gene and returns the associated transcription factor as well as the effect size for each associated transcription factor with respect to the comparison group. 

```R
gene_vector <- c("CHD8","DOK2","RGS4","FOCAD-AS1","PYGB") # these genes were generated from a random gene set generator 
gene_associations <- find_associated_TFs(gene_vector=gene_vector, tf_data_bygroup=dorothea_scores$proportionadjusted_scores_bygroup, effect_size_data=cohend_dorothea)
```

### Correlation Analysis
Another interesting way to visualize the results is to look at correlations between TF activity and pathway activities. Similar to the `downstream_heatmap` function, this function allows a user to quickly and easily output a basic correlation matrix, but the source code can easily be modified to produce matrices more aligned to the user's preferences. There is also the option not to render the plot and just return the correlations in data frame form. 

```R
corr <- correlation_analysis(tf_data = tf_data ,pathway_data = pathway_data, return_corr_data = TRUE, render_plot = TRUE)
```
