# scRegNetwrapper

Wrapper functions to simplify popular single cell regulatory network analysis ([DoRothEA](https://github.com/saezlab/dorothea), [PROGENy](https://github.com/saezlab/progeny), [pySCENIC](https://github.com/aertslab/pySCENIC)) as well as tools for downstream analysis and method comparisons.

## Installation
 ```R 
library(devtools)
library(reticulate)
use_python("/usr/bin/python3") # specify whatever python path you'd prefer
library(progeny) # loading progeny loads global variables that are required for some functions
```
In order to use the Python functions that are contained in this package (pySCENIC analysis), it is necessary to ensure that the proper python dependencies are installed. At the root of this package is a `requirements.txt` file with the Python packages and versions required. Download that file and use the following shell script to download those packages -- note: make sure that these packages are being installed in the python environment that you specified above. 

```bash
pip install -r requirements.txt
```
Finally, install and load the package:

```R
install_github("zcollest/scRegNetwrapper")
library(scRegNetwrapper)
```

## Input 
The calculation of TF/pathway scores using DoRothEA and PROGENy requires a Seurat object of scRNAseq gene expression data. Calculating TF scores using pySCENIC,  requires a path to an anndata object - to convert a Seurat object to anndata, it is recommended to use [sceasy](https://github.com/cellgeni/sceasy). 

## Tutorial

To obtain more information about any specific function, run `?"function_name"()` in your R session.

### Quantifying transcription factor and pathway activity using DoRothEA & PROGENy
The first step is to quantify the activity of transcription factors and/or pathways from the scRNAseq gene expression data. To better understand the DoRothEA confidence scores used in the `run_dorothea()` function, please refer to the [DoRothEA](https://github.com/saezlab/dorothea/) package from the Saez Lab. To better understand the inputs for the `run_progeny()` function, please refer to the [PROGENy](https://github.com/saezlab/progeny/) package from the Saez Lab. If interested in inserting custom regulons or pathways, the format must be the same as the regulons and pathways from DoRothEA and PROGENy, respectively. 

```R
# Running dorothea and progeny
pbmc <- run_dorothea(seurat_obj = pbmc, conf_scores = c("A","B","C","D"), cores = 16)
pbmc <- run_progeny(seurat_obj = pbmc, num_genes = 500, organism = "Human")

# Using custom transcription factor regulons and pathway information 
pbmc <- custom_pathways_calc(seurat_obj = pbmc, regulons = regulons, assay_name = "custom_regulon_scores")
pbmc <- custom_regulons_calc(seurat_obj = pbmc, pathways = pathways, num_genes = 100, organism = "Human", assay_name = "custom_pathway_scores") 
```
### Constructing regulons and quantifying transcription factor activity using pySCENIC
This package contains functions to easily run the pySCENIC pipeline to quantify TF activity from scRNAseq gene expression data. While there is an [R implementation of the SCENIC pipeline](https://github.com/aertslab/SCENIC), the pySCENIC implementation is faster and more updated. To better undertand the required inputs and computational steps involved in the pySCENIC pipeline, please refer to the [protocols paper](https://www.nature.com/articles/s41596-020-0336-2) and the [tutorial website](https://pyscenic.readthedocs.io/en/latest/index.html). 

```R
# Run each step of the pySCENIC pipeline
run_loom_setup(anndata_path = "pbmc.h5ad", loom_path = "pbmc.loom")
run_grn(dir = "path/to/output/results", loom_path = "pbmc.loom" , tfs_path = "path/to/input/TFlist.txt", adjacencies_fname = "pbmc_adj.csv")
run_cistarget(dir = "path/to/output/results", adjacencies_fname = "pbmc_adj.csv", loom_path = "pbmc.loom", rank_db_path ="path/to/rankdb.feather", motif_path = "path/to/motif.tbl", regulons_fname = "pbmc_reg.csv")
run_aucell(dir = "path/to/output/results",regulons_fname = "pbmc_reg.csv",loom_path = "pbmc.loom",output_loom_path = "pbmc_output.loom")
```

### Handling results for downstream analysis
With respect to the outputs from DoRothEA and PROGENy, there are functions in this package that organize the cell-wise transcription factor and pathway activity scores into various data frames that can be used for downstream analyses. The functions require a `comparison_feature` argument, which is important for two reasons: <br>
1. It defines the output of the heatmap used for downstram analyses
2. It adjusts the cell-wise scores based on the proportion of cells in each category of the comparison feature (i.e. healthy vs disease).

The outputs from the pySCENIC pipeline are slightly different. The function also requires a `comparison_feature` argument, and the outputted results include: 
1. Data frame of constructed regulons
2. Matrix of cell-wise AUC scores for each TF 
3. A data frame of annotations based on `comparison_feature`.
4. Data frame of regulon specificity scores.

```R
# Handling DoRothEA scores 
tf_scores <- handle_dorothea_scores(seurat_obj = pbmc, comparison_feature = pbmc@meta.data$indication, topTFs = 30)

# Handling PROGENy scores 
pathway_scores <- handle_progeny_scores(seurat_obj = pbmc, comparison_feature = pbmc@meta.data$indication)

# Handling pySCENIC results
pyscenic_results <- handle_pyscenic_results(dir = "path/to/output", output_loom = "pbmc_output.loom", anndata_path "anndata.h5ad", comparison_feature = "cell_type", regulon_path = "path/to/pySCENIC/ regulons.csv")
```

### Visualizing DoRothEA/PROGENy results with heatmaps 
There are multiple ways that you can visualize these scores. Perhaps one of the best ways is to look at a heatmap of transcription factor or pathway activity for a given comparison (cluster-wise comparisons, healthy vs disease comparisons, etc.). For example, this can be useful for visualizing cell-type heterogeneity from the perspective of transcription factors and pathways. The `downstream_heatmap` allows a user to quickly and easily output a basic heatmap, but the source code can easily be modified to produce heatmaps more aligned to the user's preferences.

```R
downstream_heatmap(data = dorothea_scores$proportionadjusted_tfs_bygroup, title = "progeny pathways, by indication (healthy vs disease)")
downstream_heatmap(data = progeny_scores$proportionadjusted_pathways_bygroup, title = "progeny pathways, by indication (healthy vs disease)")
```

### Visualizing pySCENIC results with RSS Plots
pySCENIC utilizes Regulon Specificity Scores (RSS) to determine how specific each regulon is to a given cell type / comparison feature. As of right now, we cannot render the plot in R since it uses matplotlib (we are working on this!). However, you can save the plot and view it as an image or a PDF.

```R
plot_RSS(dir = '/data', pyscenic_results$cellAnnot, pyscenic_results$RSS, title = "RSS_plot.png")
```


### Effect size calculations (Cohen's D) for DoRothEA and PROGENy results
To determine statistical significance for comparisons between two groups (i.e. healthy vs disease), the user can calculate effect sizes (Cohen's D scores) for each transcription factor and pathway. The output is a data frame with columns that include transcription factor / pathway name and its corresponding Cohen's D score. 

```R
cohend_dorothea <- tf_effsize_calc(data = dorothea_scores$proportionadjusted_scores_bycell)
cohend_progeny <- pathway_effsize_calc(data = progeny_scores$proportionadjusted_scores_bycell)
```

### Correlation Analysis with DoRothEA and PROGENy results
Another interesting way to visualize the results is to look at correlations between TF activity and pathway activities. Similar to the `downstream_heatmap` function, this function allows a user to quickly and easily output a basic correlation matrix, but the source code can easily be modified to produce matrices more aligned to the user's preferences. There is also the option not to render the plot and just return the correlations in data frame form. 

```R
corr <- correlation_analysis(tf_data = tf_data ,pathway_data = pathway_data, return_corr_data = TRUE, render_plot = TRUE)
```


### Finding transcription factors associated with a vector of genes  
Note: right now, this only works with the output from DoRothEA...working to add support for pySCENIC regulons too.  <br> <br>
Given an input vector of gene names, this function searches the regulons of the transcriptions for those genes. Optional arguments include effect size data as well as transcription factor activity summarized by comparison group (one of the outputs from the `handle_dorothea_scores()` function.  It returns a list of data frames in which each data frame is a target gene and returns the associated transcription factor as well as the effect size for each associated transcription factor with respect to the comparison group. 

```R
gene_vector <- c("CHD8","DOK2","RGS4","FOCAD-AS1","PYGB") # these genes were generated from a random gene set generator 
gene_associations <- find_associated_TFs(gene_vector=gene_vector, tf_data_bygroup=dorothea_scores$proportionadjusted_scores_bygroup, effect_size_data=cohend_dorothea)
```

### Note on Code
Some of the code used in this package (particularly for the steps involving quantifying TF and pathway activity) has been adapted from the vignettes for DoRothEA, PROGENy, and pySCENIC. There are links to these packages and their corresponding vignettes at the top of this README.
