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
The first step is to quantify the activity of transcription factors and/or pathways from the scRNAseq gene expression data. To better understand the DoRothEA confidence scores used in the `run_dorothea()` function, please refer to the DoRothEA package from the Saez Lab: https://github.com/saezlab/dorothea/
```R
pbmc <- run_dorothea(seurat_obj = pbmc, conf_scores = c("A","B","C","D"), cores = 16)
pbmc <- run_progeny(pbmc, c("A","B","C","D"))
