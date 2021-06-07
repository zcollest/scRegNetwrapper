# scRegNetwrapper
Wrapper functions for single cell regulatory network analysis (DoRothEA, PROGENy) as well as tools for downstream analysis.

## Installation
 ```R 
library(devtools)
install_github("zcollest/scRegNetwrapper")
library(scRegNetwrapper)
```

## Inputs
The calculation of DoRothEA/PROGENy scores requires a Seurat object of scRNAseq gene expression data. Downstrean analyses utilize data that is outputted from the `handle_dorothea_scores()` and `handle_progeny_scores()` functions. 

