import os
import csv
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from pyscenic.utils import modules_from_adjacencies, load_motifs
from arboreto.utils import load_tf_names
from pyscenic.cli.utils import load_signatures
from pyscenic.export import add_scenic_metadata
from pyscenic.rss import regulon_specificity_scores

def scenic_results_wrapper(dir, output_loom_path, anndata_path, comparison_feature, regulon_path):
  os.chdir(dir)
  lf = lp.connect( output_loom_path, mode='r', validate=False )
  adata_output = sc.read( output_loom_path, validate=False)
  adata_original = sc.read(anndata_path,         
    var_names='gene_symbols',
    cache=False)
  adata_original.var_names_make_unique()
  adata_output.obs['cell_type'] = adata_original.obs[comparison_feature] 
  sig = load_signatures(regulon_path)
  # getting dict of regulons and saving it
  regulons = {}
  for i,r in pd.DataFrame(lf.ra.Regulons,index=lf.ra.Gene).iteritems():
      regulons[i] =  list(r[r==1].index.values)
  a_file = open('regulons.csv', "w")
  writer = csv.writer(a_file)
  for key, value in regulons.items():
      writer.writerow([key, value])
  a_file.close()
  
  # getting AUC matrix and saving it
  auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
  auc_mtx.to_csv('auc_mtx.csv')
  adata_output = add_scenic_metadata(adata_output, auc_mtx, sig)
  cellAnnot = pd.concat(
      [
          pd.DataFrame( adata_output.obs['cell_type'], index=lf.ca.CellID )
      ],
      axis=1
  )
  cellAnnot.to_csv('cellAnnot.csv')
  
  # calculate RSS
  rss_cellType = regulon_specificity_scores( auc_mtx, cellAnnot['cell_type'] )
  rss_cellType.to_csv('RSS.csv')
  
