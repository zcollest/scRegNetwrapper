import os
import subprocess
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
import glob
import json
import base64
import zlib
import matplotlib.pyplot as plt
import matplotlib as mpl
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from dask.diagnostics import ProgressBar
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2
from pyscenic.cli.utils import load_signatures
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
 
def pyscenic_wrapper(dir, anndata_path, loom_path, tfs_path, rank_db_path, motif_path, output_loom_path):
    anndata = sc.read(anndata_path,         
    var_names='gene_symbols',
    cache=False)
    anndata.var_names_make_unique()
    anndata = anndata[anndata.obs['donor'] == 'nst9']
    # create basic row and column attributes for the loom file:
    row_attrs = {
        "Gene": np.array(anndata.var_names) ,
    }
    col_attrs = {
        "CellID": np.array(anndata.obs_names) ,
        "nGene": np.array( np.sum(anndata.X.transpose()>0 , axis=0)).flatten() ,
        "nUMI": np.array( np.sum(anndata.X.transpose() , axis=0)).flatten() ,
    }
    lp.create( loom_path, anndata.X.transpose(), row_attrs, col_attrs)
    
    os.chdir(dir)
    ##### SCENIC #### 
    # STEP 1: GRN BOOST
    cmd1 = "/usr/local/bin/pyscenic grn " + loom_path + " " + tfs_path + " -o anndata_adj.csv --num_workers 20"
    os.system(cmd1)

    # STEP 2/3: Regulon prediction aka cisTarget from CL
    # ranking databases
    db_names = ' '.join( glob.glob(rank_db_path) )
    
    # running it (prunes regulons for which < 80% genes can be matched to the ranking database)
    cmd2 = "/usr/local/bin/pyscenic grn anndata_adj.csv " + db_names + " --annotations_fname " + motif_path + " --expression_mtx_fname " + loom_path + " --output anndata_reg.csv --mask_dropouts --num_workers 20"
    os.system(cmd2)

    # STEP 4: AUCell
    cmd3 = "/usr/local/bin/pyscenic aucell " + loom_path + " anndata_reg.csv --output " + output_loom_path + " --num_workers 20"
    os.system(cmd3)
  
  
def setup_loom_wrapper(anndata_path, loom_path):
    anndata = sc.read(anndata_path,         
    var_names='gene_symbols',
    cache=False)
    anndata.var_names_make_unique()
    anndata = anndata[anndata.obs['donor'] == 'nst9']
    # create basic row and column attributes for the loom file:
    row_attrs = {
        "Gene": np.array(anndata.var_names) ,
    }
    col_attrs = {
        "CellID": np.array(anndata.obs_names) ,
        "nGene": np.array( np.sum(anndata.X.transpose()>0 , axis=0)).flatten() ,
        "nUMI": np.array( np.sum(anndata.X.transpose() , axis=0)).flatten() ,
    }
    lp.create( loom_path, anndata.X.transpose(), row_attrs, col_attrs)
    
def grn_wrapper(dir, loom_path, tfs_path):
    os.chdir(dir)
    ##### SCENIC #### 
    # STEP 1: GRN BOOST
    cmd1 = "/usr/local/bin/pyscenic grn " + loom_path + " " + tfs_path + " -o anndata_adj.csv --num_workers 20"
    os.system(cmd1)
    
def cistarget_wrapper(dir, loom_path, rank_db_path, motif_path):
    os.chdir(dir)
   # STEP 2/3: Regulon prediction aka cisTarget from CL
    # ranking databases
    db_names = ' '.join( glob.glob(rank_db_path) )
    
    # running it (prunes regulons for which < 80% genes can be matched to the ranking database)
    cmd2 = "/usr/local/bin/pyscenic grn anndata_adj.csv " + db_names + " --annotations_fname " + motif_path + " --expression_mtx_fname " + loom_path + " --output anndata_reg.csv --mask_dropouts --num_workers 20"
    os.system(cmd2)
    
def aucell_wrapper(dir, loom_path, output_loom_path):
    # STEP 4: AUCell
    os.chdir(dir)
    cmd3 = "/usr/local/bin/pyscenic aucell " + loom_path + " anndata_reg.csv --output " + output_loom_path + " --num_workers 20"
    os.system(cmd3)
    
    
