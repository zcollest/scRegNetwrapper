from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_rss
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text
import os

def plot_RSS_wrapper(dir, cellAnnot, rss_cellType, title):
  os.chdir(dir)
  cats = sorted(list(set(cellAnnot['cell_type'])))
  
  fig = plt.figure(figsize=(15, 8))
  for c,num in zip(cats, range(1,len(cats)+1)):
      x=rss_cellType.T[c]
      ax = fig.add_subplot(2,5,num)
      plot_rss(rss_cellType, c, top_n=5, max_n=None, ax=ax)
      ax.set_ylim( x.min()-(x.max()-x.min())*0.05 , x.max()+(x.max()-x.min())*0.05 )
      for t in ax.texts:
          t.set_fontsize(12)
      ax.set_ylabel('')
      ax.set_xlabel('')
      adjust_text(ax.texts, autoalign='xy', ha='right', va='bottom', arrowprops=dict(arrowstyle='-',color='lightgrey'), precision=0.001 )
   
  fig.text(0.5, 0.0, 'Regulon', ha='center', va='center', size='x-large')
  fig.text(0.00, 0.5, 'Regulon specificity score (RSS)', ha='center', va='center', rotation='vertical', size='x-large')
  plt.tight_layout()
  plt.rcParams.update({
      'figure.autolayout': True,
          'figure.titlesize': 'large' ,
          'axes.labelsize': 'medium',
          'axes.titlesize':'large',
          'xtick.labelsize':'medium',
          'ytick.labelsize':'medium'
          })
  plt.savefig(title, dpi=600, bbox_inches = "tight")

    
