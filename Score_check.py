import numpy as np
import scanpy as sc
import pandas as pd
import seaborn as sns
import decoupler as dc
import matplotlib.pyplot as plt
from matplotlib_venn import venn2,venn3
from mpl_toolkits.axes_grid1 import make_axes_locatable
import anndata as ad
import scipy
import os
import warnings
from collections import Counter

from G03_Annotation_Functions import *

warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings("ignore", message="No data for colormapping provided via 'c'.*")

sc.set_figure_params(figsize=(6, 6), dpi_save=300)
print("Scanpy version:",sc.__version__)

### TO CHECK WHEN CHANGING SAMPLES ###
Experiment='Wouter21_SING_CB'
### TO CHECK WHEN CHANGING SAMPLES ###

refDir = '/mnt/etemp/ahrmad/wouter/refs'
#refDir = './refs'

if not os.path.exists(refDir):
    raise ValueError('Check ref folder paths')
gene_set_table = f'{refDir}/table_s8_summary.txt'

#Epi subtypes Supp gene lists
gset_df = pd.read_csv(gene_set_table, sep='\t')

post_ext = '_post.h5ad'
lum_ext = '_annot_Lum.h5ad'

#research score
adata_epi_NoSV = sc.read_h5ad(f'{Experiment}{lum_ext}')

cor = scipy.stats.pearsonr(adata_epi_NoSV.obs['log1p_n_genes_by_counts'],adata_epi_NoSV.obs['log1p_total_counts']).statistic
sc.pl.scatter(adata_epi_NoSV,x='log1p_n_genes_by_counts',y='log1p_total_counts',size=400000/adata_epi_NoSV.n_obs,alpha=0.8,color="Epi_Luminal_1",save=f'{Experiment}_NgeneNcount_L1GeneScore',title=f'{Experiment}\nPearsonCorr:{cor:.2f}',show=False)

#for score in ['Epi_Basal_1','Epi_Luminal_1','Epi_Luminal_2Psca']:
#    print(scipy.stats.pearsonr(adata_epi_NoSV.obs[score],adata_epi_NoSV.obs['log1p_n_genes_by_counts']))
#    print(scipy.stats.pearsonr(adata_epi_NoSV.obs[score],adata_epi_NoSV.obs['log1p_total_counts']))
#print(scipy.stats.pearsonr(adata_epi_NoSV.obs['Epi_Luminal_1'],adata_epi_NoSV.obs['Epi_Luminal_2Psca']))

# batch
for sctype in [('L1','Epi_Luminal_1'),('L2','Epi_Luminal_2Psca'),('B','Epi_Basal_1')]:
    for qctype in [('Ngene','log1p_n_genes_by_counts'),('Ncounts','log1p_total_counts')]:
        sctyp,scty = sctype[1],sctype[0]
        qctyp,qcty = qctype[1],qctype[0]
        cor = scipy.stats.pearsonr(adata_epi_NoSV.obs[sctyp],adata_epi_NoSV.obs[qctyp]).statistic
        sc.pl.scatter(adata_epi_NoSV,x=sctyp,y=qctyp,size=400000/adata_epi_NoSV.n_obs,alpha=0.8,color="batch",save=f'{Experiment}_{scty}{qcty}_batch_GeneScore',title=f'{Experiment}\nPearsonCorr:{cor:.2f}',show=False)

# TypeDiffAmbig
for sctype in [('L1','Epi_Luminal_1'),('L2','Epi_Luminal_2Psca'),('B','Epi_Basal_1')]:
    for qctype in [('Ngene','log1p_n_genes_by_counts'),('Ncounts','log1p_total_counts')]:
        sctyp,scty = sctype[1],sctype[0]
        qctyp,qcty = qctype[1],qctype[0]
        cor = scipy.stats.pearsonr(adata_epi_NoSV.obs[sctyp],adata_epi_NoSV.obs[qctyp]).statistic
        sc.pl.scatter(adata_epi_NoSV,x=sctyp,y=qctyp,size=400000/adata_epi_NoSV.n_obs,alpha=0.8,color="TypeDiffAmbig",save=f'{Experiment}_{scty}{qcty}_TypeDiffAmbig_GeneScore',title=f'{Experiment}\nPearsonCorr:{cor:.2f}',show=False)



plt.hist(np.abs(adata_epi_NoSV.obs['Epi_Luminal_1']-adata_epi_NoSV.obs['Epi_Luminal_2Psca']),bins=100)
plt.xlabel('L1/L2 Difference')
plt.ylabel('Frequency')
plt.title(f'Histogram of L1/L2 Differences {Experiment}')
plt.savefig(f'TypeDiffAmbigDistribution_{Experiment}.pdf')
plt.close()
#research score ORA

adata_epi_NoSV = sc.read_h5ad(f'{Experiment}{epi_ext[:-5]}_Subtypes_ORA.h5ad')
del adata_epi_NoSV.obs['Epi_Basal_1']
del adata_epi_NoSV.obs['Epi_Luminal_1']
del adata_epi_NoSV.obs['Epi_Luminal_2Psca']

#for score in ['Epi_Basal_1','Epi_Luminal_1','Epi_Luminal_2Psca']:
#    print(scipy.stats.pearsonr(adata_epi_NoSV.obs[score],adata_epi_NoSV.obs['log1p_n_genes_by_counts']))
#    print(scipy.stats.pearsonr(adata_epi_NoSV.obs[score],adata_epi_NoSV.obs['log1p_total_counts']))
#print(scipy.stats.pearsonr(adata_epi_NoSV.obs['Epi_Luminal_1'],adata_epi_NoSV.obs['Epi_Luminal_2Psca']))

# batch
for sctype in [('L1','Epi_Luminal_1'),('L2','Epi_Luminal_2Psca'),('B','Epi_Basal_1')]:
    for qctype in [('Ngene','log1p_n_genes_by_counts'),('Ncounts','log1p_total_counts')]:
        sctyp,scty = sctype[1],sctype[0]
        qctyp,qcty = qctype[1],qctype[0]
        cor = scipy.stats.pearsonr(adata_epi_NoSV.obsm['ora_estimate'][sctyp],adata_epi_NoSV.obs[qctyp]).statistic
        sc.pl.scatter(adata_epi_NoSV,x=sctyp,y=qctyp,size=400000/adata_epi_NoSV.n_obs,alpha=0.8,color="batch",save=f'{Experiment}_{scty}{qcty}_batch_ORA',title=f'{Experiment}\nPearsonCorr:{cor:.2f}',show=False)



