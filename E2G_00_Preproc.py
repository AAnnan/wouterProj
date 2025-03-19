#Pseudobulk fragment files and their corresponding *.tbi index files in the same directory 

#get cell names for L1/L2 cd28 and intact from sc rna  (4 files)
# export raw exp matrix as csv for these cells (4 files)
# import tsv files for intact


import matplotlib.pyplot as plt
import snapatac2 as snap
import scanpy as sc
import numpy as np
import polars as pl
import pandas as pd
import anndata as ad
from matplotlib_venn import venn3
import os
import csv

Experiment='Wouter21_SING_CB'
lum_ext = '_annot_Lum_newSigs.h5ad'
base_dir="/mnt/etemp/ahrmad/wouter/"

adata = sc.read_h5ad(base_dir + '/RNA_CB/annot_other_norm/' + Experiment + lum_ext)

adata_L1 = adata[np.isin(adata.obs['Annotation'],['L1'])].copy()
adata_L1_CD28 = adata_L1[np.isin(adata_L1.obs['timePoint'],['Day28'])].copy()
adata_L1_Intact = adata_L1[np.isin(adata_L1.obs['timePoint'],['Intact'])].copy()

adata_L2 = adata[np.isin(adata.obs['Annotation'],['L2'])].copy()
adata_L2_CD28 = adata_L2[np.isin(adata_L2.obs['timePoint'],['Day28'])].copy()
adata_L2_Intact = adata_L2[np.isin(adata_L2.obs['timePoint'],['Intact'])].copy()

for name,ad_obj in zip(['L1_CD28','L2_CD28','L2_Intact','L1_Intact'],[adata_L1_CD28,adata_L2_CD28,adata_L2_Intact,adata_L1_Intact]):
  print(f'Exporting {name}')
  cell_names = ad_obj.obs.index.tolist()

  with open(f'{name}_CellNames.lst', 'w') as f:
      for cell in cell_names:
          f.write(f"{cell}\n")
  adata_reduced = ad.AnnData(X=ad_obj.layers['No_norm'])
  # Save the reduced AnnData object to a new .h5ad file
  adata_reduced.write(f'{name}_RawCountMatrix.h5ad')

# GET INTACT/CD28 tsv fragments from QC
adata_CD28 = snap.read(f'{base_dir}/Wouter21_SING/WK-1585_Castrate_Day28_AP_BL6_qcTOREMOVE.h5ad').to_memory()

adata_L1_CD28 = adata_L1[np.isin(adata_L1.obs['timePoint'],'Day28')].copy()
adata_CD28_L1 = adata_CD28[np.isin(adata_CD28.obs.index,adata_L1_CD28.obs.index)].copy()
adata_CD28_L1.obs['TYPE'] = ['CD28_L1']*adata_CD28_L1.n_obs
snap.ex.export_fragments(adata_CD28_L1, groupby='TYPE', prefix='CD28_L1_FRAGS', suffix='.tsv.gz')

adata_L2_CD28 = adata_L2[np.isin(adata_L2.obs['timePoint'],'Day28')].copy()
adata_CD28_L2 = adata_CD28[np.isin(adata_CD28.obs.index,adata_L2_CD28.obs.index)].copy()
adata_CD28_L2.obs['TYPE'] = ['CD28_L2']*adata_CD28_L2.n_obs
snap.ex.export_fragments(adata_CD28_L2, groupby='TYPE', prefix='CD28_L2_FRAGS', suffix='.tsv.gz')

adata_L1_I = adata_L1[np.isin(adata_L1.obs['timePoint'],'Intact')].copy()
adata_L2_I = adata_L2[np.isin(adata_L2.obs['timePoint'],'Intact')].copy()

for a,b in [('WK-1501_BL6_INTACT_AP_Test3','WK-1501-Test3'), ('WK-1585_INTACT_AP_BL6_Citrate','WK-1585-Citrate'), ('WK-1585_INTACT_AP_BL6_Contrl','WK-1585-Contrl')]:
	adata_I = snap.read(f'../../Wouter21_SING/{a}_qcTOREMOVE.h5ad').to_memory()
	adata_I_L1 = adata_I[np.isin(adata_I.obs.index,adata_L1_I.obs.index)].copy()
	adata_I_L1.obs['TYPE'] = ['I_L1']*adata_I_L1.n_obs
	snap.ex.export_fragments(adata_I_L1, groupby='TYPE', prefix=f'{b}_I_L1_FRAGS', suffix='.tsv.gz')

	adata_I_L2 = adata_I[np.isin(adata_I.obs.index,adata_L2_I.obs.index)].copy()
	adata_I_L2.obs['TYPE'] = ['I_L2']*adata_I_L2.n_obs
	snap.ex.export_fragments(adata_I_L2, groupby='TYPE', prefix=f'{b}_I_L2_FRAGS', suffix='.tsv.gz')

zcat *_I_*L1*.tsv.gz | pigz -p 24 --fast --stdout > All_Intact_Samples_L1.tsv.gz
zcat *_I_*L2*.tsv.gz | pigz -p 12 --fast --stdout > All_Intact_Samples_L2.tsv.gz

mv All_Intact_Samples_L1.tsv.gz All_Intact_Samples_I_L1.tsv.gz
mv All_Intact_Samples_L2.tsv.gz All_Intact_Samples_I_L2.tsv.gz

# REIMPORT FOR FILTERING

import matplotlib.pyplot as plt
import snapatac2 as snap
import scanpy as sc
import numpy as np
import polars as pl
import pandas as pd
import anndata as ad
from matplotlib_venn import venn3
import os
import csv

for name,fgt in zip(['L1_CD28','L2_CD28','L2_Intact','L1_Intact'],["CD28_L1_FRAGSCD28_L1.tsv.gz","CD28_L2_FRAGSCD28_L2.tsv.gz","All_Intact_Samples_I_L2.tsv.gz","All_Intact_Samples_I_L1.tsv.gz"]):
  print(f'Exporting {name} ATAC')
  atac_ad = snap.pp.import_data(fgt, snap.genome.mm10, min_num_fragments=0, sorted_by_barcode=True,
   chrM=['chrM', 'M'], shift_left=0, shift_right=0, chunk_size=50000, tempdir='.')
  cell_names = pd.read_csv(f'{name}_CellNames.lst', header=None)
  atac_ad = atac_ad[np.isin(atac_ad.obs.index,cell_names)].copy()
  atac_ad.obs['Exp'] = [name]*atac_ad.n_obs
  snap.ex.export_fragments(atac_ad, groupby='Exp', prefix='', suffix='_Unsorted.tsv')

for name in L1_CD28 L2_CD28 L2_Intact L1_Intact; do
    sortBed -i "${name}_Unsorted.tsv" > "${name}.tsv"
    bgzip -@ 30 "${name}.tsv"
    tabix -@ 30 -p bed "${name}.tsv.gz"
done

conda create --prefix /mnt/dataFast/ahrmad/E2G/scE2G/envs/sc_e2g --clone /home/annan/miniconda3/envs/sc_e2g

snakemake -j 64
