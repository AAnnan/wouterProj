##Preprocessing

# Purpose:
#    Import fragments into Snapatac2
#    Perform QC

from scipy.stats import median_abs_deviation
import snapatac2 as snap
import pandas as pd
import numpy as np
import os
print('snapatac2 version:',snap.__version__)

# Input files
path10x = '/mnt/ndata/daniele/wouter/Processed/CellRangerArc/'
pathRes = '/mnt/etemp/ahrmad/wouter/ATAC'
exp10x = [d for d in os.listdir(path10x) if d.startswith('WK') and os.path.isdir(os.path.join(path10x, d))]

# QC filters
minFrags = 1000
minTSSe = 5

# Hard Threshold (cells are filtered out): maxFrags = max(MADs_maxFrags,Alt_maxFrags)
# Min Threshold (cells are labeled): maxFrags = Alt_maxFrags
Alt_maxFrags = 80000 
nmads = 6 

for exp in exp10x:
	print(f'IMPORT AND FILTER')
	print(f'{exp}')
	#Import ATAC fragments from 10x pipeline
	data = snap.pp.import_data(
		fragment_file=os.path.join(path10x,exp,'outs/atac_fragments.tsv.gz'),
		chrom_sizes=snap.genome.mm10,
		whitelist=os.path.join(pathRes,"whitelists",exp),
		sorted_by_barcode=False,min_num_fragments=0,
		tempdir=pathRes
		)

	#Filter cells
	snap.metrics.tsse(data, gene_anno=snap.genome.mm10, inplace=True)
	MADs_maxFrags = np.median(data.obs["n_fragment"]) + nmads * median_abs_deviation(data.obs["n_fragment"])
	maxFrags = max(MADs_maxFrags,Alt_maxFrags)
	snap.pp.filter_cells(data, min_counts=minFrags,max_counts=maxFrags, min_tsse=minTSSe, inplace=True)
	
	# QC metrics plots
	snap.metrics.frag_size_distr(data,add_key='frag_size_distr', inplace=True)
	snap.pl.tsse(data, min_fragment=0, width=750, height=600, interactive=False, show=False, out_file=exp+'_tsseFiltered.pdf')
	
	#Export fragments and anndata
	data.obs['Exp'] = pd.Categorical([exp]*data.n_obs)
	snap.ex.export_fragments(data, groupby='Exp', prefix='', suffix='.tsv.gz')
	data.write(filename=f'{exp}_filt.h5ad')
	del data