##Processing

# Purpose:
#    Doublet Analysis based on peaks from MACS3

import snapatac2 as snap
import pandas as pd
import numpy as np
import os
print('snapatac2 version:',snap.__version__)

# Input files
path10x = '/mnt/ndata/daniele/wouter/Processed/CellRangerArc/'
pathRes = '/mnt/etemp/ahrmad/wouter/batch_ATAC'
refDir = '/mnt/etemp/ahrmad/wouter/refs'
exp10x = [d for d in os.listdir(path10x) if d.startswith('WK') and os.path.isdir(os.path.join(path10x, d))]

if not os.path.exists(refDir+'/csv'):
    os.makedirs(refDir+'/csv')
    print(f"Directory {refDir+'/csv'} created.")
else:
    print(f"Directory {refDir+'/csv'} already exists.")

for exp in exp10x:
	print(f'DOUBLET ANALYSIS')

	# Read in QC'd Anndata file
	pm = snap.read(exp +'_filt.h5ad').to_memory()
	# Create cell-by-peak matrix from it with IOM'd Macs3 peak file
	data = snap.pp.make_peak_matrix(pm, peak_file=f'{exp}_macs3_Q01/{exp}_ITMPeaks.bed')
	# Copy obsm and ref seqs
	data.obsm = pm.obsm.copy()
	data.uns['reference_sequences'] = pm.uns['reference_sequences'].copy()
	del pm
	
	# Feature selection 
	snap.pp.select_features(data, n_features=80000, inplace=True, blacklist="mm10-blacklist.v2.bed.gz")
	# Doublet Detection
	snap.pp.scrublet(data, features='selected', n_comps=15, sim_doublet_ratio=2.0, expected_doublet_rate=0.1)	
	# Doublet Annotation
	doub = snap.pp.filter_doublets(data, probability_threshold=0.7, inplace=False)
	data.obs['doublet_class'] = np.where(doub, 'singlet', 'doublet')
	
	# Export Anndata
	data.write(filename=f'{exp}_filt_PEAKS_Self.h5ad')
	# Export CSV containing Doublet scores and annotation
	df = pd.DataFrame({'sample': exp, 'obs_names': data.obs_names, 'doublet_score': data.obs['doublet_probability'], 'doublet_class': data.obs['doublet_class']})
	df.to_csv(f'{refDir}/csv/{exp}_doublet_scores_ATACPeaks_Self.csv', index=False)

	print(f'{exp} QC DONE and files created')