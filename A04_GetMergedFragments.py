##Get Merged Fragments

# /!\ Can only be performed after GEX Processing

# Purpose:
#   Merge Selected preprocessed samples in ATAC modality


import snapatac2 as snap
import pandas as pd
import numpy as np
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
print('Snap: ',snap.__version__)

# To change according to samples
# SING for singulator, ENZ for enzymatic digestion
Experiments=['Wouter21_ENZ','Wouter21_SING']

# Input Files
DirATAC = '/mnt/etemp/ahrmad/wouter/ATAC'
Dir10x = '/mnt/ndata/daniele/wouter/Processed/CellRangerArc/'
qc_ext = '_filt.h5ad'
refDir = '/mnt/etemp/ahrmad/wouter/refs'
resDir = '/mnt/etemp/ahrmad/wouter/'
All_Samples = [d for d in os.listdir(Dir10x) if d.startswith('WK')]

for Experiment in Experiments:
	if not os.path.exists(resDir+Experiment):
	    os.makedirs(resDir+Experiment)
	    print(f"Directory {resDir+Experiment} created.")
	else:
	    print(f"Directory {resDir+Experiment} already exists.")

	# Create a sample list based on the Experiment name
	Samples = []
	for sample in All_Samples:
		if '1350' in sample:
			if 'ENZ' in Experiment:
				Samples.append(sample)
		else:
			if 'SING' in Experiment:
				Samples.append(sample)

	print(f'Loading Metadata...')
	# Doublet Analysis from RNA modality is required
	# Create Whitelist Barcode Dictionary containing singlet cells by exp
	BC_dict = {}
	for sample in Samples:
		atacdf = pd.read_csv(f'{refDir}/csv/{sample}_doublet_scores_ATAC_AMU.csv')
		gexdf = pd.read_csv(f'{refDir}/csv/{sample}_doublet_scores_CB_GEX.csv')
		print(sample)
		print(f'Pre-merging: ATAC: {atacdf.shape[0]} GEX: {gexdf.shape[0]} cells')
		merged_df_all = atacdf.merge(gexdf, on='obs_names', how='inner')
		merged_df_all['doublet_class'] = 'WHATAMI'
		merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'singlet') & (merged_df_all['doublet_class_y'] == 'singlet'), 'doublet_class'] = 'Singlet Only'
		merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'doublet') & (merged_df_all['doublet_class_y'] == 'singlet'), 'doublet_class'] = 'Singlet GEX Only'
		merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'singlet') & (merged_df_all['doublet_class_y'] == 'doublet'), 'doublet_class'] = 'Singlet ATAC Only'
		merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'doublet') & (merged_df_all['doublet_class_y'] == 'doublet'), 'doublet_class'] = 'Doublet Both'
		
		print(f'Post-merging: {merged_df_all.shape[0]} cells')
	    # Retain QC passing cells (present in the CSV) 
	    # that were called singlets by at least 1 modality
	    merged_df_all_singlet = merged_df_all[merged_df_all['doublet_class'].str.contains('Singlet')]
		print(f'Singlets: {merged_df_all_singlet.shape[0]} cells')
		BC_dict[sample] = list(merged_df_all_singlet['obs_names'])

	# Load metadata
	tissueProv_dict = {}
	for sample in Samples:
		if 'AP' in sample:
			prov='AP'
		elif 'LP' in sample:
			prov='LP'
		elif 'DP' in sample:
			prov='DP'
		elif 'VP' in sample:
			prov='VP'
		else:
			raise ValueError				
		tissueProv_dict[sample] = prov

	isoMeth_dict = {}
	for sample in Samples:
		if '1350' in sample:
			iso = 'Enzymatic Digestion'
		else:
			iso = 'Singulator'
		isoMeth_dict[sample] = iso

	timePoint_dict = {}
	for sample in Samples:
		if 'I' in sample:
			timeP='Intact'
		elif 'R3' in sample:
			timeP='RegenDay3'
		elif 'Regen_Day3' in sample:
			timeP='RegenDay3'
		elif 'RegenDay2' in sample:
			timeP='RegenDay2'
		elif 'RegenDay1' in sample:
			timeP='RegenDay1'
		elif 'Day28' in sample:
			timeP='Day28'
		else:
			raise ValueError
		timePoint_dict[sample] = timeP

	mouseID_dict = {}
	for sample in Samples:
		if 'I-1' in sample:
			mID = 'WK-1350-Intact-1'
		elif 'I-2' in sample:
			mID = 'WK-1350-Intact-2'
		elif 'R3-1' in sample:
			mID = 'WK-1350-R3-1'
		elif 'R3-2' in sample:
			mID = 'WK-1350-R3-2'
		else:
			mID = sample
		mouseID_dict[sample] = mID

	# Concatenation
	#Create concatenation list
	adata_list = []
	print(f'Loading AnnData...')
	for sample in Samples:
		print(f'Reading {sample}')
		a = snap.read(os.path.join(DirATAC,sample+qc_ext)).to_memory()
		#a.obs['seqDate'] =  pd.Categorical([seqDate_dict[sample]]*a.n_obs)
		# Add metadata
		a.obs['Experiment'] = pd.Categorical([Experiment]*a.n_obs)
		a.obs['tissueProv'] = pd.Categorical([tissueProv_dict[sample]]*a.n_obs)
		a.obs['timePoint'] = pd.Categorical([timePoint_dict[sample]]*a.n_obs)
		a.obs['isoMeth'] = pd.Categorical([isoMeth_dict[sample]]*a.n_obs)
		a.obs['mouseID'] = pd.Categorical([mouseID_dict[sample]]*a.n_obs)

		# Filter cells that arent in the WL BC dict
		a = a[a.obs.index.isin(BC_dict[sample])].copy()
		# Export
		a_fname = f'{resDir}{Experiment}/{sample}_qcTOREMOVE.h5ad'
		a.write(a_fname)
		# Append to list for concatenation
		adata_list.append((sample,a_fname))
		del a

	# Exporting fragments
	print(f'Exporting fragments from all samples in {Experiment}...')
	adata = snap.AnnDataSet(adata_list,filename=f'{Experiment}_ATAC.h5ad', add_key='sample')
	adata.obs['Exp'] = [Experiment]*adata.n_obs
	snap.ex.export_fragments(adata, groupby='Exp', prefix='', suffix='.tsv.gz',out_dir=f'{resDir}{Experiment}')
	adata.close()








