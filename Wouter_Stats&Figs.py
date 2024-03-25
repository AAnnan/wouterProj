import numpy as np
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from scipy.stats import gaussian_kde
import csv
import os

path10x = '/mnt/ndata/daniele/wouter/Processed/CellRangerArc/'
exp10x = [d for d in os.listdir(path10x) if d.startswith('WK') and os.path.isdir(os.path.join(path10x, d))]
Samples = exp10x
pathATAC = '/mnt/etemp/ahrmad/wouter/ATAC'
pathGEX = '/mnt/etemp/ahrmad/wouter/RNA_CB'
qc_ext = '_qc.h5ad'
wldir=os.path.join(pathATAC, "whitelists")
ATACcsv = "_doublet_scores_ATACPeaks_Self.csv"
GEXcsv = "_doublet_scores_CB_GEX.csv"
#GEXcsv = "_doublet_scores_ATAC_AMU.csv"
refDir = '/mnt/etemp/ahrmad/wouter/refs'


#ATAC STATS gathering
#n_frag_in_cell_ATACfilt = []
#n_frag_in_cell_ATACfiltQC = []
#n_cell_filt = []
#for exp in exp10x:
#	print(exp)
#	data = snap.read(pathRes + exp +'_RAW.h5ad')
#	n_frag_in_cell_ATACfilt.append(np.sum(data.obs['n_fragment']))
#
#	data = snap.read(pathRes + exp +'_filt.h5ad')
#	n_frag_in_cell_ATACfiltQC.append(np.sum(data.obs['n_fragment']))
#	n_cell_filt.append(data.n_obs)
#ATACstats = np.vstack((exp10x,n_frag_in_cell_ATACfilt,n_frag_in_cell_ATACfiltQC,n_cell_filt))
#np.savetxt('ATACstats.csv', ATACstats, fmt='%s', delimiter=',')
at1 = np.genfromtxt(os.path.join(pathATAC, "ATACstats.csv"), delimiter=',', dtype='str')

numpy_arrays = []
for exp in exp10x:
	with open(os.path.join(wldir,exp+'_LOG'), 'r') as file:
        # Create a CSV reader for the current file
        csv_reader = csv.reader(file)
        
        # Read the CSV data and convert it into a NumPy array
        data = np.array(list(csv_reader), dtype=str)
        
        # Append the NumPy array to the list
        numpy_arrays.append(data)

lis_log = []
with open(os.path.join(pathGEX, "Wouter_scRNA_log_QC"), 'r') as file:
    # Create a CSV reader for the current file
    csv_reader = csv.reader(file)
    
    # Read the CSV data and convert it into a NumPy array
    data = list(csv_reader)
    # Append the NumPy array to the list
    lis_log.append(data)
lis_log = lis_log[0][1:]
lis_log = [el[0].split() for el in lis_log]

def extract_elements(input_list, x, y,start_line=1):
    result = []
    for i in range(start_line - 1, len(input_list), x):
        if y - 1 < len(input_list[i]):
            result.append(input_list[i][y - 1])
    return result
lis_names = [el.rstrip('.') for el in extract_elements(lis_log,14,2,2)]
log1p_total_counts_outliers = [int(el) for el in extract_elements(lis_log,14,2,5)]
log1p_n_genes_by_counts_outliers = [int(el) for el in extract_elements(lis_log,14,2,6)]
pct_counts_in_top_20_genes_outliers = [int(el) for el in extract_elements(lis_log,14,2,7)]
mt_outliers = [int(el) for el in extract_elements(lis_log,14,2,8)]
outliers = [int(el.rstrip(')').lstrip('()')) for el in extract_elements(lis_log,14,2,10)]
n_cells = [int(el) for el in extract_elements(lis_log,14,5,4)]

print(sum(log1p_total_counts_outliers)/sum(n_cells))
print(sum(log1p_n_genes_by_counts_outliers)/sum(n_cells))
print(sum(pct_counts_in_top_20_genes_outliers)/sum(n_cells))
print(sum(mt_outliers)/sum(n_cells))



n_cell_filt_atac = np.array(at1[3],dtype=np.int32)
(sum(n_cells) - sum(n_cell_filt_atac))/sum(n_cells)
names = np.array(at1[0])
n_frag_in_cell_ATACraw = np.array([arr[1][-2] for arr in numpy_arrays], dtype=np.int32)
n_frag_in_cell_GEXraw = np.array([arr[1][-1] for arr in numpy_arrays], dtype=np.int32)

n_frag_in_cell_ATACfilt = np.array(at1[1],dtype=np.int32) #atac
n_frag_in_cell_ATACfiltQC = np.array(at1[2],dtype=np.int32) #atac

n_frag_in_cell_GEXfilt = []
n_cell = []
for exp in exp10x:
	print(exp)
	adata2 = sc.read_10x_h5(os.path.join(path10x,exp+'/outs/filtered_feature_bc_matrix.h5')).to_memory()
	n_frag_in_cell_GEXfilt.append(np.sum(adata2.X))
	n_cell.append(adata2.n_obs)
n_frag_in_cell_GEXfilt = np.array(n_frag_in_cell_GEXfilt)
n_cell = np.array(n_cell)

n_frag_in_cell_GEXfiltQC = []
n_cell_filt_rna = []
for exp in exp10x:
	print(exp)
	adata = sc.read_h5ad(os.path.join(pathGEX, exp + qc_ext)).to_memory()
	n_frag_in_cell_GEXfiltQC.append(np.sum(adata.X))
	n_cell_filt_rna.append(adata.n_obs)
n_frag_in_cell_GEXfiltQC = np.array(n_frag_in_cell_GEXfiltQC)
n_cell_filt_rna = np.array(n_cell_filt_rna)


####################
# Create a histogram
plt.bar(names, n_frag_in_cell_ATACraw)
plt.xticks(rotation=45, ha='right')
# Adding labels and title
plt.xlabel('')
plt.ylabel('Total number of reads \n (raw ATAC experiment)')
plt.title('')
plt.tight_layout()
# Display the histogram
plt.savefig('ATAC_totRead_num_raw.pdf')
plt.close()
####################

####################
# Create a histogram
plt.bar(names, n_frag_in_cell_GEXraw)
plt.xticks(rotation=45, ha='right')
# Adding labels and title
plt.xlabel('')
plt.ylabel('Total number of reads \n (raw GEX experiment)')
plt.title('')
plt.tight_layout()
# Display the histogram
plt.savefig('GEX_totRead_num_raw.pdf')
plt.close()
####################

####################
# Create a histogram
plt.bar(names, (n_frag_in_cell_ATACfilt/n_frag_in_cell_ATACraw)*100)
plt.xticks(rotation=45, ha='right')
# Adding labels and title
plt.xlabel('')
plt.ylabel('Amount of reads in 10x called cells\n (% raw ATAC reads)')
plt.title('')
plt.tight_layout()
# Display the histogram
plt.savefig('ATAC_read_pct_CC.pdf')
plt.close()
####################

####################
# Create a histogram
plt.bar(names, (n_frag_in_cell_GEXfilt/n_frag_in_cell_GEXraw)*100)
plt.xticks(rotation=45, ha='right')
# Adding labels and title
plt.xlabel('')
plt.ylabel('Amount of reads in 10x called cells\n (% total GEX reads)')
plt.title('')
plt.tight_layout()
# Display the histogram
plt.savefig('GEX_read_pct_CC.pdf')
plt.close()
####################

####################
# Create a histogram
plt.bar(names, n_frag_in_cell_ATACfilt)
plt.xticks(rotation=45, ha='right')
# Adding labels and title
plt.xlabel('')
plt.ylabel('Number of reads \n (in called cells ATAC)')
plt.title('')

plt.tight_layout()
# Display the histogram
plt.savefig('ATAC_read_num_CC.pdf')
plt.close()
####################

####################
# Create a histogram
plt.bar(names, n_frag_in_cell_GEXfilt)
plt.xticks(rotation=45, ha='right')
# Adding labels and title
plt.xlabel('')
plt.ylabel('Number of reads \n (in called cells GEX)')
plt.title('')

plt.tight_layout()
# Display the histogram
plt.savefig('GEX_read_num_CC.pdf')
plt.close()
####################

####################
# Create a histogram
plt.bar(names, (n_frag_in_cell_ATACfiltQC/n_frag_in_cell_ATACfilt)*100)
plt.xticks(rotation=45, ha='right')
# Adding labels and title
plt.xlabel('')
plt.ylabel('Amount of reads remaining after QC\n (% called cells ATAC reads)')
plt.ylim(85, 100)
plt.title('')
plt.tight_layout()
# Display the histogram
plt.savefig('ATAC_read_pct_QC.pdf')
plt.close()
####################

####################
# Create a histogram
plt.bar(names, (n_frag_in_cell_GEXfiltQC/n_frag_in_cell_GEXfilt)*100)
plt.xticks(rotation=45, ha='right')
# Adding labels and title
plt.xlabel('')
plt.ylabel('Amount of reads remaining after QC\n (% called cells GEX reads)')
plt.title('')
plt.tight_layout()
# Display the histogram
plt.savefig('GEX_read_pct_QC.pdf')
plt.close()
####################

####################
# Create a histogram
plt.bar(names, n_cell)
plt.xticks(rotation=45, ha='right')
# Adding labels and title
plt.xlabel('')
plt.ylabel('Number of called cells \n (jointly called by 10x ARC)')
plt.title('')

plt.tight_layout()
# Display the histogram
plt.savefig('Joint_cell_num_CC.pdf')
plt.close()
####################

####################
# Create a histogram
plt.bar(names, (n_cell_filt_atac/n_cell)*100)
plt.xticks(rotation=45, ha='right')
# Adding labels and title
plt.xlabel('')
plt.ylabel('Amount of cells remaining after QC\n (% jointly called cells)')
plt.title('')

plt.tight_layout()
# Display the histogram
plt.savefig('ATAC_cell_pct_QC.pdf')
plt.close()
####################

####################
# Create a histogram
plt.bar(names, (n_cell_filt_rna/n_cell)*100)
plt.xticks(rotation=45, ha='right')
# Adding labels and title
plt.xlabel('')
plt.ylabel('Amount of cells remaining after QC\n (% jointly called cells)')
plt.title('')

plt.tight_layout()
# Display the histogram
plt.savefig('GEX_cell_pct_QC.pdf')
plt.close()
####################


##############################
#DOUBLET ANALYSIS
##############################

#Per Sample:
#SCATTER Score
for sample in Samples:
	print(sample)
	# Step 1: Read data from the two CSV files
	atacfile = f'{refDir}/csv/{sample}{ATACcsv}'
	gexfile = f'{refDir}/csv/{sample}{GEXcsv}'

	atacdf = pd.read_csv(atacfile)
	gexdf = pd.read_csv(gexfile)

	merged_df = atacdf.merge(gexdf, on='obs_names', how='outer').fillna(-0.05)

	# Calculate the average doublet score
	merged_df['average_doublet_score'] = (merged_df['doublet_score_x'] + merged_df['doublet_score_y']) / 2

	# Create a scatter plot with colors based on the average doublet score
	sc = plt.scatter(merged_df['doublet_score_x'], merged_df['doublet_score_y'],
	c=merged_df['average_doublet_score'],
	cmap='viridis',
	s=6,alpha=0.5)

	# Label axes and set a title
	plt.xlabel('Doublet Probability ATAC (scrublet)')
	plt.ylabel('Doublet Probability GEX (scDblFinder)')
	plt.title('Doublet Probability in ATAC and GEX')
	plt.suptitle(sample)

	## Add a color bar to indicate the average doublet scores
	cbar = plt.colorbar(sc)
	cbar.set_label('Average Doublet Probability')

	# Add a legend
	#plt.legend(loc='best')

	# Display the plot
	plt.tight_layout()

	# Save the scatter plot
	plt.savefig(f'{sample}_scatter_doublet.pdf')
	plt.close()

#SCATTER Class
for sample in Samples:
	print(sample)
	# Step 1: Read data from the two CSV files
	atacfile = f'{refDir}/csv/{sample}{ATACcsv}'
	gexfile = f'{refDir}/csv/{sample}{GEXcsv}'

	atacdf = pd.read_csv(atacfile)
	gexdf = pd.read_csv(gexfile)

	# Merge the dataframes on 'obs_names'
	merged_df = atacdf.merge(gexdf, on='obs_names', how='outer').fillna(-0.05)
	merged_df.loc[(merged_df['doublet_class_x'] == -0.05), 'doublet_class_x'] = 'singlet'
	merged_df.loc[(merged_df['doublet_class_y'] == -0.05), 'doublet_class_y'] = 'singlet'

	# Determine the 'doublet_class' in both dataframes
	merged_df['doublet_class'] = 'Singlet Both'
	merged_df.loc[(merged_df['doublet_class_x'] == 'singlet') & (merged_df['doublet_score_y'] == -0.05), 'doublet_class'] = 'Singlet 1 Mod'
	merged_df.loc[(merged_df['doublet_score_x'] == -0.05) & (merged_df['doublet_class_y'] == 'singlet'), 'doublet_class'] = 'Singlet 1 Mod'

	merged_df.loc[(merged_df['doublet_class_x'] == 'doublet') & (merged_df['doublet_class_y'] == 'doublet'), 'doublet_class'] = 'Doublet Both'
	merged_df.loc[(merged_df['doublet_class_x'] == 'doublet') & (merged_df['doublet_class_y'] == 'singlet'), 'doublet_class'] = 'Doublet ATAC Only'
	merged_df.loc[(merged_df['doublet_class_x'] == 'singlet') & (merged_df['doublet_class_y'] == 'doublet'), 'doublet_class'] = 'Doublet GEX Only'

	# Create a scatter plot with colors based on doublet class
	plt.figure(figsize=(8, 6))  # Set the figure size

	# Define colors for each class
	colors = {'Singlet 1 Mod': 'grey', 'Singlet Both': '#355070', 'Doublet Both': '#e56b6f', 'Doublet ATAC Only': '#b56576', 'Doublet GEX Only': '#6d597a'}
	colors = {'Singlet 1 Mod': 'grey', 'Singlet Both': '#2970A0', 'Doublet Both': '#FF5A63', 'Doublet ATAC Only': '#CC607B', 'Doublet GEX Only': '#7A4E8E'}

	# Map the colors based on the doublet class
	scatter_colors = merged_df['doublet_class'].map(colors)

	# Create a scatter plot with colors based on the average doublet score
	sc = plt.scatter(merged_df['doublet_score_x'], merged_df['doublet_score_y'],
	c=scatter_colors,
	s=12,alpha=0.6)

	# Label axes and set a title
	plt.xlabel('Doublet Probability ATAC (scrublet)')
	plt.ylabel('Doublet Probability GEX (scDblFinder)')
	plt.title('Doublet Classes in ATAC and GEX')
	plt.suptitle(sample)

	# Add a legend
	legend_labels = ['Singlet 1 Mod','Singlet Both', 'Doublet Both', 'Doublet ATAC Only', 'Doublet GEX Only']
	legend_handles = [plt.Line2D([0], [0], marker='o', color='w', label=label, markersize=12,alpha=0.6, markerfacecolor=colors[label]) for label in legend_labels]
	plt.legend(handles=legend_handles, bbox_to_anchor=(1.05, 1.0), loc='upper left', title='')


	# Display the plot
	plt.tight_layout()

	# Save the scatter plot
	plt.savefig(f'{sample}_scatter_doublet_class.pdf')
	plt.close()

#All Samples merged:
scatter = False
density = False
barplot = True
frames = []
for sample in Samples:
	print(sample)
	atacfile = f'{refDir}/csv/{sample}{ATACcsv}'
	gexfile = f'{refDir}/csv/{sample}{GEXcsv}'
	atacdf = pd.read_csv(atacfile)
	gexdf = pd.read_csv(gexfile)
	# Merge the dataframes on 'obs_names'
	merged_df = atacdf.merge(gexdf, on='obs_names', how='outer').fillna(-0.05)
	frames.append(merged_df)
merged_df_all = pd.concat(frames)


merged_df_all['doublet_class'] = 'QC_FILTERED'
merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'singlet') & (merged_df_all['doublet_class_y'] == 'singlet'), 'doublet_class'] = 'Singlet Both'
merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'doublet') & (merged_df_all['doublet_class_y'] == 'doublet'), 'doublet_class'] = 'Doublet Both'
merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'doublet') & (merged_df_all['doublet_class_y'] == 'singlet'), 'doublet_class'] = 'Doublet ATAC Only'
merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'singlet') & (merged_df_all['doublet_class_y'] == 'doublet'), 'doublet_class'] = 'Doublet GEX Only'

merged_df_all['average_doublet_score'] = (merged_df_all['doublet_score_x'] + merged_df_all['doublet_score_y']) / 2
merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'singlet') & (merged_df_all['doublet_score_y'] == -0.05), 'doublet_class'] = 'Singlet Unimod (ATAC)'
merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'doublet') & (merged_df_all['doublet_score_y'] == -0.05), 'doublet_class'] = 'Doublet Unimod (ATAC)'
merged_df_all.loc[(merged_df_all['doublet_score_x'] == -0.05) & (merged_df_all['doublet_class_y'] == 'singlet'), 'doublet_class'] = 'Singlet Unimod (GEX)'
merged_df_all.loc[(merged_df_all['doublet_score_x'] == -0.05) & (merged_df_all['doublet_class_y'] == 'doublet'), 'doublet_class'] = 'Doublet Unimod (GEX)'

if barplot:
	# Create a bar plot
	sns.set(style="whitegrid")
	plt.figure(figsize=(8, 6))

	# Count the occurrences of each 'doublet_class'
	plot_data = merged_df_all['doublet_class'].value_counts()

	# Create the bar plot
	#ax = sns.barplot(x=plot_data.index, y=plot_data.values, palette="Set2")
	ax = sns.barplot(x=plot_data.index, y=plot_data.values, hue=plot_data.index, palette="Set2", dodge=False)
	ax.set(xlabel=None, ylabel="Cell Count")
	plt.title("Distribution of Doublet Classes")
	plt.xticks(rotation=45)  # Rotate x-axis labels for better readability

	# Add percentages on top of the bars
	total = len(merged_df_all)
	for p in ax.patches:
	    height = p.get_height()
	    ax.annotate(f'{height/total:.2%}', (p.get_x() + p.get_width() / 2., height), ha="center")

	plt.tight_layout()
	# Save the scatter plot
	plt.savefig(f'ALL_BarPlot_Dblt_Classes.pdf')
	plt.close()
if scatter:
	# Create a scatter plot with colors based on the average doublet score
	sc = plt.scatter(merged_df_all['doublet_score_x'], merged_df_all['doublet_score_y'],
	c=merged_df_all['average_doublet_score'],
	cmap='viridis',
	s=6,alpha=0.5)

	# Label axes and set a title
	plt.xlabel('Doublet Probability ATAC (scrublet)')
	plt.ylabel('Doublet Probability GEX (scDblFinder)')
	plt.title('Doublet Probability in ATAC and GEX')
	plt.suptitle('All Samples')

	## Add a color bar to indicate the average doublet scores
	cbar = plt.colorbar(sc)
	cbar.set_label('Average Doublet Probability')

	# Add a legend
	#plt.legend(loc='best')

	# Display the plot
	plt.tight_layout()

	# Save the scatter plot
	plt.savefig(f'ALL_scatter_doublet.pdf')
	plt.close()
if density:
	# Calculate the kernel density estimate
	kde = gaussian_kde([merged_df_all['doublet_score_x'], merged_df_all['doublet_score_y']], bw_method=0.1)

	# Create a grid of points for the density estimation
	x, y = np.mgrid[merged_df_all['doublet_score_x'].min():merged_df_all['doublet_score_x'].max():100j, merged_df_all['doublet_score_y'].min():merged_df_all['doublet_score_y'].max():100j]
	positions = np.vstack([x.ravel(), y.ravel()])
	t = np.reshape(kde(positions).T, x.shape)
	#z = np.log2(t)
	z = np.sqrt(t)

	plt.figure(figsize=(8, 6))
	plt.contourf(x, y, z, cmap='plasma', levels=30, antialiased=True)  # Create a filled contour plot for density
	plt.ylim(0,1)
	plt.xlim(0,1)
	# Create the scatterplot on top of the density plot
	#plt.scatter(merged_df_all['doublet_score_x'], merged_df_all['doublet_score_y'], c='k', marker='o', s=100, edgecolors='k')

	# Add color bar
	cbar = plt.colorbar()
	cbar.set_label('Density', rotation=90)

	# Label axes and set a title
	plt.xlabel('Doublet Probability ATAC (scrublet)')
	plt.ylabel('Doublet Probability GEX (scDblFinder)')
	plt.title('Doublet Probability in ATAC and GEX')
	plt.suptitle('All Samples')
	# Display the plot
	plt.tight_layout()
	# Save the scatter plot
	plt.savefig(f'ALL_scatter_doublet_density.pdf')
	plt.close()

##############################
#VENN DIAGRAMS QC
##############################

enz = [0,0,0]
sing = [0,0,0]
for sample in Samples:
	print(sample)
	# Step 1: Read data from the two CSV files
	atacfile = f'{refDir}/csv/{sample}{ATACcsv}'
	gexfile = f'{refDir}/csv/{sample}{GEXcsv}'

	atacdf = pd.read_csv(atacfile)
	gexdf = pd.read_csv(gexfile)

	atacset = set(atacdf['obs_names'])
	gexset = set(gexdf['obs_names'])

	# Elements only in set1
	only_in_atacset = list(atacset - gexset)
	# Elements only in set2
	only_in_gexset = list(gexset - atacset)
	# Elements in both sets
	in_both = list(atacset.intersection(gexset))

	if '1350' in sample:
		enz = [i + j for i, j in zip(enz, [len(only_in_atacset),len(only_in_gexset),len(in_both)])]
	else:
		sing = [i + j for i, j in zip(sing, [len(only_in_atacset),len(only_in_gexset),len(in_both)])]
	
	#Create a Venn diagram
	plt.figure()
	total_count = len(only_in_atacset) + len(only_in_gexset) + len(in_both)
	venn2([atacset, gexset], ('ATAC', 'GEX'),subset_label_formatter=lambda x: f"{x} ({(x/total_count):1.0%})")
	plt.title("QC-passing cells in ATAC and GEX")
	plt.suptitle(sample)
	plt.tight_layout()
	# Display the histogram
	plt.savefig(f'{sample}_venn_cells.pdf')
	plt.close()

#ENZYMATIC DIGESTION VENNS
plt.figure()
total_count = sum(enz)
# Create a Venn diagram with the numbers
venn2(subsets=(enz[0], enz[1], enz[2]), set_labels=('ATAC', 'GEX'),subset_label_formatter=lambda x: f"{x} ({(x/total_count):1.0%})")
plt.title("QC-passing cells in ATAC and GEX\nEnzymatic Digestion")
plt.tight_layout()
# Display the histogram
plt.savefig(f'ENZ_venn_cells.pdf')
plt.close()

#SINGULATOR VENNS
plt.figure()
# Create a Venn diagram with the numbers
total_count = sum(sing)
venn2(subsets=(sing[0], sing[1], sing[2]), set_labels=('ATAC', 'GEX'),subset_label_formatter=lambda x: f"{x} ({(x/total_count):1.0%})")
plt.title("QC-passing cells in ATAC and GEX\nSingulator")

plt.tight_layout()
# Display the histogram
plt.savefig(f'SING_venn_cells.pdf')
plt.close()

#ALL VENNS
plt.figure()
# Create a Venn diagram with the numbers
total_count = sum(sing) + sum(enz)
venn2(subsets=(sing[0]+enz[0], sing[1]+enz[1], sing[2]+enz[2]), set_labels=('ATAC', 'GEX'),subset_label_formatter=lambda x: f"{x} ({(x/total_count):1.0%})")
plt.title("QC-passing cells in ATAC and GEX")

plt.tight_layout()
# Display the histogram
plt.savefig(f'ALL_venn_cells.pdf')
plt.close()

# Isolation Method BARPLOTs
frames = []
ENZSamples = [sample for sample in Samples if '1350' in sample]
for sample in ENZSamples:
	print(sample)
	atacfile = f'{refDir}/csv/{sample}{ATACcsv}'
	gexfile = f'{refDir}/csv/{sample}{GEXcsv}'
	atacdf = pd.read_csv(atacfile)
	gexdf = pd.read_csv(gexfile)
	# Merge the dataframes on 'obs_names'
	merged_df = atacdf.merge(gexdf, on='obs_names', how='outer').fillna(-0.05)
	frames.append(merged_df)
merged_df_all = pd.concat(frames)


merged_df_all['doublet_class'] = 'QC_FILTERED'
merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'singlet') & (merged_df_all['doublet_class_y'] == 'singlet'), 'doublet_class'] = 'Singlet Both'
merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'doublet') & (merged_df_all['doublet_class_y'] == 'doublet'), 'doublet_class'] = 'Doublet Both'
merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'doublet') & (merged_df_all['doublet_class_y'] == 'singlet'), 'doublet_class'] = 'Doublet ATAC Only'
merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'singlet') & (merged_df_all['doublet_class_y'] == 'doublet'), 'doublet_class'] = 'Doublet GEX Only'

# Create a bar plot
sns.set(style="whitegrid")
plt.figure(figsize=(8, 6))
#ENZdf = merged_df_all.loc[(merged_df_all['sample_x'].str.contains('1350')) | (merged_df_all['sample_y'].str.contains('1350'))]
ENZdf = merged_df_all.copy()
total = len(ENZdf)
# Count the occurrences of each 'doublet_class'
plot_data = ENZdf['doublet_class'].value_counts()

# Create the bar plot
#ax = sns.barplot(x=plot_data.index, y=plot_data.values, palette="Set2")
ax = sns.barplot(x=plot_data.index, y=plot_data.values, hue=plot_data.index, palette="Set2", dodge=False)
ax.set(xlabel=None, ylabel="Cell Count")
plt.title("Distribution of Doublet Classes")
plt.xticks(rotation=45)  # Rotate x-axis labels for better readability

# Add percentages on top of the bars
for p in ax.patches:
    height = p.get_height()
    ax.annotate(f'{height/total:.2%}', (p.get_x() + p.get_width() / 2., height), ha="center")

plt.tight_layout()
# Save the scatter plot
plt.savefig(f'ENZ_BarPlot_Dblt_Classes.pdf')
plt.close()

frames = []
SingSamples = [sample for sample in Samples if '1350' not in sample]
for sample in SingSamples:
	print(sample)
	atacfile = f'{refDir}/csv/{sample}{ATACcsv}'
	gexfile = f'{refDir}/csv/{sample}{GEXcsv}'
	atacdf = pd.read_csv(atacfile)
	gexdf = pd.read_csv(gexfile)
	# Merge the dataframes on 'obs_names'
	merged_df = atacdf.merge(gexdf, on='obs_names', how='outer').fillna(-0.05)
	frames.append(merged_df)
merged_df_all = pd.concat(frames)


merged_df_all['doublet_class'] = 'QC_FILTERED'
merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'singlet') & (merged_df_all['doublet_class_y'] == 'singlet'), 'doublet_class'] = 'Singlet Both'
merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'doublet') & (merged_df_all['doublet_class_y'] == 'doublet'), 'doublet_class'] = 'Doublet Both'
merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'doublet') & (merged_df_all['doublet_class_y'] == 'singlet'), 'doublet_class'] = 'Doublet ATAC Only'
merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'singlet') & (merged_df_all['doublet_class_y'] == 'doublet'), 'doublet_class'] = 'Doublet GEX Only'

	# Create a bar plot
sns.set(style="whitegrid")
plt.figure(figsize=(8, 6))
#SINGdf = merged_df_all.loc[~((merged_df_all['sample_x'].str.contains('1350')) | (merged_df_all['sample_y'].str.contains('1350')))]
SINGdf = merged_df_all.copy()
total = len(SINGdf)
# Count the occurrences of each 'doublet_class'
plot_data = SINGdf['doublet_class'].value_counts()

# Create the bar plot
#ax = sns.barplot(x=plot_data.index, y=plot_data.values, palette="Set2")
ax = sns.barplot(x=plot_data.index, y=plot_data.values, hue=plot_data.index, palette="Set2", dodge=False)
ax.set(xlabel=None, ylabel="Cell Count")
plt.title("Distribution of Doublet Classes")
plt.xticks(rotation=45)  # Rotate x-axis labels for better readability

# Add percentages on top of the bars

for p in ax.patches:
    height = p.get_height()
    ax.annotate(f'{height/total:.2%}', (p.get_x() + p.get_width() / 2., height), ha="center")

plt.tight_layout()
# Save the scatter plot
plt.savefig(f'SING_BarPlot_Dblt_Classes.pdf')
plt.close()

#ALL VENN SINGLET
from collections import Counter
merged_df_all['doublet_class'].value_counts()
cnt = Counter([el for el in merged_df_all['doublet_class'] if 'inglet' in el])
total_count = sum(cnt.values())
plt.figure()
# Create a Venn diagram with the numbers
venn2(subsets=(cnt['Singlet Unimod (ATAC)'],cnt['Singlet Unimod (GEX)'],cnt['Singlet Both']), set_labels=('ATAC','GEX'),subset_label_formatter=lambda x: f"{x} ({(x/total_count):1.0%})")
plt.title("Singlets in ATAC and GEX")

plt.tight_layout()
# Display the histogram
plt.savefig(f'ALL_venn_singlet.pdf')
plt.close()

#ALL VENN DOUBLET
from collections import Counter
merged_df_all['doublet_class'].value_counts()
cnt = Counter([el for el in merged_df_all['doublet_class'] if 'oublet' in el if 'Unimod' not in el])
total_count = sum(cnt.values())
plt.figure()
# Create a Venn diagram with the numbers
venn2(subsets=(cnt['Doublet ATAC Only'],cnt['Doublet GEX Only'],cnt['Doublet Both']), set_labels=('ATAC','GEX'),subset_label_formatter=lambda x: f"{x} ({(x/total_count):1.0%})")
plt.title("Doublets in ATAC and GEX")

plt.tight_layout()
# Display the histogram
plt.savefig(f'ALL_venn_doublet.pdf')
plt.close()

#### CellBender stats

# Cell Calling Hist
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

#gexdf = pd.read_csv(gexfile)
stats_dict = {}

CB_stats = np.unique([d for d in os.listdir('.') if d.startswith('WK') if d.endswith('_CellBender_metrics.csv')])

for CB_stat in CB_stats:
	stat_df = pd.read_csv(CB_stat,header=None)

	batch = CB_stat.split('_CellBender_metrics.csv')[0]
	val = [int(stat_df[stat_df[0] == 'expected_cells'][1].iloc[0]),int(stat_df[stat_df[0] == 'found_cells'][1].iloc[0])]
	stats_dict[batch] = val

# Extract keys, first values, and second values from the dictionary
keys = list(stats_dict.keys())
value_1 = [item[0] for item in stats_dict.values()]
value_2 = [item[1] for item in stats_dict.values()]

# Plotting the histogram
bar_width = 0.35
index = np.arange(len(keys))
fig, ax = plt.subplots()

# Plot the first values (10x) in blue
bar1 = ax.bar(index, value_1, bar_width, label='10x Joint Cell Calling', color='#fc8d59',alpha=0.8)

# Plot the second values (CellBender) in orange
bar2 = ax.bar(index + bar_width, value_2, bar_width, label='CellBender Cell Calling on RNA', color='#91bfdb',alpha=0.8)

# Customize the plot
ax.set_ylabel('Cell Count')
ax.set_xticks(index + bar_width / 2)
ax.set_xticklabels(keys, rotation=45, ha='right')
ax.legend()


# Show the plot
plt.tight_layout()
plt.savefig('CB_Cell_Calling_Comp_Hist.pdf')
plt.close()


# Cell Calling Venn
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles
import numpy as np
import pandas as pd
import os

CB_stats = np.unique([d.split('_CellBender_metrics.csv')[0] for d in os.listdir('.') if d.startswith('WK') if d.endswith('_CellBender_metrics.csv')])

for CB_stat in CB_stats:
	print(f'Venn {CB_stat}')
	CB_bcs = pd.Series(pd.read_csv(CB_stat+"_CellBender_cell_barcodes.csv",header=None)[0], name='CellBender')
	ten_bcs = pd.read_csv('raws/'+CB_stat+"_per_barcode_metrics.csv",usecols=['barcode','is_cell'])
	ten_bcs = pd.Series(ten_bcs[ten_bcs['is_cell']==1]['barcode'],name="10x")

	# Create a Venn diagram
	venn_labels = {'14564360': 'CellBender', '246234461': '10x', '1436234641': 'Both'}
	venn2(subsets=(set(CB_bcs), set(ten_bcs)), set_labels=('CellBender', '10x'), ax=plt.gca(), set_colors=('skyblue', 'green'))
	venn2_circles(subsets=(set(CB_bcs), set(ten_bcs)), ax=plt.gca(), linewidth=1.5)

	# Add labels
	for text in plt.gca().texts:
	    if text.get_text() in venn_labels:
	        text.set_text(venn_labels[text.get_text()])

	# Show the plot
	plt.title(CB_stat)
	plt.tight_layout()
	plt.savefig(f'{CB_stat}_Venn_BCs.pdf')
	plt.close()

#pdfjam --nup 4x3 --landscape --a4paper *_Venn*.pdf --outfile CB_Cell_Calling_Comp_Venn.pdf

# Read drop Hist
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

#gexdf = pd.read_csv(gexfile)
stats_dict = {}

CB_stats = np.unique([d for d in os.listdir('.') if d.startswith('WK') if d.endswith('_CellBender_metrics.csv')])

for CB_stat in CB_stats:
	stat_df = pd.read_csv(CB_stat,header=None)

	batch = CB_stat.split('_CellBender_metrics.csv')[0]
	val = [int(stat_df[stat_df[0] == 'total_raw_counts'][1].iloc[0]),int(stat_df[stat_df[0] == 'total_output_counts'][1].iloc[0])]
	stats_dict[batch] = val

# Extract keys, first values, and second values from the dictionary
keys = list(stats_dict.keys())
value_1 = [item[0] for item in stats_dict.values()]
value_2 = [item[1] for item in stats_dict.values()]

# Plotting the histogram
bar_width = 0.35
index = np.arange(len(keys))
fig, ax = plt.subplots()

# Plot the first values (10x) in blue
bar1 = ax.bar(index, value_1, bar_width, label='Total Raw Counts', color='#fc8d59',alpha=0.8)

# Plot the second values (CellBender) in orange
bar2 = ax.bar(index + bar_width, value_2, bar_width, label='CellBender Counts', color='#91bfdb',alpha=0.8)

# Customize the plot
ax.set_ylabel('Read Count')
ax.set_xticks(index + bar_width / 2)
ax.set_xticklabels(keys, rotation=45, ha='right')
ax.legend()


# Show the plot
plt.tight_layout()
plt.savefig('CB_Read_Drop_Hist.pdf')
plt.close()

# Read drop ratio Hist
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

stats_dict = {}

CB_stats = np.unique([d for d in os.listdir('.') if d.startswith('WK') if d.endswith('_CellBender_metrics.csv')])

for CB_stat in CB_stats:
	stat_df = pd.read_csv(CB_stat,header=None)

	batch = CB_stat.split('_CellBender_metrics.csv')[0]
	val = float(stat_df[stat_df[0] == 'fraction_counts_removed'][1].iloc[0]*100)
	stats_dict[batch] = val

# Extract keys, first values, and second values from the dictionary
keys = list(stats_dict.keys())
value_1 = stats_dict.values()

# Plotting the histogram
#bar_width = 0.35
index = np.arange(len(keys))
fig, ax = plt.subplots()

# Plot the first values (10x) in blue
bar1 = ax.bar(index, value_1, bar_width, color='#91bfdb',alpha=0.8)

# Customize the plot
ax.set_ylabel('Fraction Reads removed (%)')
ax.set_xticks(index)
ax.set_xticklabels(keys, rotation=45, ha='right')


# Show the plot
plt.tight_layout()
plt.savefig('CB_Read_Drop_Ratio_Hist.pdf')
plt.close()
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import snapatac2 as snap
import pandas as pd
import numpy as np
import csv
import os

Path10x = '/mnt/ndata/daniele/wouter/Processed/CellRangerArc/'
PathGEX = '/mnt/etemp/ahrmad/wouter/batch_scanpy_FiltTop'
PathATAC = '/mnt/etemp/ahrmad/wouter/snap'
Samples = [d for d in os.listdir(Path10x) if d.startswith('WK') and os.path.isdir(os.path.join(Path10x, d))]

ext_GEX = '_qcFiltTop.h5ad'
ext_ATAC = '_filt.h5ad'

lis_log = []
with open(os.path.join(PathGEX, "Wouter_scRNA_log_QC"), 'r') as file:
    # Create a CSV reader for the current file
    csv_reader = csv.reader(file)
    
    # Read the CSV data and convert it into a NumPy array
    data = list(csv_reader)
    # Append the NumPy array to the list
    lis_log.append(data)
lis_log = lis_log[0][1:]
lis_log = [el[0].split() for el in lis_log]

def extract_elements(input_list, x, y,start_line=1):
    result = []
    for i in range(start_line - 1, len(input_list), x):
        if y - 1 < len(input_list[i]):
            result.append(input_list[i][y - 1])
    return result
n_cells = [int(el) for el in extract_elements(lis_log,14,5,4)]
sample_names = [el.rstrip('.') for el in extract_elements(lis_log,14,2,2)]

cell_info = {}
cell_info['names'] = sample_names
cell_info['n_cell'] = n_cells
cell_info['n_cell_RNAQC'] = []
cell_info['n_cell_ATACQC'] = []
cell_info['n_cell_BOTHQC'] = []
cell_info['n_cell_Singlets'] = []

for sample in Samples:
	print(sample)
	# Step 1: Read data from the two CSV files
	atacfile = f'{refDir}/csv/{sample}{ATACcsv}'
	gexfile = f'{refDir}/csv/{sample}{GEXcsv}'

	atacdf = pd.read_csv(atacfile)
	gexdf = pd.read_csv(gexfile)

	atacset = set(atacdf['obs_names'])
	gexset = set(gexdf['obs_names'])

	cell_info['n_cell_RNAQC'].append(len(gexset))
	cell_info['n_cell_ATACQC'].append(len(atacset))

	# Elements only in set1
	only_in_atacset = list(atacset - gexset)
	# Elements only in set2
	only_in_gexset = list(gexset - atacset)
	# Elements in both sets
	in_both = list(atacset.intersection(gexset))
	cell_info['n_cell_BOTHQC'].append(len(in_both))

	# Merge the dataframes on 'obs_names'
	merged_df_all = atacdf.merge(gexdf, on='obs_names', how='outer').fillna(-0.05)
	merged_df_all['doublet_class'] = 'WHAAT'
	merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'singlet') & (merged_df_all['doublet_class_y'] == 'singlet'), 'doublet_class'] = 'Singlet Only'
	merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'singlet') & (merged_df_all['doublet_score_y'] == -0.05), 'doublet_class'] = 'Singlet Unimod (ATAC)'
	merged_df_all.loc[(merged_df_all['doublet_score_x'] == -0.05) & (merged_df_all['doublet_class_y'] == 'singlet'), 'doublet_class'] = 'Singlet Unimod (GEX)'
	merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'doublet') & (merged_df_all['doublet_score_y'] == -0.05), 'doublet_class'] = 'Doublet Unimod (ATAC)'
	merged_df_all.loc[(merged_df_all['doublet_score_x'] == -0.05) & (merged_df_all['doublet_class_y'] == 'doublet'), 'doublet_class'] = 'Doublet Unimod (GEX)'
	merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'doublet') & (merged_df_all['doublet_class_y'] == 'doublet'), 'doublet_class'] = 'Doublet Both'
	merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'doublet') & (merged_df_all['doublet_class_y'] == 'singlet'), 'doublet_class'] = 'Doublet ATAC Only'
	merged_df_all.loc[(merged_df_all['doublet_class_x'] == 'singlet') & (merged_df_all['doublet_class_y'] == 'doublet'), 'doublet_class'] = 'Doublet GEX Only'
	cell_info['n_cell_Singlets'].append(np.sum(merged_df_all['doublet_class'].str.contains('Only')))

# Create an array for the x-axis positions
x = np.arange(len(cell_info['names']))

# Set the width of the bars
bar_width = 0.18

# Create a figure and axis
fig, ax = plt.subplots(figsize=(10, 6))

bar1 = ax.bar(x - 2*bar_width, cell_info['n_cell'], width=bar_width, label='Total Cells', color='#66c2a5')
bar3 = ax.bar(x - bar_width, cell_info['n_cell_ATACQC'], width=bar_width, label='Passed ATAC QC', color='#8da0cb')
bar2 = ax.bar(x, cell_info['n_cell_RNAQC'], width=bar_width, label='Passed GEX QC', color='#fc8d62')
bar4 = ax.bar(x + bar_width, cell_info['n_cell_BOTHQC'], width=bar_width, label='Passed Both QC', color='#e78ac3')
bar5 = ax.bar(x + 2*bar_width, cell_info['n_cell_Singlets'], width=bar_width, label='QC-passing Singlets', color='#a6d854')

#hist
#bar1 = ax.bar(x - 2*bar_width, cell_info['n_cell'], width=bar_width, label='Total Cells')
#bar3 = ax.bar(x - bar_width, cell_info['n_cell_ATACQC'], width=bar_width, label='Passed ATAC QC')
#bar2 = ax.bar(x, cell_info['n_cell_RNAQC'], width=bar_width, label='Passed GEX QC')
#bar4 = ax.bar(x + bar_width, cell_info['n_cell_BOTHQC'], width=bar_width, label='Passed Both QC')
#bar5 = ax.bar(x + 2*bar_width, cell_info['n_cell_Singlets'], width=bar_width, label='QC-passing Singlets')

# Customize the plot
ax.set_xlabel('')
ax.set_ylabel('Number of Cells', fontsize=14)
ax.set_title('Quality Control Summary', fontsize=16)
ax.set_xticks(x)
ax.set_xticklabels(cell_info['names'], rotation=45, ha='right', fontsize=12)
ax.legend(fontsize=12)

plt.xticks(rotation=45, ha='right')
plt.tight_layout()
# Display the histogram
plt.savefig('QC_summary.pdf')
plt.close()
####################




