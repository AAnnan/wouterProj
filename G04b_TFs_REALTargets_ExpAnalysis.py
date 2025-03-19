import pandas as pd
import numpy as np
import decoupler as dc
import scanpy as sc
import matplotlib.pyplot as plt

adata_dir = '/mnt/etemp/ahrmad/wouter/RNA_CB/annot_other_norm/'

def make_target_genes_tsv_from_tf_and_net(tf_of_interest, net):
    # Convert the transcription factor to title case for consistent filtering
    tf_of_interest = tf_of_interest.title()
    
    # Filter 'net' DataFrame rows where 'source' matches the TF of interest
    filtered_df = net[net['source'].str.contains(tf_of_interest, case=False)]
    
    # Select 'target' and 'mor' columns only and rename them with the TF name
    filtered_df = filtered_df[['target', 'weight']]
    filtered_df.columns = [tf_of_interest, f"{tf_of_interest}_weight"]
    
    # Check if the filtered dataframe is empty
    if filtered_df.empty:
        print(f"No matches found for transcription factor: {tf_of_interest}")
        return
    
    # Write the filtered dataframe to a TSV file named after the TF of interest
    output_filename = f"{tf_of_interest}_Targets.tsv"
    filtered_df.to_csv(output_filename, sep='\t', index=False, quoting=0)
    
    # Print the resulting dataframe
    print(filtered_df)
    return 0


net = dc.get_collectri(organism='mouse', split_complexes=False)

# Read the list of TFs from a file
tf_of_interest_list = pd.read_csv('CleanTFs_Hv12_L1_I-CD28.tsv', header=None)[0].tolist()

# Read the list of TFs from filtered_genes
#pd.DataFrame(filtered_genes['Gene'], columns=['Gene']).to_csv('CleanTFs_Hv12_L1_I-CD28.tsv', sep='\t', index=False,header=False)

# Loop through each TF and call the make_target_genes_tsv_from_tf_and_net function
for tf in tf_of_interest_list:
    make_target_genes_tsv_from_tf_and_net(tf, net)

def fast_dotplot(adat_f,TF_targets_f,TF_ana_f,cat_order_f,suf_f,TF_Source,Experiment,positions_f,values_f):

        if 'Rel' in suf_f:
            std_scale = 'var'
        elif 'Abs' in suf_f:
            std_scale = None

        info = TF_ana_f.split('_')
        sc.pl.DotPlot.DEFAULT_LEGENDS_WIDTH = 2.2  # inches

        sc.pl.DotPlot(adat_f, var_names=TF_targets_f, groupby='timePoint', use_raw=False,
            mean_only_expressed=False, expression_cutoff=0,standard_scale=std_scale,#'var', #defaults stuff
            categories_order=cat_order_f,cmap='bwr',
            var_group_positions=positions_f,var_group_labels=values_f,
            title=f'{Experiment}\n{info[1]} vs {info[2]} - {TF_Source} TARGETS').savefig(f'./figures/{TF_Source}_TARGETS_{Experiment}_{TF_ana_f}_{suf_f}.pdf')
        plt.close()
        return 0


Experiment2020='Wouter2020'
Experiment21='Wouter21_SING'
#all2020
#cat_order_all2020 = pd.Categorical(['intact1','intact2','CastDay1', 'CastDay7', 'CastDay14', 'CastDay28', 'RegenDay1', 'RegenDay2', 'RegenDay3','RegenDay7', 'RegenDay14', 'RegenDay28', 'Epi1', 'Epi2','NonEpi1', 'NonEpi2', 'Unsorted1'],categories=['intact1','intact2','CastDay1', 'CastDay7', 'CastDay14', 'CastDay28', 'RegenDay1', 'RegenDay2', 'RegenDay3','RegenDay7', 'RegenDay14', 'RegenDay28', 'Epi1', 'Epi2','NonEpi1', 'NonEpi2', 'Unsorted1'])
cat_order_all2020 = pd.Categorical(['intact1','intact2','CastDay1', 'CastDay7', 'CastDay14', 'CastDay28', 'RegenDay1', 'RegenDay2', 'RegenDay3','RegenDay7', 'RegenDay14', 'RegenDay28'],categories=['intact1','intact2','CastDay1', 'CastDay7', 'CastDay14', 'CastDay28', 'RegenDay1', 'RegenDay2', 'RegenDay3','RegenDay7', 'RegenDay14', 'RegenDay28'])
cat_order_all2020 = pd.Categorical(['intact1','intact2','CastDay28','RegenDay3'],categories=['intact1','intact2','CastDay28','RegenDay3'])
cat_order_all21 = pd.Categorical([ 'Intact','Day28', 'RegenDay1', 'RegenDay2', 'RegenDay3'], categories=['Intact','Day28', 'RegenDay1', 'RegenDay2', 'RegenDay3'])
cat_order_all21 = pd.Categorical([ 'Intact','Day28', 'RegenDay3'], categories=['Intact','Day28', 'RegenDay3'])

adata2020 = sc.read_h5ad(adata_dir + 'Wouter2020.h5ad')
adata2020 = adata2020[np.isin(adata2020.obs['timePoint'],cat_order_all2020)].copy()
adata2020.obs['timePoint'].cat.remove_unused_categories()

adata21 = sc.read_h5ad(adata_dir + 'Wouter21_SING_CB_annot_Lum_newSigs.h5ad')
adata_L1 = adata21[np.isin(adata21.obs['Annotation'],['L1'])].copy()
adata_L1 = adata_L1[np.isin(adata_L1.obs['timePoint'],cat_order_all21)].copy()
adata_L1.obs['timePoint'].cat.remove_unused_categories()

TF_ana = 'L1_I_CD28'
tf_of_interest_list = pd.read_csv('CleanTFs_Hv12_L1_I-CD28.tsv', header=None)[0].tolist()

too_many_genes = []
for tf_process in tf_of_interest_list:
    print(f"\n\nPROCESSING {tf_process} TARGETS\n\n")
    try:
        TF_targets = pd.read_csv(f"{tf_process}_Targets.tsv",sep = '\t')[tf_process].tolist()
    except FileNotFoundError:
        # Handle the error if the file does not exist
        print(f"The file was not found. The TF {tf_process} had no found target in OmniPath")
        continue
    else:
        TF_targets_in21 = set([g for g in TF_targets if g in adata_L1.var.index])
        print(f'Lost {set(TF_targets) - set(TF_targets_in21)}: NOT in Wouter21')
        TF_targets_in20 = set([g for g in TF_targets if g in adata2020.var.index])
        print(f'Lost {set(TF_targets) - set(TF_targets_in20)}: NOT in Wouter2020')
        TF_targets = list(TF_targets_in21.intersection(TF_targets_in20))

        if len(TF_targets) == 0:
            print(f"The TF {tf_process} was found, but none of the targets found in ADATAs")
            continue

        elif len(TF_targets) > 130:
            print(f"Too many genes regulated to plot {tf_process} regulated genes dotplot")
            too_many_genes.append(tf_process)
            continue
        else:
            TF_targets_df = pd.read_csv(f"{tf_process}_Targets.tsv",sep = '\t')
            TF_targets_df = TF_targets_df[TF_targets_df[tf_process].isin(TF_targets)]
            df_sorted = TF_targets_df.sort_values(by=tf_process+'_weight').reset_index(drop=True)
            TF_name = tf_process
            
            # Initializing lists to store positions and corresponding values
            positions = []
            values = []

            # Identifying consecutive same values
            start_idx = 0

            for i in range(1, len(df_sorted)):
                if df_sorted[f'{TF_name}_weight'][i] != df_sorted[f'{TF_name}_weight'][start_idx]:
                    if i - start_idx > 1:  # There are consecutive values
                        positions.append((start_idx, i - 1))
                        values.append(round(df_sorted[f'{TF_name}_weight'][start_idx], 1))
                    start_idx = i

            # Checking the last group of consecutive values
            if len(df_sorted) - start_idx > 1:
                positions.append((start_idx, len(df_sorted) - 1))
                values.append(round(df_sorted[f'{TF_name}_weight'][start_idx], 1))

            # Handling single entries (not consecutive but unique values)
            for i in range(len(df_sorted)):
                if i == 0 or df_sorted[f'{TF_name}_weight'][i] != df_sorted[f'{TF_name}_weight'][i - 1]:
                    if i == len(df_sorted) - 1 or df_sorted[f'{TF_name}_weight'][i] != df_sorted[f'{TF_name}_weight'][i + 1]:
                        positions.append((i, i))
                        values.append(round(df_sorted[f'{TF_name}_weight'][i], 1))

            # Sorting positions and values to ensure correct order after adding single entries
            sorted_positions_values = sorted(zip(positions, values))
            positions, values = zip(*sorted_positions_values)

            positions = list(positions)
            values = [str(va) for va in values]

            fast_dotplot(adata2020,TF_targets,TF_ana,cat_order_all2020,'Abs_AllTP',tf_process,Experiment2020,positions,values) 
            fast_dotplot(adata2020,TF_targets,TF_ana,cat_order_all2020,'Rel_AllTP',tf_process,Experiment2020,positions,values)
            fast_dotplot(adata_L1,TF_targets,TF_ana,cat_order_all21,'Abs_AllTP',tf_process,Experiment21,positions,values)
            fast_dotplot(adata_L1,TF_targets,TF_ana,cat_order_all21,'Rel_AllTP',tf_process,Experiment21,positions,values)

print(f'Those TFs had too many targets (>130 genes) to be plotted:{too_many_genes}')

