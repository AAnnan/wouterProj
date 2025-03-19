##Downstream Analysis

# Purpose:
#   define open chromatin and TF’s binding motifs enriched
#   enriched in one vs all 

#   for example in the singular: enriched in D28/RD3 vs INTACT 
#   or pairwise comparision (castrationD28 vs intact and RD3 vs intact etc)

import snapatac2 as snap
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt
import os
import warnings
from scipy.stats import ttest_1samp,ttest_ind, wilcoxon
from collections import Counter
warnings.simplefilter(action='ignore', category=FutureWarning)
sc.settings.verbosity = 0
sc.settings.set_figure_params(
	figsize=(6, 6),
    dpi_save=300,
    fontsize=12,
    facecolor="white",
    frameon=True
    ,
)
print(f'Anndata: {ad.__version__}\nSnapatac2: {snap.__version__}\nScanpy: {sc.__version__}')

# To change according to samples
# SING for singulator, ENZ for enzymatic digestion
Experiment='Wouter21_SING'
temp_macs3 = '/mnt/etemp/ahrmad/wouter/Wouter21_SING/annot/temp'

atac_h5ad = '/mnt/etemp/ahrmad/wouter/Wouter21_SING/annot/Wouter21_SING_Post.h5ad'
rna_h5ad = "/mnt/etemp/ahrmad/wouter/RNA_CB/annot_other_norm/Wouter21_SING_CB_annot_All_newSigs.h5ad"


###########
# FUNCTIONS
###########

def output_TFs_tsv(motifs, log_fc_min, max_adjPval,output_file):
    # Apply the conditions to each DataFrame in the dictionary
    filtered_motifs = {
        key: df.filter(
            (
                (df["log2(fold change)"] >= log_fc_min) | 
                (df["log2(fold change)"] <= -log_fc_min)
            ) & 
            (df["adjusted p-value"] < max_adjPval)
        ).with_columns([pl.lit(key).alias("provenance")])  # Add provenance column
        for key, df in motifs.items()
    }

    # Concatenate all filtered DataFrames
    combined_df = pl.concat(
        [
            df.select(["id", "log2(fold change)", "adjusted p-value", "provenance"])
            for df in filtered_motifs.values()
        ]
    )

    # Pivot the DataFrame to have separate columns for each provenance
    pivot_df = combined_df.pivot(
        values=["log2(fold change)", "adjusted p-value"], 
        index="id", 
        on="provenance"
    )

    # Save the result to a CSV file
    pivot_df.write_csv(f'./figures/{output_file}',separator='\t')
    
    return pivot_df.shape[0]

# Volcano Plot
def plot_volcano(data,name="VolcanoPlot_DiffTest",fc_threshold = 1,p_threshold = 0.01,whats_tested='Peaks',pval_lim=0,lege=True):
    # Extracting necessary columns
    log2_fc = data["log2(fold_change)"].to_numpy()
    adj_p_not = data["adjusted p-value"].to_numpy()
    adj_p = adj_p_not.copy()

    adj_p[adj_p < pval_lim] = pval_lim

    # Compute -log10 of the adjusted p-value
    neg_log10_p = -np.log10(adj_p)

    # Colors
    colors = np.where((log2_fc > fc_threshold) & (adj_p < p_threshold), "red", 
                      np.where((log2_fc < -fc_threshold) & (adj_p < p_threshold), "blue", "gray"))

    # Create the plot
    sns.set(style="whitegrid")
    plt.figure(figsize=(10, 8))
    
    # Scatter plot
    scatter = plt.scatter(log2_fc, neg_log10_p, c=colors, alpha=0.7, linewidths=0.1, edgecolor='k', s=40)

    # Add thresholds for fold change and significance
    plt.axvline(x=fc_threshold, color="red", linestyle="--", linewidth=1)
    plt.axvline(x=-fc_threshold, color="blue", linestyle="--", linewidth=1)
    plt.axhline(y=-np.log10(p_threshold), color="green", linestyle="--", linewidth=1)

    # Annotate top points
    #top_features = data.sort("adjusted p-value").head(10)
    #for _, row in top_features.to_pandas().iterrows():
    #    plt.text(row[1], -np.log10(row[3]), row[0], fontsize=8, ha='right')

    # Labels and Title
    plt.title(f"{wanted_tp_comp[1]} vs {wanted_tp_comp[0]}\n{data.shape[0]} {whats_tested} tested", fontsize=16)
    plt.xlabel("Log2(Fold Change)", fontsize=14)
    plt.ylabel("-Log10(Adjusted P-Value)", fontsize=14)

    if lege:
        # Custom Legend
        plt.legend(handles=[
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, label='Upregulated'),
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=10, label='Downregulated'),
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=10, label='Non-significant\n(AdjPval>0.01)')
        ], loc='upper right')

    plt.tight_layout()
    plt.savefig(f'./figures/{name}.pdf')
    plt.close()
    return 0

# Volcano Plot
def plot_volcano_natio(data):
    # Extracting necessary columns
    log2_fc = data["log2(fold_change)"].to_numpy()
    adj_p = data["adjusted p-value"].to_numpy()

    # Compute -log10 of the adjusted p-value
    neg_log10_p = -np.log10(adj_p)

    # Define thresholds for coloring
    fc_threshold = 1  # Fold change threshold
    p_threshold = 0.01  # P-value threshold (adjusted p-value)

    # Colors
    colors = np.where((log2_fc > fc_threshold) & (adj_p < p_threshold), "red", 
                      np.where((log2_fc < -fc_threshold) & (adj_p < p_threshold), "blue", "white"))

    # Create the plot
    sns.set(style="dark", rc={"axes.facecolor": "black", "figure.facecolor": "black"})
    plt.figure(figsize=(10, 8))
    
    # Scatter plot
    scatter = plt.scatter(log2_fc, neg_log10_p, c=colors, alpha=0.7, linewidths=0.05, edgecolor='w', s=40)

    # Add thresholds for fold change and significance
    plt.axvline(x=fc_threshold, color="red", linestyle="--", linewidth=1)
    plt.axvline(x=-fc_threshold, color="blue", linestyle="--", linewidth=1)
    # Add thresholds for significance with color transitions
    x_vals = np.linspace(min(log2_fc) - 1, max(log2_fc) + 1, 1000)
    y_threshold = -np.log10(p_threshold)
    plt.plot(x_vals[x_vals < -fc_threshold], [y_threshold] * sum(x_vals < -fc_threshold), color="blue", linestyle="--", linewidth=1)
    plt.plot(x_vals[(x_vals >= -fc_threshold) & (x_vals <= fc_threshold)], 
             [y_threshold] * sum((x_vals >= -fc_threshold) & (x_vals <= fc_threshold)), 
             color="white", linestyle="--", linewidth=1)
    plt.plot(x_vals[x_vals > fc_threshold], [y_threshold] * sum(x_vals > fc_threshold), color="red", linestyle="--", linewidth=1)

    # Annotate top points
    #top_features = data.sort("adjusted p-value").head(10)
    #for _, row in top_features.to_pandas().iterrows():
    #    plt.text(row[1], -np.log10(row[3]), row[0], fontsize=8, ha='right')

    # Labels and Title
    plt.title(f"{wanted_tp_comp[1]} vs {wanted_tp_comp[0]}\n{data.shape[0]} Peaks tested", fontsize=16,color='w')
    plt.xlabel("Log2(Fold Change)", fontsize=14, color='w')
    plt.ylabel("-Log10(Adjusted P-Value)", fontsize=14, color='w')

    # Customize ticks to white
    plt.tick_params(axis='x', colors='white', labelsize=12)
    plt.tick_params(axis='y', colors='white', labelsize=12)

    # Customize axis spines to match background or remain visible
    plt.gca().spines['bottom'].set_color('white')
    plt.gca().spines['left'].set_color('white')
    plt.gca().spines['top'].set_color('white')
    plt.gca().spines['right'].set_color('white')

    # Custom Legend
    plt.legend(handles=[
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, label='Upregulated'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=10, label='Downregulated'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='white', markersize=10, label='Non-significant\n(AdjPval>0.01)')
    ], loc='upper right', fontsize=12, labelcolor="white", facecolor="black")

    plt.tight_layout()
    plt.savefig(f'./figures/VolcanoPlotNatio_DiffTest.pdf')
    plt.close()

    return 0


def test_overall_accessibility(data, significance_level=0.05):
    data = data.filter(pl.col("adjusted p-value") < 0.001)
    # Extract log2 fold changes
    log2_fc = data["log2(fold_change)"].to_numpy()    
    
    # Perform one-sample t-test
    t_stat, p_value_ttest = ttest_1samp(log2_fc, popmean=0, alternative='less')

    # Perform Wilcoxon signed-rank test
    w_stat, p_value_wilcoxon = wilcoxon(log2_fc - 0, alternative='less', mode='approx')

    # Summarize results
    result = {
        "mean_log2_fc": np.mean(log2_fc),
        "median_log2_fc": np.median(log2_fc),
        "t-test": {"t_stat": t_stat, "p_value": p_value_ttest},
        "wilcoxon": {"w_stat": w_stat, "p_value": p_value_wilcoxon},
        "interpretation_ttest": "Significant" if p_value_ttest < significance_level else "Not Significant",
        "interpretation_wilcoxon": "Significant" if p_value_wilcoxon < significance_level else "Not Significant"
    }
    return result

def test_overall_accessibility_PosNeg(data, significance_level=0.05):
    data = data.filter(pl.col("adjusted p-value") < 0.001)
    # Extract log2 fold changes
    log2_fc = data["log2(fold_change)"].to_numpy()    
    negative_fc = data.filter(pl.col("log2(fold_change)") < 0)
    positive_fc = data.filter(pl.col("log2(fold_change)") > 0)


    # Extract values for plotting
    negative_fc_values = negative_fc["log2(fold_change)"].to_numpy()
    positive_fc_values = positive_fc["log2(fold_change)"].to_numpy()
    negative_fc_values = - negative_fc_values
    
    # Perform one-sample t-test
    t_stat, p_value_ttest = ttest_ind(negative_fc_values,positive_fc_values, equal_var=False, alternative='two-sided')

    # Summarize results
    result = {
        "mean_log2_fc": np.mean(log2_fc),
        "median_log2_fc": np.median(log2_fc),
        "t-test": {"t_stat": t_stat, "p_value": p_value_ttest},
        "interpretation_ttest": "Significant" if p_value_ttest < significance_level else "Not Significant",
    }
    return result


def visualize_accessibility_testPosNeg(data, test_results):
    data = data.filter(pl.col("adjusted p-value") < 0.001)

    # Simulated data for negative_fc_values and positive_fc_values
    negative_fc = data.filter(pl.col("log2(fold_change)") < 0)
    positive_fc = data.filter(pl.col("log2(fold_change)") > 0)


    # Extract values for plotting
    negative_fc_values = negative_fc["log2(fold_change)"].to_numpy()
    positive_fc_values = positive_fc["log2(fold_change)"].to_numpy()
    negative_fc_values = - negative_fc_values
    
    # Calculate statistics for the two groups
    mean_negative = np.mean(negative_fc_values)
    median_negative = np.median(negative_fc_values)
    mean_positive = np.mean(positive_fc_values)
    median_positive = np.median(positive_fc_values)

    # Create the figure and axes
    fig, ax = plt.subplots(figsize=(10, 6))

    # Boxplot data and customization
    box_data = [negative_fc_values, positive_fc_values]
    box_colors = ['#1f77b4', '#ff7f0e']
    box = ax.boxplot(box_data, patch_artist=True, widths=0.6, vert=True, labels=["Negative FC", "Positive FC"])

    # Style the box plots
    for patch, color in zip(box['boxes'], box_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    # Add means and medians as points and lines
    ax.scatter(1, mean_negative, color='black', zorder=3, label='Mean', marker='o', s=100)
    ax.scatter(2, mean_positive, color='black', zorder=3, marker='o', s=100)
    ax.hlines([median_negative, median_positive], [0.75, 1.75], [1.25, 2.25], colors='darkred', label='Median', linestyles='--', linewidth=2)

    # Add test results as text
    test_text = (
        f"T-Test Results:\n"
        f"t-statistic: {test_results['t-test']['t_stat']}\n"
        f"p-value: {test_results['t-test']['p_value']}"
    )
    ax.text(1.5, -1.5, test_text, fontsize=12, bbox=dict(facecolor='white', alpha=0.8, boxstyle='round'))

    # Beautify the plot
    ax.set_title("Comparison of Negative and Positive Chromatin Accessibility FC Values", fontsize=16, weight='bold')
    ax.set_ylabel("FC Values", fontsize=14)
    ax.tick_params(axis='both', labelsize=12)
    ax.grid(axis='y', linestyle='--', alpha=0.7)
    ax.legend(loc='upper left', fontsize=12, frameon=False)

    # Show the plot
    plt.tight_layout()
    plt.savefig(f'./figures/BoxPlot_DiffTestStatsPosNeg.pdf')
    plt.close()

    return 0

def visualize_accessibility_distrib(data):
    # Split the DataFrame into negative and positive log fold changes
    negative_fc = data.filter(pl.col("log2(fold_change)") < 0)
    positive_fc = data.filter(pl.col("log2(fold_change)") > 0)


    # Extract values for plotting
    negative_fc_values = negative_fc["log2(fold_change)"].to_numpy()
    positive_fc_values = positive_fc["log2(fold_change)"].to_numpy()
    negative_fc_values = - negative_fc_values

    # Create the plot
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(10, 8))

    # Plot density for negative and positive log fold changes
    sns.kdeplot(negative_fc_values, color="blue", fill=True, alpha=0.6, linewidth=2, label=f"Negative FC: {len(negative_fc_values)} Regions")
    sns.kdeplot(positive_fc_values, color="red", fill=True, alpha=0.6, linewidth=2, label=f"Positive FC: {len(positive_fc_values)} Regions")

    # Get the maximum y-value for the density curves to position the text
    negative_peak = sns.kdeplot(negative_fc_values, color="blue").get_lines()[0].get_data()[1].max()
    positive_peak = sns.kdeplot(positive_fc_values, color="red").get_lines()[0].get_data()[1].max()

    # Add labels and legend
    plt.title("Distribution of Log2(Fold Changes)", fontsize=16, color="k")
    plt.xlabel("Log2(Fold Change)", fontsize=14, color="k")
    plt.ylabel("Density", fontsize=14, color="k")
    plt.legend(fontsize=12, frameon=False, loc="upper right", labelcolor="k")

    # Tight layout and show the plot
    plt.tight_layout()
    plt.savefig(f'./figures/Distrib_DiffTest_FCs_Signi_SameSide.pdf')
    plt.close()
    return 0


def visualize_accessibility_distrib_S(data):

    # Extract values for plotting
    fc_values = data["log2(fold_change)"].to_numpy()

    # Create the plot
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(10, 8))

    # Plot density for negative and positive log fold changes
    sns.kdeplot(fc_values, color="grey", fill=True, alpha=0.6, linewidth=2)

    # Add labels and legend
    plt.title(f"Distribution of Log2(Fold Changes)\n{len(fc_values)} Peaks", fontsize=16, color="k")
    plt.xlabel("Log2(Fold Change)", fontsize=14, color="k")
    plt.ylabel("Density", fontsize=14, color="k")

    # Tight layout and show the plot
    plt.tight_layout()
    plt.savefig(f'./figures/Distrib_FCs_Signi_All.pdf')
    plt.close()
    return 0


def visualize_accessibility_distribPval(data):
    # Split the DataFrame into negative and positive log fold changes
    negative_fc = data.filter(pl.col("log2(fold_change)") < 0)
    positive_fc = data.filter(pl.col("log2(fold_change)") > 0)

    # Extract values for plotting
    negative_fc_values = negative_fc.filter(pl.col("adjusted p-value") != 0)["adjusted p-value"].to_numpy()
    positive_fc_values = positive_fc.filter(pl.col("adjusted p-value") != 0)["adjusted p-value"].to_numpy()

    # Create the plot
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(10, 8))

    # Plot density for negative and positive log fold changes
    sns.kdeplot(negative_fc_values, color="blue", fill=True, alpha=0.6, linewidth=2, label=f"Negative FC: {len(negative_fc_values)} Regions")
    sns.kdeplot(positive_fc_values, color="red", fill=True, alpha=0.6, linewidth=2, label=f"Positive FC: {len(positive_fc_values)} Regions")

    # Get the maximum y-value for the density curves to position the text
    negative_peak = sns.kdeplot(negative_fc_values, color="blue").get_lines()[0].get_data()[1].max()
    positive_peak = sns.kdeplot(positive_fc_values, color="red").get_lines()[0].get_data()[1].max()

    # Add labels and legend
    plt.title("Distribution of adjusted p-value", fontsize=16, color="k")
    plt.xlabel("adjusted p-value", fontsize=14, color="k")
    plt.ylabel("Density", fontsize=14, color="k")
    plt.legend(fontsize=12, frameon=False, loc="upper right", labelcolor="k")

    # Tight layout and show the plot
    plt.tight_layout()
    plt.savefig(f'./figures/Distrib_PVals.pdf')
    plt.close()
    return 0

def treat_motifs(moti):
    moti = moti.with_columns(
        pl.col("id")
        .str.replace(r"\.", " ")  # Replace ".", " ", and "::" with a space
        .str.replace(r"::", " ")
        .str.split(" ")  # Split by space
        .list.first()
        .str.to_uppercase()  # Get the first element from the list
        .alias("id_group")  # Rename the new column
    )

    moti = moti.with_columns(
        pl.col("log2(fold change)").abs().alias("abs_log2_fold_change")
    )

    # Sort the DataFrame first
    moti = moti.sort("abs_log2_fold_change", descending=True)

    # Use unique to get the first occurrence of each group
    moti = moti.unique(subset="id_group", keep="first")

    moti.drop("id")
    moti = moti.with_columns(
    pl.col("id_group").alias("id")  # Replace 'id' with values from 'id_group'
    )

    # Drop unnecessary columns
    moti = moti.drop(["id_group","name", "family", "abs_log2_fold_change"], strict=False)

    return moti

def scatter_accessibility(polars_df_with_counts):
    polars_df_with_counts = polars_df_with_counts.filter(pl.col("adjusted p-value") < 0.001)
    # Extract relevant columns
    day28 = np.log(polars_df_with_counts["sum_counts_Day28"].to_numpy()+1)
    intact = np.log(polars_df_with_counts["sum_counts_Intact"].to_numpy()+1)
    fc = polars_df_with_counts["log2(fold_change)"].to_numpy()

    # Create a colormap for fold change (blue -> white -> red)
    cmap = sns.diverging_palette(240, 10, as_cmap=True)

    # Scatter plot
    plt.figure(figsize=(10, 7))
    sc = plt.scatter(day28, intact, c=fc, cmap=cmap, edgecolor='k', linewidth=0.1 )

    # Add colorbar for fold change
    cb = plt.colorbar(sc)
    cb.set_label('log2(Fold Change)', fontsize=12)
    
    # Set axis limits (same for x and y)
    #min_val = min(day28.min(), intact.min())
    #max_val = max(day28.max(), intact.max())
    #plt.xlim(min_val#, max_val)
    #plt.ylim(min_val, max_val)

    # Add a diagonal dotted line (y = x)
    #plt.plot([min_val, max_val], [min_val, max_val], linestyle='--', color='gray', linewidth=1)

    # Customize axes
    #plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel('Counts per region in CD28 (log scale)', fontsize=14)
    plt.ylabel('Counts per region in Intact (log scale)', fontsize=14)

    # Add grid and improve aesthetics
    plt.grid(visible=True, which='both', linestyle='--', linewidth=0.4, alpha=0.8)
    plt.title('Regions Accessibility Comparison: CD28 vs Intact', fontsize=16)
    sns.despine()

    # Save or display the plot
    plt.tight_layout()
    plt.savefig("scatterplot_accessibility.png", dpi=300)  # Save as high-resolution file
    plt.close()

    return 0

# Function to create a plot
def plot_delta_rank(data, delta_col, title,num_ex=20):
    x_label="Index"
    y_label=delta_col

    data = data.sort_values(by=delta_col, ascending=False).reset_index()
    extreme_points = pd.concat([data.iloc[:num_ex],data.iloc[-num_ex:]])

    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(8, 12))  # Create figure and axis
    
    # Scatterplot
    scatter = sns.scatterplot(
        ax=ax,
        x=data.index,
        y=data[delta_col],
        hue=data[delta_col],
        palette="coolwarm",
        s=70,
        edgecolor="k",
        linewidth=0.2,
        legend=None
    )

    texts = []  # Store text objects for adjustment
    if 'JASPAR' in title:
        # Annotate the 10 most extreme points
        for _, row in extreme_points.iterrows():
            texts.append(
                ax.text(
                    row.name, row[delta_col], row["id"].split(' ')[1],
                    ha='right' if row[delta_col] > 0 else 'left' , va='bottom' if row[delta_col] > 0 else 'top',
                    fontsize=10, color='black',weight="semibold"
                )
            )
    elif 'H12' in title:
        # Annotate the 10 most extreme points
        for _, row in extreme_points.iterrows():
            texts.append(
                ax.text(
                    row.name, row[delta_col], row["id"].split('.')[0],
                    ha='center' if row[delta_col] > 0 else 'left' , va='bottom' if row[delta_col] > 0 else 'top',
                    fontsize=10, color='black',weight="semibold"
                )
            )
    
    # Title and labels
    ax.set_title(title, fontsize=16, weight="bold")
    ax.set_xlabel(x_label, fontsize=14)
    ax.set_ylabel(y_label, fontsize=14)
    ax.set_xlim(-200,data.shape[0]+200)

    # Adjust text to avoid overlap
    adjust_text(
        texts, ax=ax, avoid_self=True, arrowprops=dict(arrowstyle="-", color='gray', lw=0.5),
        expand_points=(5, 5),  # Maximum expansion from data points
        expand_text=(4, 5),    # Maximum spacing among text labels
        force_text=5,          # Extremely strong repelling force between labels
        force_points=5,        # Extremely strong repelling force from points
        force_objects=3        # Strong force to avoid overlapping other plot elements
    )

    # Customize the color bar
    #norm = plt.Normalize(data[delta_col].min(), data[delta_col].max())
    #sm = plt.cm.ScalarMappable(cmap="coolwarm", norm=norm)
    #sm.set_array([])
    #cbar = fig.colorbar(sm, ax=ax, label="Delta Rank", orientation="vertical", pad=0.02)
    
    # Add gridlines for clarity
    ax.grid(color="gray", linestyle="--", linewidth=0.5, alpha=0.7)
    plt.tight_layout()
    plt.savefig(f'./figures/{title}.pdf')
    plt.close()

    return 0

# Function to create a scatter plot for delta_log2_fc vs delta_rank_adj_pval
def plot_delta_comparison(data, title,col_adjPval,num_ex=20, ranklogfc=True):
    
    if ranklogfc:
        x_label="Delta Rank log2(Fold Change)"
        deltacol_fc = "delta_rank_log2_fc"
        add_title = "rankFC"
    else:
        x_label="Delta log2(Fold Change)"
        deltacol_fc = "delta_log2_fc"
        add_title = "FC"

    #y_label="Delta Rank Adjusted P-Value"
    
    # Sort by absolute delta rank for each metric
    df_log2_fc = data.sort_values(by=deltacol_fc, ascending=False).reset_index()
    #df_adj_pval = data.sort_values(by="delta_rank_adj_pval", ascending=False).reset_index()
    df_adj_pval = data.sort_values(by=col_adjPval, ascending=False).reset_index()

    # Identify the top 10 extreme points
    extreme_log2_fc = pd.concat([df_log2_fc.iloc[:num_ex],df_log2_fc.iloc[-num_ex:]])
    extreme_adj_pval = pd.concat([df_adj_pval.iloc[:num_ex]])

    extreme_points=pd.concat([extreme_log2_fc, extreme_adj_pval])
    
    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(12, 8))  # Create figure and axis
    
    # Scatterplot
    scatter = sns.scatterplot(
        ax=ax,
        x=data[deltacol_fc],
        y=data[col_adjPval],
        hue=data[deltacol_fc],
        palette="coolwarm",
        s=70,
        edgecolor="k",
        linewidth=0.2,
        legend=None
    )

    # Annotate extreme points
    texts = []
    for _, row in extreme_points.iterrows():
        if 'JASPAR' in title:
            label = row["id"].split(' ')[1]  # Use the second part of the ID for JASPAR
        elif 'H12' in title:
            label = row["id"].split('.')[0]  # Use the first part of the ID for HOCOMOCO
        else:
            label = row["id"]
        
        texts.append(
            ax.text(
                row[deltacol_fc], row[col_adjPval], label,
                ha='center', va='bottom', fontsize=10, color='black', weight="semibold"
            )
        )
    
    # Title and labels
    ax.set_title(title, fontsize=16, weight="semibold")
    ax.set_xlabel(x_label, fontsize=14)
    ax.set_ylabel(col_adjPval, fontsize=14)
    
    # Add gridlines
    ax.grid(color="gray", linestyle="--", linewidth=0.5, alpha=0.7)

    # Adjust text to avoid overlap
    adjust_text(
        texts, ax=ax, avoid_self=True, 
        arrowprops=dict(arrowstyle="-", color='gray', lw=0.5),
        expand_points=(5, 5), expand_text=(4, 5),
        force_text=6, force_points=6, force_objects=4
    )

    plt.tight_layout()
    plt.savefig(f'./figures/{title}_{add_title}.pdf')
    plt.close()

    return 0
###########
# CODE
###########

# Load ATAC data (and RNA for annot)
adata_ATAC = snap.read(atac_h5ad, backed=None)
adata_AllAnnot = sc.read_h5ad(rna_h5ad)

# Add RNA All Annot to adata_ATAC
adata_ATAC = adata_ATAC[np.isin(adata_ATAC.obs.index,adata_AllAnnot.obs.index)]
adata_ATAC.obs = adata_ATAC.obs.reindex(index=sorted(adata_ATAC.obs.index))
adata_AllAnnot.obs = adata_AllAnnot.obs.reindex(index=sorted(adata_AllAnnot.obs.index))

assert (adata_ATAC.obs.index == adata_AllAnnot.obs.index).all()
adata_ATAC.obs['Annotation'] = adata_AllAnnot.obs['Annotation']
if np.sum(adata_ATAC.obs['Annotation'] == 'Other')==0:
    adata_ATAC.obs['Annotation'] = adata_ATAC.obs['Annotation'].cat.remove_categories('Other')

del adata_AllAnnot

desired_order = ['Intact', 'Day28', 'RegenDay1', 'RegenDay2', 'RegenDay3']
adata_ATAC.obs['timePoint'] = adata_ATAC.obs['timePoint'].astype(pd.CategoricalDtype(categories=desired_order, ordered=True))

# Plotting
sc.pl.umap(adata_ATAC,size=400000/adata_ATAC.n_obs,alpha=0.8,color=['Annotation'],save=Experiment+'_TF_Annot',show=False)
sc.pl.umap(adata_ATAC,size=400000/adata_ATAC.n_obs,alpha=0.8,color=['timePoint'],save=Experiment+'_TF_timePoints',show=False)
sc.pl.umap(adata_ATAC,size=400000/adata_ATAC.n_obs,alpha=0.8,color=['batch'],save=Experiment+'_TF_batch',show=False)
sc.pl.umap(adata_ATAC,size=400000/adata_ATAC.n_obs,alpha=0.8,color=['batch'],save=Experiment+'_TF_batch',show=False)
sc.pl.umap(adata_ATAC,size=400000/adata_ATAC.n_obs,alpha=0.8,color=["n_fragment", "frac_dup", "frac_mito"],save=Experiment+'_Counts_QC',show=False)
sc.pl.umap(adata_ATAC,size=400000/adata_ATAC.n_obs,alpha=0.8,color=["FRiP", "tsse"],save=Experiment+'_PctCounts_QC',show=False)
leiR = 1
leiRes = f'leiden_res{leiR}'
sc.pl.umap(adata_ATAC,size=400000/adata_ATAC.n_obs,alpha=0.8,color=leiRes,legend_loc="on data",save=f"{Experiment}_cluster{leiRes}",title=Experiment,show=False)


# Performing differential testing of accessibility

wanted_tp_comp = ('Intact','Day28') #('Intact','Day28') #  ('Day28','RegenDay1') ('Day28','RegenDay2')
luminal_celltype = ['L2', 'Ambig_L1|L2'] #['L1'] # ['L2', 'Ambig_L1|L2']# ],['Basal']]

#CHECK luminal_celltype AND wanted_tp_comp
print('Cell type(s):',luminal_celltype)
print('Time Points Compared:'," / ".join(list(wanted_tp_comp)))
adata_ATAC_CellTypeSubset = adata_ATAC[np.isin(adata_ATAC.obs['Annotation'],luminal_celltype)].copy()
#Kick citrate out
adata_ATAC_CellTypeSubset = adata_ATAC_CellTypeSubset[~np.isin(adata_ATAC_CellTypeSubset.obs['batch'],['WK-1585_INTACT_AP_BL6_Citrate'])].copy()

qcatac = adata_ATAC[np.isin(adata_ATAC.obs['batch'],['WK-1501_BL6_INTACT_AP_Test3','WK-1585_INTACT_AP_BL6_Contrl','WK-1585_Castrate_Day28_AP_BL6'])].copy()
qcatac.obs['log1p_n_fragment'] = np.log1p(adata_ATAC_CellTypeSubset.obs['n_fragment'])
sc.pl.umap(qcatac,size=400000/adata_ATAC.n_obs,alpha=0.8,color=["n_fragment","log1p_n_fragment","FRiP", "tsse"],save=Experiment+'_QC_Intact',show=False)
sc.pl.umap(qcatac,size=400000/adata_ATAC.n_obs,alpha=0.8,color=["Annotation","batch"],save=Experiment+'_QC_Intact_post',show=False)
sc.pl.violin(qcatac, ["n_fragment","log1p_n_fragment","FRiP", "tsse"],multi_panel=True,groupby='batch',rotation=90,save=Experiment+'_QC_Intact.pdf')
sc.pl.violin(qcatac, ["n_fragment","log1p_n_fragment","FRiP", "tsse"],multi_panel=True,rotation=90,save=Experiment+'_QC_Intact_Samplestogether.pdf')

#group2 is Intact
cell_group2 = adata_ATAC_CellTypeSubset[np.isin(adata_ATAC_CellTypeSubset.obs['timePoint'],wanted_tp_comp[0])].obs.index
#group1 is CD28
cell_group1 = adata_ATAC_CellTypeSubset[np.isin(adata_ATAC_CellTypeSubset.obs['timePoint'],wanted_tp_comp[1])].obs.index

#positive FCs: features enriched in group1 (CD28). 
#negative FCs: features enriched in group2 (Intact)
diff_test_DF = snap.tl.diff_test(adata_ATAC_CellTypeSubset, cell_group1, cell_group2, covariates=None, direction='both', min_log_fc=0, min_pct=0.05)

# Plotting results
plot_volcano(diff_test_DF)
plot_volcano_natio(diff_test_DF)
diff_test_DF.write_csv(f'diff_test_DF_{luminal_celltype[0]}.tsv',separator='\t')

# Select Significant regions and check if there's stat upregulation
diff_test_DF_signi = diff_test_DF.filter(pl.col("adjusted p-value") < 0.001)
diff_test_DF_signi = diff_test_DF_signi.filter(pl.col("log2(fold_change)").abs() > 0.5)

visualize_accessibility_distribPval(diff_test_DF_signi) # All signi PVal values
visualize_accessibility_distrib_S(diff_test_DF_signi) # All signi FC values
visualize_accessibility_distrib(diff_test_DF_signi) # All signi FC values BUT negative and positive on the same side

#test and plot boxplots on difference of FC values between neg and pos
visualize_accessibility_testPosNeg(diff_test_DF_signi, test_overall_accessibility_PosNeg(diff_test_DF_signi))

# Separate into pos and neg regulation regions and motif analysis
regions_signi_negative_fc = diff_test_DF_signi.filter(pl.col("log2(fold_change)") < 0)["feature name"].to_list()
regions_signi_positive_fc = diff_test_DF_signi.filter(pl.col("log2(fold_change)") > 0)["feature name"].to_list()
regions_background = diff_test_DF["feature name"].to_list()

regions = {
    "CD28-Enriched": regions_signi_positive_fc,
    "CD28-Depleted": regions_signi_negative_fc,
}

# motif analysis with H12
TFmotifsHC = '/mnt/etemp/ahrmad/wouter/motif_databases/20241214_Redundant_Motifs/H12CORE.meme'
motifsHC = snap.tl.motif_enrichment(
    motifs=snap.read_motifs(TFmotifsHC),
    regions=regions,
    background=regions_background,
    genome_fasta=snap.genome.mm10,
    method="hypergeometric"
)

# motif analysis with JASPAR
TFmotifsJA = '/mnt/etemp/ahrmad/wouter/motif_databases/20241214_Redundant_Motifs/JASPAR2024_CORE_vertebrates_redundant_pfms.meme'
motifsJA = snap.tl.motif_enrichment(
    motifs=snap.read_motifs(TFmotifsJA),
    regions=regions,
    background=regions_background,
    genome_fasta=snap.genome.mm10,
    method="hypergeometric"
)


#Plot and save results
for nameDB,motifs in [("JASPAR2024",motifsJA),("HOCOMOCOv12",motifsHC)]:
    for name_key,PL_value in motifs.items():
        PL_value.write_csv(f'Motifs_{nameDB}_{name_key}_{luminal_celltype[0]}.tsv',separator='\t')

    max_adjPval = 0.01
    for log_fc_min in [0.5,1,2]:
        print(f'PROCESSING {nameDB} - {log_fc_min}')
        num_TFs = output_TFs_tsv(motifs, log_fc_min, max_adjPval,f'{Experiment}_TF_motifs{nameDB}_{"_".join(luminal_celltype)}_TP_{"_".join(wanted_tp_comp)}_minFC{str(log_fc_min).replace(".", "")}.tsv')
        if num_TFs == 0:
            break
        snap.pl.motif_enrichment(motifs, max_fdr=max_adjPval, min_log_fc=log_fc_min, height=150 + num_TFs*25,width=1000, out_file=f'figures/heatmap{Experiment}_TF_motifs{nameDB}_{"_".join(luminal_celltype)}_TP_{"_".join(wanted_tp_comp)}_minFC{str(log_fc_min).replace(".", "")}.png', scale=1,show=False)



polars_df = pl.read_csv(f"diff_test_DF_{luminal_celltype[0]}.tsv",separator='\t')  # Replace with your actual file or DataFrame
for wtp in wanted_tp_comp:
    adata_inUse = adata_ATAC_CellTypeSubset[np.isin(adata_ATAC_CellTypeSubset.obs['timePoint'],wtp)].copy()
    # Load the Polars DataFrame and AnnData object
    # Ensure variable names in the Polars DataFrame match those in the AnnData object
    variable_list = polars_df["feature name"].to_list()  # Replace with your actual column name in Polars DataFrame
    variables_to_use = [var for var in variable_list if var in adata_inUse.var_names]
    # Sum counts for each variable across all cells in the AnnData object
    var_indices = np.where(np.isin(adata_inUse.var_names.values, np.array(variables_to_use)))[0]
    sums = adata_inUse[:, var_indices].X.sum(axis=0).A1
    # Convert the counts_sum dictionary to a Polars DataFrame
    counts_df = pl.DataFrame({"feature name": variables_to_use, f"sum_counts_{wtp}": sums})
    # Join the two DataFrames on the variable column
    polars_df = polars_df.join(counts_df, on="feature name", how="left",validate='1:1')
    # Save or use the resulting DataFrame
polars_df.write_csv(f"diff_test_DF_{luminal_celltype[0]}_withCounts.tsv",separator='\t')  # Replace with desired output file path


polars_df = pl.read_csv(f"diff_test_DF_{luminal_celltype[0]}_withCounts.tsv",separator='\t')  # Replace with your actual file or DataFrame
scatter_accessibility(polars_df)


#os.chdir('../figL2')
#luminal_celltype[0] = 'L2'
#polars_df = pl.read_csv(f"diff_test_DF_{luminal_celltype[0]}_withCounts.tsv",separator='\t')  # Replace with your actual file or DataFrame
#scatter_accessibility(polars_df)



##############
# Marker Regions analysis
##############

adata_ATAC_subset_comp = adata_ATAC_CellTypeSubset[np.isin(adata_ATAC_CellTypeSubset.obs['timePoint'],wanted_tp_comp)].copy()

for ba,nm in dict(Counter(adata_ATAC_subset_comp.obs['batch'])).items():
    print(f'{ba}: {nm} cells')

marker_peaks = snap.tl.marker_regions(adata_ATAC_subset_comp, groupby='timePoint', pvalue=0.16)
print(marker_peaks)

# Resulting heatmap
snap.pl.regions(adata_ATAC_subset_comp, groupby='timePoint', width=1200, height=800, peaks=marker_peaks, out_file=f'figures/heatmap{Experiment}_TF_{luminal_celltype[0]}_regionsTP_{"_".join(wanted_tp_comp)}.png',show=False)

regions_background_MR = adata_ATAC_subset_comp.var.index.to_list()

# motif analysis with H12
TFmotifsHC = '/mnt/etemp/ahrmad/wouter/motif_databases/20241214_Redundant_Motifs/H12CORE.meme'
motifsHC_MR = snap.tl.motif_enrichment(
    motifs=snap.read_motifs(TFmotifsHC),
    regions=marker_peaks,
    background=regions_background_MR,
    genome_fasta=snap.genome.mm10,
    method="hypergeometric"
)

# motif analysis with JASPAR
TFmotifsJA = '/mnt/etemp/ahrmad/wouter/motif_databases/20241214_Redundant_Motifs/JASPAR2024_CORE_vertebrates_redundant_pfms.meme'
motifsJA_MR = snap.tl.motif_enrichment(
    motifs=snap.read_motifs(TFmotifsJA),
    regions=marker_peaks,
    background=regions_background_MR,
    genome_fasta=snap.genome.mm10,
    method="hypergeometric"
)

#Plot and save results
for nameDB,motifs in [("JASPAR2024",motifsJA_MR),("HOCOMOCOv12",motifsHC_MR)]:
    for name_key,PL_value in motifs.items():
        PL_value.write_csv(f'Motifs_{nameDB}_{name_key}_{luminal_celltype[0]}_MarkerRegions.tsv',separator='\t')

    max_adjPval = 0.01
    for log_fc_min in [0.3,0.5,1]:
        print(f'PROCESSING {nameDB} - {log_fc_min}')
        num_TFs = output_TFs_tsv(motifs, log_fc_min, max_adjPval,f'{Experiment}_TF_motifs{nameDB}_{"_".join(luminal_celltype)}_TP_{"_".join(wanted_tp_comp)}_minFC{str(log_fc_min).replace(".", "")}_MarkerRegions.tsv')
        if num_TFs == 0:
            break
        snap.pl.motif_enrichment(motifs, max_fdr=max_adjPval, min_log_fc=log_fc_min, height=150 + num_TFs*25,width=1000, out_file=f'figures/heatmap{Experiment}_TF_motifs{nameDB}_{"_".join(luminal_celltype)}_TP_{"_".join(wanted_tp_comp)}_minFC{str(log_fc_min).replace(".", "")}_MarkerRegions.png', scale=1,show=False)




import polars as pl
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text

Motif_res_list = [("Motifs_HOCOMOCOv12_CD28-Enriched_L1.tsv","Motifs_HOCOMOCOv12_CD28-Depleted_L1.tsv"),("Motifs_JASPAR2024_CD28-Enriched_L1.tsv","Motifs_JASPAR2024_CD28-Depleted_L1.tsv"),("Motifs_JASPAR2024_Day28_L1_MarkerRegions.tsv","Motifs_JASPAR2024_Intact_L1_MarkerRegions.tsv"),("Motifs_HOCOMOCOv12_Day28_L1_MarkerRegions.tsv","Motifs_HOCOMOCOv12_Intact_L1_MarkerRegions.tsv")]
Motif_res_list = [("Motifs_HOCOMOCOv12_CD28-Enriched_L2.tsv","Motifs_HOCOMOCOv12_CD28-Depleted_L2.tsv"),("Motifs_JASPAR2024_CD28-Enriched_L2.tsv","Motifs_JASPAR2024_CD28-Depleted_L2.tsv"),("Motifs_JASPAR2024_Day28_L2_MarkerRegions.tsv","Motifs_JASPAR2024_Intact_L2_MarkerRegions.tsv"),("Motifs_HOCOMOCOv12_Day28_L2_MarkerRegions.tsv","Motifs_HOCOMOCOv12_Intact_L2_MarkerRegions.tsv")]
rank_method = 'dense'

for pos_df_path,neg_df_path in Motif_res_list:
    pos_df = pl.read_csv(pos_df_path,separator='\t')
    neg_df = pl.read_csv(neg_df_path,separator='\t')

    if "HOCOMOCOv12" in pos_df_path:
        tt = 'H12'
    elif "JASPAR2024" in pos_df_path:
        tt = "JASPAR2024"
    if "_MarkerRegions" in pos_df_path:
        ttt = "_MarkerRegions"
        rank_method='random'
    else:
        ttt=''

    # Adding rank columns
    pos_df = pos_df.with_columns([
        pl.col("log2(fold change)").rank(method=rank_method,descending=False).alias("rank_log2_fold_change"),
        pl.col("adjusted p-value").rank(method=rank_method,descending=True).alias("rank_adjusted_p_value")
    ])

    neg_df = neg_df.with_columns([
        pl.col("log2(fold change)").rank(method=rank_method,descending=False).alias("rank_log2_fold_change"),
        pl.col("adjusted p-value").rank(method=rank_method,descending=True).alias("rank_adjusted_p_value")
    ])

    df = pd.DataFrame({"id":pos_df["id"],"delta_log2_fc":pos_df["log2(fold change)"].cast(pl.Float64) - neg_df["log2(fold change)"].cast(pl.Float64),
                            "delta_rank_log2_fc":pos_df["rank_log2_fold_change"].cast(pl.Int32) - neg_df["rank_log2_fold_change"].cast(pl.Int32),
                            "delta_rank_adj_pval":pos_df["rank_adjusted_p_value"].cast(pl.Int32) - neg_df["rank_adjusted_p_value"].cast(pl.Int32),
                            "Positive FC adj. p-values":-pos_df["adjusted p-value"].log10().clip(-30,0),
                            "Negative FC adj. p-values":-neg_df["adjusted p-value"].log10().clip(-30,0)
                            })
    df.to_csv(f'Motifs_Delta_{tt}{ttt}.tsv',sep='\t')

    top_points=15
    # Plot delta_log2_fc
    plot_delta_rank(df, "delta_log2_fc",title=f"Delta_log2FC_{tt}{ttt}",num_ex=top_points)
    # Plot delta_rank_log2_fc
    plot_delta_rank(df, "delta_rank_log2_fc",title=f"Delta_Rank_log2FC_{tt}{ttt}",num_ex=top_points)
    # Plot delta_rank_adj_pval
    #plot_delta_rank(df, "delta_rank_adj_pval",title=f"Delta_Rank_AdjPValue_{tt}{ttt}",num_ex=top_points)
    
    plot_delta_rank(df, "Negative FC adj. p-values",title=f"NegativeFC_adj_pvalues_{tt}{ttt}",num_ex=top_points)
    plot_delta_rank(df, "Positive FC adj. p-values",title=f"PositiveFC_adj_pvalues_{tt}{ttt}",num_ex=top_points)

    # Plot delta_log2_fc vs delta_rank_adj_pval
    #plot_delta_comparison(df,title=f"Delta_Comparison_{tt}{ttt}",col_adjPval="delta_rank_adj_pval",ranklogfc=False,num_ex=top_points)
    
    plot_delta_comparison(df,title=f"Delta_Comparison_{tt}{ttt}_PosadjPval",col_adjPval="Positive FC adj. p-values",ranklogfc=False,num_ex=top_points)
    plot_delta_comparison(df,title=f"Delta_Comparison_{tt}{ttt}_NegadjPval",col_adjPval="Negative FC adj. p-values",ranklogfc=False,num_ex=top_points)
    # Plot delta_rank_log2_fc vs delta_rank_adj_pval
    #plot_delta_comparison(df,title=f"Delta_Comparison_{tt}{ttt}",col_adjPval="delta_rank_adj_pval",ranklogfc=True,num_ex=top_points)
    plot_delta_comparison(df,title=f"Delta_Comparison_{tt}{ttt}_PosadjPval",col_adjPval="Positive FC adj. p-values",ranklogfc=True,num_ex=top_points)
    plot_delta_comparison(df,title=f"Delta_Comparison_{tt}{ttt}_NegadjPval",col_adjPval="Negative FC adj. p-values",ranklogfc=True,num_ex=top_points)



import numpy as np
Motif_res_list = [("Motifs_HOCOMOCOv12_CD28-Enriched_L1.tsv","Motifs_HOCOMOCOv12_CD28-Depleted_L1.tsv"),("Motifs_JASPAR2024_CD28-Enriched_L1.tsv","Motifs_JASPAR2024_CD28-Depleted_L1.tsv"),("Motifs_JASPAR2024_Day28_L1_MarkerRegions.tsv","Motifs_JASPAR2024_Intact_L1_MarkerRegions.tsv"),("Motifs_HOCOMOCOv12_Day28_L1_MarkerRegions.tsv","Motifs_HOCOMOCOv12_Intact_L1_MarkerRegions.tsv")]
Motif_res_list = [("Motifs_HOCOMOCOv12_CD28-Enriched_L2.tsv","Motifs_HOCOMOCOv12_CD28-Depleted_L2.tsv"),("Motifs_JASPAR2024_CD28-Enriched_L2.tsv","Motifs_JASPAR2024_CD28-Depleted_L2.tsv"),("Motifs_JASPAR2024_Day28_L2_MarkerRegions.tsv","Motifs_JASPAR2024_Intact_L2_MarkerRegions.tsv"),("Motifs_HOCOMOCOv12_Day28_L2_MarkerRegions.tsv","Motifs_HOCOMOCOv12_Intact_L2_MarkerRegions.tsv")]

for pos_df_path,neg_df_path in Motif_res_list:
    pos_df = pl.read_csv(pos_df_path,separator='\t')
    pos_df=pos_df.rename({"log2(fold change)": "log2(fold_change)"})
    neg_df = pl.read_csv(neg_df_path,separator='\t')
    neg_df=neg_df.rename({"log2(fold change)": "log2(fold_change)"})

    if "HOCOMOCOv12" in pos_df_path:
        tt = 'H12'
    elif "JASPAR2024" in pos_df_path:
        tt = "JASPAR2024"
    if "_MarkerRegions" in pos_df_path:
        ttt = "_MarkerRegions"
        rank_method='random'
    else:
        ttt=''

    plot_volcano(pos_df,f"Volcano_MotifEnrichment_Pos_{tt}{ttt}",fc_threshold = 1,p_threshold = 0.01,whats_tested='Motifs',pval_lim=1e-30,lege=False)
    plot_volcano(neg_df,f"Volcano_MotifEnrichment_Neg_{tt}{ttt}",fc_threshold = 1,p_threshold = 0.01,whats_tested='Motifs',pval_lim=1e-30,lege=False)











































































# TESTING without opposite group in background
            regions_background_for_negFC = list(set(regions_background) - set(regions_signi_positive_fc))
            regions_background_for_posFC = list(set(regions_background) - set(regions_signi_negative_fc))

            regions_pos = {
                "I_CD28_positive_fc": regions_signi_positive_fc,
            }

            regions_neg = {
                "I_CD28_negative_fc": regions_signi_negative_fc,
            }

            # motif analysis with H12
            TFmotifsHC = '/mnt/etemp/ahrmad/wouter/motif_databases/H12CORE_meme_format.meme'
            motifsHC_pos = snap.tl.motif_enrichment(
                motifs=snap.read_motifs(TFmotifsHC),
                regions=regions_pos,
                background=regions_background_for_posFC,
                genome_fasta=snap.genome.mm10,
                method="hypergeometric"
            )

            motifsHC_neg = snap.tl.motif_enrichment(
                motifs=snap.read_motifs(TFmotifsHC),
                regions=regions_neg,
                background=regions_background_for_negFC,
                genome_fasta=snap.genome.mm10,
                method="hypergeometric"
            )
            motifsHC = {}
            motifsHC.update(motifsHC_pos)
            motifsHC.update(motifsHC_neg)

            # motif analysis with JASPAR
            TFmotifsJA = '/mnt/etemp/ahrmad/wouter/motif_databases/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme_NameFirst.meme'
            motifsJA_pos = snap.tl.motif_enrichment(
                motifs=snap.read_motifs(TFmotifsJA),
                regions=regions_pos,
                background=regions_background_for_posFC,
                genome_fasta=snap.genome.mm10,
                method="hypergeometric"
            )

            motifsJA_neg = snap.tl.motif_enrichment(
                motifs=snap.read_motifs(TFmotifsJA),
                regions=regions_neg,
                background=regions_background_for_negFC,
                genome_fasta=snap.genome.mm10,
                method="hypergeometric"
            )
            motifsJA = {}
            motifsJA.update(motifsJA_pos)
            motifsJA.update(motifsJA_neg)













# try and build a consensus Venn

#L1
motifs_old = pl.read_csv("/mnt/etemp/ahrmad/wouter/Wouter21_SING/annot/Wouter21_SING_TF_motifsHOCOMOCOv12_L1_TP_Intact_Day28_minFC0.3.csv")
motifs_old = motifs_old.filter(pl.col("adjusted p-value_Day28") < 0.001)
motifs_old = motifs_old.filter(pl.col("log2(fold change)_Day28") > 0.5)
motifs_old = motifs_old.filter(~motifs_old["id"].str.contains(r"H12CORE\.[^0]"))

motifs_new = pl.read_csv("/mnt/etemp/ahrmad/wouter/Wouter21_SING/annot/figL1/figures/Wouter21_SING_TF_motifsHOCOMOCOv12_L1_TP_Intact_Day28_minFC1.tsv",separator='\t')
motifs_new = motifs_new.filter(pl.col("adjusted p-value_I_CD28_positive_fc") < 0.001)
motifs_new = motifs_new.filter(pl.col("log2(fold change)_I_CD28_positive_fc") > 0.5)
motifs_new = motifs_new.filter(~motifs_new["id"].str.contains(r"H12CORE\.[^0]"))

# Extract unique IDs as sets
ids_old = set(motifs_old["id"].to_list())
ids_new = set(motifs_new["id"].to_list())

# Create the Venn diagram
plt.figure(figsize=(8, 8))
venn = venn2([ids_old, ids_new], ("Motifs Old", "Motifs New"))

# Customize appearance
for text in venn.set_labels:
    if text:
        text.set_fontsize(14)
        text.set_color("navy")
for text in venn.subset_labels:
    if text:
        text.set_fontsize(12)
        text.set_color("darkgreen")

# Add title
plt.title("Venn Diagram of IDs in Motifs Old (minFC0.5) and New (minFC1)", fontsize=16, color="darkred")

# Save to file
plt.tight_layout()
plt.savefig("TF_Old05_New1L1_Venn.pdf", format="pdf")
plt.close()




# try and build a consensus TABLE/HEATMAP

#Start with HC and JA made above

motifsHC_P = motifsHC['I_CD28_positive_fc']
motifsHC_P = treat_motifs(motifsHC_P)
motifsHC_P = motifsHC_P.rename({col: f"{col}_HOCOMOCOv12" for col in motifsHC_P.columns})

motifsHC_N = motifsHC['I_CD28_negative_fc']
motifsHC_N = treat_motifs(motifsHC_N)
motifsHC_N = motifsHC_N.rename({col: f"{col}_HOCOMOCOv12" for col in motifsHC_N.columns})

motifsJA_P = motifsJA['I_CD28_positive_fc']
motifsJA_P = treat_motifs(motifsJA_P)
motifsJA_P = motifsJA_P.rename({col: f"{col}_JASPAR2024" for col in motifsJA_P.columns})

motifsJA_N = motifsJA['I_CD28_negative_fc']
motifsJA_N = treat_motifs(motifsJA_N)
motifsJA_N = motifsJA_N.rename({col: f"{col}_JASPAR2024" for col in motifsJA_N.columns})

motifs_P = motifsHC_P.join(motifsJA_P,
left_on='id_HOCOMOCOv12',
right_on='id_JASPAR2024',
how='inner',
validate='1:1',
suffix=''
)

motifs_N = motifsHC_N.join(motifsJA_N,
left_on='id_HOCOMOCOv12',
right_on='id_JASPAR2024',
how='inner',
validate='1:1',
suffix=''
)

# Filter the data based on adjusted p-values
motifs_filtered = motifs_P.filter(
    (pl.col("adjusted p-value_HOCOMOCOv12") <= 0.01) & 
    (pl.col("adjusted p-value_JASPAR2024") <= 0.01)
)

# Define thresholds for fold change
thresholds = [0.5, 1, 2]

# Loop through the thresholds to create individual heatmaps for each
for threshold in thresholds:
    # Add threshold columns for both HOCOMOCOv12 and JASPAR2024 fold changes
    motifs_filtered_threshold = motifs_filtered.with_columns([
        (pl.col("log2(fold change)_HOCOMOCOv12").abs() > threshold).alias(f"log2fc_HOCOMOCOv12_{threshold}"),
        (pl.col("log2(fold change)_JASPAR2024").abs() > threshold).alias(f"log2fc_JASPAR2024_{threshold}")
    ])
    
    # Keep only the columns you need (id and the thresholded columns)
    motifs_filtered_heatmap = motifs_filtered_threshold.select([
        "id_HOCOMOCOv12",  # Keep the "id" column
        f"log2fc_HOCOMOCOv12_{threshold}",
        f"log2fc_JASPAR2024_{threshold}"
    ])

    # Convert boolean values to 1/0 for heatmap visualization
    motifs_filtered_heatmap = motifs_filtered_heatmap.select(
        ["id_HOCOMOCOv12"] + [
            pl.col(col).cast(pl.Int8) for col in motifs_filtered_heatmap.columns if col != "id_HOCOMOCOv12"
        ]
    )
    motifs_filtered_heatmap = motifs_filtered_heatmap.rename({"id_HOCOMOCOv12": "id"})
    # Convert to pandas for seaborn heatmap
    heatmap_df = motifs_filtered_heatmap.to_pandas()
    heatmap_df = heatmap_df[~((heatmap_df[f"log2fc_HOCOMOCOv12_{threshold}"] == 0) & (heatmap_df[f"log2fc_JASPAR2024_{threshold}"] == 0))]

    # Set the "id" column as the index
    heatmap_df.set_index("id", inplace=True)

    # Create the heatmap for the current threshold
    plt.figure(figsize=(8, 8))  # Adjust size as needed
    try:
        sns.heatmap(heatmap_df, cmap="RdYlGn", annot=False, cbar=False, linewidths=0, yticklabels=True)
    except ValueError:
        continue
    # Adjust y-axis labels to be smaller, without rotation and with proper alignment
    plt.yticks(fontsize=2)  # Adjust font size for y labels

    # Add title for the heatmap
    plt.title(f"Heatmap of log2 Fold Change Threshold {threshold} (Red = False, Green = True)")

    # Save the heatmap as a PDF
    plt.savefig(f'./figures/HEATMAP_{threshold}_P.pdf')

    # Close the plot to avoid overlap
    plt.close()

# Add to the consensus the old analysis

motifs_old = pl.read_csv("./figures/Wouter21_SING_TF_motifsHOCOMOCOv12_L1_TP_Intact_Day28_minFC0.3.csv")
motifs_old_P = motifs_old.select(
    [
        pl.col("id").alias("id"),
        pl.col("log2(fold change)_Intact").alias("log2(fold change)"),
        pl.col("adjusted p-value_Intact").alias("adjusted p-value")
    ]
)

motifs_old_N = motifs_old.select(
    [
        pl.col("id").alias("id"),
        pl.col("log2(fold change)_Day28").alias("log2(fold change)"),
        pl.col("adjusted p-value_Day28").alias("adjusted p-value")
    ]
)

motifs_old_P = treat_motifs(motifs_old_P)
motifs_old_N = treat_motifs(motifs_old_N)

motifs_P = motifs_P.join(motifs_old_P,
left_on='id_HOCOMOCOv12',
right_on='id',
how='inner',
validate='1:1',
suffix='')

motifs_N = motifs_N.join(motifs_old_N,
left_on='id_HOCOMOCOv12',
right_on='id',
how='inner',
validate='1:1',
suffix='')

motifs_N = motifs_N.rename({
    "log2(fold change)": "log2(fold change)_old", 
    "adjusted p-value": "adjusted p-value_old"
})
motifs_P = motifs_P.rename({
    "log2(fold change)": "log2(fold change)_old", 
    "adjusted p-value": "adjusted p-value_old"
})

motifs_N = motifs_N.with_columns(
    (-pl.col("log2(fold change)_HOCOMOCOv12")).alias("log2(fold change)_HOCOMOCOv12")
)
motifs_N = motifs_N.with_columns(
    (-pl.col("log2(fold change)_JASPAR2024")).alias("log2(fold change)_JASPAR2024")
)
motifs_N.write_csv(f'./figures/motifs_N_L1_Intact_CD28.tsv',separator='\t')
motifs_P.write_csv(f'./figures/motifs_P_L1_Intact_CD28.tsv',separator='\t')






# Assume `motifs_N` is your original Polars DataFrame

# Filter the data based on adjusted p-values
motifs_filtered = motifs_P.filter(
    (pl.col("adjusted p-value_HOCOMOCOv12") <= 0.01) & 
    (pl.col("adjusted p-value_JASPAR2024") <= 0.01)
)

# Define thresholds for fold change
thresholds = [0.5, 1, 2]

# Loop through the thresholds to create individual heatmaps for each
for threshold in thresholds:
    # Add threshold columns for both HOCOMOCOv12 and JASPAR2024 fold changes
    motifs_filtered_threshold = motifs_filtered.with_columns([
        (pl.col("log2(fold change)_HOCOMOCOv12").abs() > threshold).alias(f"log2fc_HOCOMOCOv12_{threshold}"),
        (pl.col("log2(fold change)_JASPAR2024").abs() > threshold).alias(f"log2fc_JASPAR2024_{threshold}")
    ])
    
    # Keep only the columns you need (id and the thresholded columns)
    motifs_filtered_heatmap = motifs_filtered_threshold.select([
        "id_HOCOMOCOv12",  # Keep the "id" column
        f"log2fc_HOCOMOCOv12_{threshold}",
        f"log2fc_JASPAR2024_{threshold}"
    ])

    # Convert boolean values to 1/0 for heatmap visualization
    motifs_filtered_heatmap = motifs_filtered_heatmap.select(
        ["id_HOCOMOCOv12"] + [
            pl.col(col).cast(pl.Int8) for col in motifs_filtered_heatmap.columns if col != "id_HOCOMOCOv12"
        ]
    )
    motifs_filtered_heatmap = motifs_filtered_heatmap.rename({"id_HOCOMOCOv12": "id"})
    # Convert to pandas for seaborn heatmap
    heatmap_df = motifs_filtered_heatmap.to_pandas()
    heatmap_df = heatmap_df[~((heatmap_df[f"log2fc_HOCOMOCOv12_{threshold}"] == 0) & (heatmap_df[f"log2fc_JASPAR2024_{threshold}"] == 0))]

    # Set the "id" column as the index
    heatmap_df.set_index("id", inplace=True)

    # Create the heatmap for the current threshold
    plt.figure(figsize=(8, 8))  # Adjust size as needed
    try:
        sns.heatmap(heatmap_df, cmap="RdYlGn", annot=False, cbar=False, linewidths=0, yticklabels=True)
    except ValueError:
        continue
    # Adjust y-axis labels to be smaller, without rotation and with proper alignment
    plt.yticks(fontsize=2)  # Adjust font size for y labels

    # Add title for the heatmap
    plt.title(f"Heatmap of log2 Fold Change Threshold {threshold} (Red = False, Green = True)")

    # Save the heatmap as a PDF
    plt.savefig(f'./figures/HEATMAP_{threshold}_P.pdf')

    # Close the plot to avoid overlap
    plt.close()














filtered_motifsJA = {
    key: df.with_columns([pl.lit(key).alias("provenance")])  # Add provenance column
    for key, df in motifsJA.items()
}


# Concatenate all filtered DataFrames
combined_dfJA = pl.concat(
    [
        df.select(["id", "log2(fold change)", "adjusted p-value", "provenance"])
        for df in filtered_motifsJA.values()
    ]
)

# Pivot the DataFrame to have separate columns for each provenance
pivot_dfJA = combined_dfJA.pivot(
    values=["log2(fold change)", "adjusted p-value"], 
    index="id", 
    on="provenance"
)




filtered_motifsHC = {
    key: df.with_columns([pl.lit(key).alias("provenance")])  # Add provenance column
    for key, df in motifsHC.items()
}


# Concatenate all filtered DataFrames
combined_dfHC = pl.concat(
    [
        df.select(["id", "log2(fold change)", "adjusted p-value", "provenance"])
        for df in filtered_motifsHC.values()
    ]
)

# Pivot the DataFrame to have separate columns for each provenance
pivot_dfHC = combined_dfHC.pivot(
    values=["log2(fold change)", "adjusted p-value"], 
    index="id", 
    on="provenance"
)




pivot_dfHC = pivot_dfHC.with_columns(
    pl.col("id").str.split(".").list.first().alias("id_group")
)

pivot_dfHC = pivot_dfHC.with_columns(
    pl.when(pl.col("log2(fold change)_I_CD28_positive_fc").abs() > pl.col("log2(fold change)_I_CD28_negative_fc").abs())
    .then(pl.col("log2(fold change)_I_CD28_positive_fc").abs())
    .otherwise(pl.col("log2(fold change)_I_CD28_negative_fc").abs())
    .alias("max_abs_log2_fold_change")
)


# Sort the DataFrame first
sorted_df = pivot_dfHC.sort("max_abs_log2_fold_change", descending=True)

# Use unique to get the first occurrence of each group
result = sorted_df.unique(subset="id_group", keep="first")

result = result.drop("abs_log2_fold_change")

print(result)














####### NETWORK STUFF
#lists of regions 
regions_signi_negative_fc
regions_signi_positive_fc
#gtf
gtf = "/mnt/ndata/daniele/wouter/Scripts/Tools/refdata-cellranger-arc-mm10-2020-A-2.0.0/genes/genes.gtf.gz"

adata_ATAC_CellTypeSubsetTP
net = snap.tl.init_network_from_annotation(regions_signi_positive_fc, gtf, upstream=250000, downstream=250000, id_type='gene_id', coding_gene_only=True)

net_tf = snap.tl.add_tf_binding(network=net, motifs=snap.read_motifs(TFmotifsHC), genome_fasta=snap.genome.mm10, pvalue=0.01)
net_tf_gene = snap.tl.link_tf_to_gene(net)

snap.pl.network_edge_stat(net_tf_gene,show=False,out_file=f"./figures/network.pdf")

# Function to print the graph
def print_graph(graph):
    print("Nodes:")
    for node in graph.node_indices():
        print(f"Node {node}: {graph[node]}")
    
    print("\nEdges:")
    for edge in graph.edge_indices():
        source, target = graph[edge]
        print(edge)

# Print the graph
print_graph(net)


genemat = snap.read("Wouter21_SING_GeneMatNOMAGIC.h5ad", backed=None)
genemat = genemat[np.isin(genemat.obs.index,adata_ATAC_CellTypeSubsetTP.obs.index)]
net_cor = snap.tl.add_cor_scores(net, gene_mat=genemat, peak_mat=adata_ATAC_CellTypeSubsetTP, overwrite=False)

import matplotlib.pyplot as plt
import networkx as nx
import rustworkx as rx

# Convert rustworkx graph to a NetworkX graph for visualization
nx_graph = nx.DiGraph(net_tf_gene)

# Draw the graph using Matplotlib
plt.figure(figsize=(8, 8))
pos = nx.spring_layout(nx_graph)  # Compute positions for nodes
nx.draw(nx_graph, pos, with_labels=True, node_color="skyblue", font_size=10, font_weight="bold", edge_color="gray")
plt.savefig(f"./figures/network.pdf")

# Assuming you have a `rustworkx.PyDiGraph` object named `graph`
pos = rx.spring_layout(net_tf_gene)

# Draw the graph
fig, ax = plt.subplots(figsize=(10, 10))
rx.plot(graph, pos, ax=ax, with_labels=True, node_color='lightblue', edge_color='gray', font_size=10)
plt.savefig(f"./figures/network.pdf")

from rustworkx.visualization import graphviz_draw
graphviz_draw(net_tf_gene, filename=f"./figures/network.pdf")












