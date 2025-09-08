# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 10:43:30 2025

@author: tik105
"""

import pandas as pd
import gseapy as gp
import numpy as np

# Load the gene expression data from Excel
excel_file = "C://Users//tik105//OneDrive - Harvard University//SOMSCAN draft//Other related files//Full_Results_One_vs_All.xlsx"  # Update with your file path
df = pd.read_excel(excel_file)

# Load custom gene set from a line-delimited text file
custom_gene_set_file = "C://Users//tik105//OneDrive - Harvard University//SOMSCAN draft//Other related files//MaturationGeneList.txt"  # Update with your file path
with open(custom_gene_set_file, "r") as f:
    custom_gene_list = [line.strip() for line in f if line.strip()]

# Convert list into a dictionary for gseapy
custom_gene_set = {"Custom_Set": custom_gene_list}

# Extract column headers to find all conditions
conditions = [col.replace("Log2FoldChange_", "").replace("_vs_All", "") 
              for col in df.columns if "Log2FoldChange_" in col]

# Output directory for results
output_dir = "C://Users//tik105//OneDrive - Harvard University//SOMSCAN draft//Other related files//Screen SC-B treated genes//GSEA-like results"

# Process each condition
for condition in conditions:
    print(f"Processing condition: {condition}")

    # Extract relevant columns
    logFC_col = f"Log2FoldChange_{condition}_vs_All"
    pval_col = f"PValue_{condition}_vs_All"

    # Check if columns exist
    if logFC_col not in df.columns or pval_col not in df.columns:
        print(f"Skipping {condition}, required columns missing.")
        continue

    # Compute ranking score: Signed -log10(P-value) * Log2FC
    df["Ranking_Score"] = -np.log10(df[pval_col]) * np.sign(df[logFC_col])

    # Drop NaNs and rank genes
    ranked_df = df[["Gene_Name.Gene_Name", "Ranking_Score"]].dropna()
    ranked_df.columns = ["Gene", "Score"]
    ranked_df = ranked_df.sort_values(by="Score", ascending=False)

    # Save ranked list to a file
    ranked_file = f"{output_dir}/ranked_{condition}.rnk"
    ranked_df.to_csv(ranked_file, sep="\t", index=False, header=False)

    # Run GSEA-like analysis
    gsea_results = gp.prerank(
        rnk=ranked_file, 
        gene_sets=custom_gene_set, 
        outdir=f"{output_dir}/{condition}",
        min_size=5,  # Minimum genes in set
        max_size=5000,  # Maximum genes in set
        permutation_num=1000,  # Number of permutations
        no_plot=True  # Disable plotting
    )

