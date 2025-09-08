# -*- coding: utf-8 -*-
"""
Created on Fri Aug 29 05:42:44 2025

@author: trk51
"""

import pandas as pd
import gseapy as gp
import numpy as np
from statsmodels.stats.multitest import multipletests

# Load the gene expression data from Excel
excel_file = "C://Users//trk51//OneDrive - Harvard University//SOMSCAN draft//Other related files//Full_Results_One_vs_All.xlsx"  # Update with your file path
df = pd.read_excel(excel_file)

# Load custom gene set from a line-delimited text file
custom_gene_set_file = "C://Users//trk51//OneDrive - Harvard University//SOMSCAN draft//Other related files//MaturationGeneList.txt"  # Update with your file path
with open(custom_gene_set_file, "r") as f:
    custom_gene_list = [line.strip() for line in f if line.strip()]

# Convert list into a dictionary for gseapy
custom_gene_set = {"Custom_Set": custom_gene_list}

# Extract column headers to find all conditions
conditions = [col.replace("Log2FoldChange_", "").replace("_vs_All", "") 
              for col in df.columns if "Log2FoldChange_" in col]

# Output directory for results
output_dir = "C://Users//trk51//OneDrive - Harvard University//SOMSCAN draft//Other related files//Screen SC-B treated genes//GSEA-like results 2"

all_rows = []  # collect per-condition results for later cross-condition FDR

for condition in conditions:
    print(f"Processing condition: {condition}")
    logFC_col = f"Log2FoldChange_{condition}_vs_All"
    pval_col = f"PValue_{condition}_vs_All"
    if logFC_col not in df.columns or pval_col not in df.columns:
        print(f"Skipping {condition}, required columns missing.")
        continue

    df["Ranking_Score"] = -np.log10(df[pval_col]) * np.sign(df[logFC_col])
    ranked_df = df[["Gene_Name.Gene_Name", "Ranking_Score"]].dropna()
    ranked_df.columns = ["Gene", "Score"]
    ranked_df = ranked_df.sort_values(by="Score", ascending=False)

    ranked_file = f"{output_dir}/ranked_{condition}.rnk"
    ranked_df.to_csv(ranked_file, sep="\t", index=False, header=False)

    gsea = gp.prerank(
        rnk=ranked_file,
        gene_sets=custom_gene_set,     # single set
        outdir=f"{output_dir}/{condition}",
        min_size=5,
        max_size=5000,
        permutation_num=1000,
        no_plot=True
    )

    # gsea.res2d is a DataFrame with the columns you mentioned
    res = gsea.res2d.copy()
    res["Condition"] = condition
    all_rows.append(res)

# Combine results across conditions
if all_rows:
    combined = pd.concat(all_rows, ignore_index=True)

    # Apply cross-condition correction per gene set “Term”
    # (You only have one Term, but this generalizes if you add more later.)
    adjusted_frames = []
    for term, sub in combined.groupby("Term", as_index=False):
        # Use the nominal p-values for correction across conditions
        p = sub["NOM p-val"].to_numpy()
        rej, q_bh, _, _ = multipletests(p, method="fdr_bh")   # Benjamini–Hochberg
        _, q_holm, _, _ = multipletests(p, method="holm")     # Holm (FWER-style)
        sub = sub.assign(FDR_q_across_conditions=q_bh,
                         Holm_p_across_conditions=q_holm,
                         Rejected_BH=rej)
        adjusted_frames.append(sub)
    combined_adj = pd.concat(adjusted_frames, ignore_index=True)

    # Save a tidy table you can sort/filter later
    combined_adj.to_csv(f"{output_dir}/GSEA_across_conditions_adjusted.csv", index=False)
