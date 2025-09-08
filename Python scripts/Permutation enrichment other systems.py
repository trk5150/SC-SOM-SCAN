import pandas as pd
import re

# Load the Excel file and specific sheet
#df = pd.read_excel("C://Users//tik105//Desktop//BigScanTest//Genes in other systems//UC//GutUCSOM//Structured_Analysis.xlsx", sheet_name="Dextran up")
#df = pd.read_excel("C://Users//tik105//Desktop//BigScanTest//Genes in other systems//UC//GutUCSOM//Structured_Analysis.xlsx", sheet_name="Dextran down")
#df = pd.read_excel("C://Users//tik105//Desktop//BigScanTest//Genes in other systems//Liver//Liver_SOM//Structured_Analysis.xlsx", sheet_name="FOXA2")
#df = pd.read_excel("C://Users//tik105//Desktop//BigScanTest//Genes in other systems//Liver//Liver_SOM//Structured_Analysis.xlsx", sheet_name="NAFLD")
#df = pd.read_excel("C://Users//tik105//Desktop//BigScanTest//Genes in other systems//S1//S1_SOM//Structured_Analysis.xlsx", sheet_name="DE")
#df = pd.read_excel("C://Users//tik105//Desktop//BigScanTest//Genes in other systems//S4//S4_UT_SOM//Structured_Analysis.xlsx", sheet_name="endocrine")
#df = pd.read_excel("D://Skin single cell//Drug and Gene Structured outputs.xlsx", sheet_name="immune genes")
df = pd.read_excel("D://Cystic Fibrosis single cell//Gene pert out//Name_cleaned_CF_SOM//Structured_Analysis.xlsx", sheet_name="CF down")
# Extract gene names
def extract_gene_name(filename):
    match = re.search(r'_([^_ ]+)(?: |$)', filename)
    return match.group(1) if match else None

df['gene'] = df['fileName'].apply(extract_gene_name)
# Drop missing values and get unique genes
foreground_genes = df['gene'].dropna().unique().tolist()
print("foreground read")

#!!! Load the background gene counts TSV file
background_df = pd.read_csv(
    "D://LINCs sets//GO background//Crispr_and_OE_combined_gene_counts.tsv",
    sep="\t",
    header=0,  # If there's no header row in the file
    names=["gene", "test_count"]  # Rename columns
)

# Drop rows with missing or zero test counts
background_df = background_df.dropna(subset=["gene", "test_count"])
background_df = background_df[background_df["test_count"] > 0]

# Normalize test counts into sampling weights
background_df["weight"] = background_df["test_count"] / background_df["test_count"].sum()

print("background read")

#!!!get a all go-terms associated with the background genes, no p-val cutoff
from gprofiler import GProfiler

# Use all background genes to fetch all relevant GO terms
gp = GProfiler(return_dataframe=True)

# Query gProfiler using your tested gene universe
gene_universe = background_df["gene"].tolist()

gprof_result = gp.convert(organism="hsapiens", query=gene_universe)
# (Optional: check ID formats if needed, but we’ll proceed assuming symbols work.)

# Now get enrichment results — with a loose threshold so we capture many sets
go_results = gp.profile(
    organism="hsapiens",
    query=gene_universe,
    sources=["GO:BP"],
    user_threshold=1.0,         # Keep all GO terms
    no_evidences=False,         # Include gene list in results
)

#!!!turn the go results into a dictionary for permutation analysis
# Build gene set dictionary: term name -> list of genes
from collections import defaultdict

gene_sets = defaultdict(set)

for _, row in go_results.iterrows():
    term_name = row["name"]
    genes_in_term = row["intersections"]  # Already a list
    gene_sets[term_name].update(g.strip() for g in genes_in_term)

# Convert to plain dict with lists
gene_sets = {k: list(v) for k, v in gene_sets.items()}
# Prune gene sets with 10 or fewer assigned genes
gene_sets = {k: v for k, v in gene_sets.items() if len(v) > 10}
print(f"Retained {len(gene_sets)} gene sets with >10 genes.")
print("GO dictionary built")
#!!!permutations from weighted background distribution
import numpy as np
from collections import defaultdict
from statsmodels.stats.multitest import multipletests

# Prepare for permutation test
n_permutations = 10000
foreground_size = len(foreground_genes)
gene_pool = background_df["gene"].tolist()
weights = background_df["weight"].values

# Track null hits per gene set
null_counts = defaultdict(list)

# Pre-convert sets for speed
gene_sets_sets = {term: set(genes) for term, genes in gene_sets.items()}
foreground_set = set(foreground_genes)

# Compute observed hits
observed_hits = {
    term: len(foreground_set & term_genes)
    for term, term_genes in gene_sets_sets.items()
}

# Permutation loop
for _ in range(n_permutations):
    sampled = np.random.choice(gene_pool, size=foreground_size, replace=False, p=weights)
    sampled_set = set(sampled)
    for term, term_genes in gene_sets_sets.items():
        null_counts[term].append(len(sampled_set & term_genes))

# Empirical p-values
p_values = {}
for term, null_dist in null_counts.items():
    obs = observed_hits.get(term, 0)
    p = sum(k >= obs for k in null_dist) / n_permutations
    p_values[term] = p

# FDR correction
terms = list(p_values.keys())
# Calculate expected hits and enrichment ratios
expected_hits = {term: np.mean(null_counts[term]) for term in terms}
enrichment_ratio = {
    term: (observed_hits[term] / expected_hits[term]) if expected_hits[term] > 0 else np.nan
    for term in terms
}

# Optional: Pre-filter weak terms before FDR
filtered_terms = [t for t in terms if observed_hits[t] >= 2 and p_values[t] < 0.05]
raw_p_filtered = [p_values[t] for t in filtered_terms]
_, adj_p_fdr, _, _ = multipletests(raw_p_filtered, method='fdr_bh')

# Map filtered FDR values back to full list (default to 1.0 if not included)
adj_p_map = {t: 1.0 for t in terms}
for t, adj in zip(filtered_terms, adj_p_fdr):
    adj_p_map[t] = adj

# Final enrichment table with extras
enrichment_df = pd.DataFrame({
    "GO term": terms,
    "term_size": [len(gene_sets[t]) for t in terms],
    "observed_hits": [observed_hits[t] for t in terms],
    "expected_hits": [expected_hits[t] for t in terms],
    "enrichment_ratio": [enrichment_ratio[t] for t in terms],
    "p_value": [p_values[t] for t in terms],
    "adj_p_value": [adj_p_map[t] for t in terms],
    "matched_genes": [
        ", ".join(sorted(set(foreground_genes) & set(gene_sets[t])))
        for t in terms
    ]
}).sort_values("p_value")

# Output preview
print(enrichment_df.head())

# Save results
enrichment_df.to_excel("D://Cystic Fibrosis single cell//Gene GO//CF_down_GO_enrichment.xlsx", index=False)





