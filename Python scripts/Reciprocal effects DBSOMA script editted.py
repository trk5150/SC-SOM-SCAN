# -*- coding: utf-8 -*-
"""
Pipeline:
- Load Excel
- Parse first column: split on "_" and spaces, strip ".txt"
- Standardize direction to {"up","down"}; drop rows w/o clear direction
- Extract: cell_line = 2nd token, drug = 5th token
- Sort by 'down prod'
- Identify matched pairs: (drug, cell_line) in top 500 'down prod' (direction=down)
  that also appear in top 500 'up prod' (direction=up)
- Save results to Excel with multiple sheets
"""

import re
import pandas as pd
from pathlib import Path

# ==== USER SETTINGS ====
FILE_PATH   = r"E:\Cystic Fibrosis single cell\Drug pert out\Name_cleaned_CF_SOM\CF_Structured_Analysis.xlsx"
SHEET_NAME  = 0          # or a sheet name like "Sheet1"
TOP_N       = 500
OUTPUT_XLSX = r"C:\Users\trk51\OneDrive - Harvard University\Dissertation\CF_structured_matches.xlsx"
# ========================

# ------------- Helpers -------------
def normalize_col_lookup(df, target_name):
    """
    Find a column in df that matches target_name case-insensitively and
    ignoring spaces/underscores/hyphens.
    """
    want = re.sub(r'[\s_\-]+', '', target_name).lower()
    for c in df.columns:
        key = re.sub(r'[\s_\-]+', '', str(c)).lower()
        if key == want:
            return c
    return None

def parse_name_cell(raw):
    """
    Example input:
        'PBIOA020_HT29_24H_G10_BRD-K55630925_0.37uM do.txt'
    - strip trailing .txt
    - split on '_' and ' '
    - interpret final token as direction: {u,up} -> up, {d,dn,down} -> down
    - extract cell_line (2nd token), drug (5th token)
    - return dict or None if direction missing
    """
    if not isinstance(raw, str):
        return None

    base = re.sub(r'\.txt$', '', raw.strip(), flags=re.IGNORECASE)
    parts = [p for p in re.split(r'[ _]+', base) if p]  # split on _ and space

    if not parts:
        return None

    last = parts[-1].lower()
    dir_map = {
        'u': 'up', 'up': 'up',
        'd': 'down', 'dn': 'down', 'down': 'down'
    }
    direction = dir_map.get(last)
    if direction is None:
        # Per spec: drop row if no clear direction token
        return None

    # Extract cell line (2nd token) and drug (5th token) if present
    cell_line = parts[1] if len(parts) > 1 else None
    drug      = parts[4] if len(parts) > 4 else None

    return {
        'parsed_parts': parts,
        'direction': direction,
        'cell_line': cell_line,
        'drug': drug
    }

# ------------- Load data -------------
df = pd.read_excel(FILE_PATH, sheet_name=SHEET_NAME)

# Identify required columns
first_col = df.columns[0]  # the "fileName ..." column to parse

down_col = normalize_col_lookup(df, 'down prod')
up_col   = normalize_col_lookup(df, 'up prod')

if down_col is None:
    raise ValueError("Couldn't find a column matching 'down prod' (case/space-insensitive).")
if up_col is None:
    raise ValueError("Couldn't find a column matching 'up prod' (case/space-insensitive).")

# ------------- Parse & clean -------------
parsed = df[first_col].apply(parse_name_cell)
parsed_df = pd.DataFrame([p for p in parsed if isinstance(p, dict)])

# Keep only rows that parsed cleanly
mask_keep = parsed.notna()
clean = df.loc[mask_keep].copy().reset_index(drop=True)
clean = pd.concat([clean.reset_index(drop=True), parsed_df.reset_index(drop=True)], axis=1)

# Drop rows lacking direction/cell_line/drug post-parse
clean = clean.dropna(subset=['direction', 'cell_line', 'drug'])

# ------------- Sort by 'down prod' -------------
# Interpreting "top 500" as highest values → descending.
clean_sorted_by_down = clean.sort_values(by=down_col, ascending=False).reset_index(drop=True)

# ------------- Top sets -------------
top_down = (clean_sorted_by_down
            .sort_values(by=down_col, ascending=False)
            .head(TOP_N)
            .copy())

top_up = (clean
          .sort_values(by=up_col, ascending=False)
          .head(TOP_N)
          .copy())

# ------------- Opposite-direction matches -------------
# Same (drug, cell_line) appears in top_down (direction=down) and top_up (direction=up)
matches = pd.merge(
    top_down[['drug', 'cell_line', down_col]].rename(columns={down_col: 'down_prod'}),
    top_up[['drug', 'cell_line', up_col]].rename(columns={up_col: 'up_prod'}),
    on=['drug', 'cell_line'],
    how='inner',
)

# Add ranks for context within each top list
top_down = top_down.assign(down_rank=top_down[down_col].rank(method='first', ascending=False))
top_up   = top_up.assign(up_rank=top_up[up_col].rank(method='first', ascending=False))

matches = (matches
           .merge(top_down[['drug','cell_line','down_rank']], on=['drug','cell_line'], how='left')
           .merge(top_up[['drug','cell_line','up_rank']],     on=['drug','cell_line'], how='left')
           .sort_values(['down_rank','up_rank'])
           .reset_index(drop=True))

# ------------- Save outputs -------------
output_path = Path(OUTPUT_XLSX)
with pd.ExcelWriter(output_path, engine='xlsxwriter') as xw:
    clean_sorted_by_down.to_excel(xw, index=False, sheet_name='clean_sorted_by_down')
    top_down.to_excel(xw, index=False, sheet_name=f'top{TOP_N}_down')
    top_up.to_excel(xw, index=False, sheet_name=f'top{TOP_N}_up')
    matches.to_excel(xw, index=False, sheet_name='matched_pairs')

print(f"Saved results to: {output_path}")
print("Sheets:")
print("- clean_sorted_by_down")
print(f"- top{TOP_N}_down")
print(f"- top{TOP_N}_up")
print("- matched_pairs")
