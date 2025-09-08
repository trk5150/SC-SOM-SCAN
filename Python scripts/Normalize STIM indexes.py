# -*- coding: utf-8 -*-
"""
Normalize 'Stim index' by DMSO mean within each 'Cell Source' and add log2 transform.
Designed to run directly in Spyder with hard-coded paths.

What it does:
1) Reads your Excel file/sheet
2) Computes DMSO mean Stim index per Cell Source
3) Adds:
   - 'Stim index (normalized)'
   - 'Stim index (log2 normalized)'
4) Writes a new Excel file with:
   - 'Normalized Data' sheet (full table + new columns)
   - 'DMSO Means' sheet (summary per Cell Source)
"""

from pathlib import Path
import pandas as pd
import numpy as np

# =========================
# ===== USER SETTINGS =====
# =========================
INPUT_PATH  = Path(r"C://Users//trk51//OneDrive - Harvard University//Drug treatment analyses//Diff 9, 12 and MM drug treatments//Aggregate Analysis.xlsx")
OUTPUT_PATH = Path(r"C://Users//trk51//OneDrive - Harvard University//Drug treatment analyses//Diff 9, 12 and MM drug treatments//NormalizedStims.xlsx")  # <-- change me
SHEET       = 0   # sheet index or name, e.g., 0 or "Sheet1"
# If your control label varies in case/spacing, leave CONTROL_LABEL = "DMSO"
CONTROL_LABEL = "DMSO"

# Required columns (edit only if your headers differ)
REQUIRED_COLUMNS = [
    "Sample",
    "Cell Source",
    "Plate",
    "Dose",
    "Low glucose",
    "High glucose",
    "KCl",
    "Stim index",
]

# =========================
# ===== MAIN SCRIPT =======
# =========================
def main():
    # --- Read
    if not INPUT_PATH.exists():
        raise FileNotFoundError(f"Input Excel not found:\n{INPUT_PATH}")

    df = pd.read_excel(INPUT_PATH, sheet_name=SHEET)

    # --- Column sanity check
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}\n"
                         f"Available columns: {list(df.columns)}")

    # --- Build a clean 'Sample' for control matching (case/space tolerant)
    sample_clean = (
        df["Sample"]
        .astype(str)
        .str.strip()
        .str.upper()
    )
    control_key = str(CONTROL_LABEL).strip().upper()
    dmso_mask = sample_clean == control_key

    if not dmso_mask.any():
        raise ValueError(
            f"No rows found where Sample == '{CONTROL_LABEL}' "
            "(case-insensitive match)."
        )

    # --- Compute DMSO means per 'Cell Source'
    dmso_df = df.loc[dmso_mask, ["Cell Source", "Stim index"]].copy()
    dmso_means = (
        dmso_df.groupby("Cell Source", dropna=False)["Stim index"]
        .mean()
        .rename("DMSO mean Stim index")
    )

    # --- Merge means back to all rows by 'Cell Source'
    out = df.copy()
    out = out.merge(
        dmso_means.reset_index(),
        on="Cell Source",
        how="left",
        validate="many_to_one"
    )

    # --- Create normalized and log2-normalized columns
    out["Stim index (normalized)"] = out["Stim index"] / out["DMSO mean Stim index"]

    # log2 only for positive normalized values; else NaN
    norm = out["Stim index (normalized)"]
    out["Stim index (log2 normalized)"] = np.where(norm > 0, np.log2(norm), np.nan)

    # --- Optional: info on Cell Sources without a control
    no_control_sources = out.loc[out["DMSO mean Stim index"].isna(), "Cell Source"].unique()
    if len([s for s in no_control_sources if pd.notna(s)]) > 0:
        print("WARNING: No control rows found for these Cell Source(s); "
              "normalized values will be NaN:")
        for s in no_control_sources:
            print(f"  - {s}")

    # --- Write output Excel
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(OUTPUT_PATH, engine="openpyxl") as writer:
        out.to_excel(writer, index=False, sheet_name="Normalized Data")
        dmso_means.to_frame().to_excel(writer, sheet_name="DMSO Means")

    print(f"Done. Wrote:\n{OUTPUT_PATH}")

if __name__ == "__main__":
    main()
