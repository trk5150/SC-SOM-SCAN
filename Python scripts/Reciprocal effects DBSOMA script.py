# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 23:40:50 2024

@author: Tim
"""

import pandas as pd
# Load the Excel file
file_path = 'C:/Users/Tim/Desktop/Gut_Structured_Analysis.xlsx'  # Replace with your file path
df = pd.read_excel(file_path)

# Step 1: Sort by "down-up product"
df_sorted = df.sort_values(by="down-up product", ascending=True)


# Define the number of top rows to consider
n = 500  # Replace with your desired number

# Step 2: Iterate over the top `n` rows
top_rows = df_sorted.head(n)

# Step 3: Find the first matching instance from the bottom for each "drug" and "line" combination
results = []
for _, row in top_rows.iterrows():
    drug = row["drug"]
    line = row["Line"]
    initial_row_down_up = row["up-down"]
    initial_row_value = row["down-up product"]
    
    # Find the first matching instance from the bottom
    matching_row = df_sorted[(df_sorted["drug"] == drug) & (df_sorted["Line"] == line)].iloc[-1]
    matching_row_index = df_sorted[(df_sorted["drug"] == drug) & (df_sorted["Line"] == line)].index[-1]
    
    # Collect data for the result
    results.append({
        "Drug": drug,
        "Line": line,
        "Initial_Row_Down_Up": initial_row_down_up,
        "Initial_Row_Value": initial_row_value,
        "Found_Row_Down_Up": matching_row["up-down"],
        "Found_Row_Value": matching_row["down-up product"],
        "Matching_Row_Number": matching_row_index
    })

# Convert results to a DataFrame for better visualization
results_df = pd.DataFrame(results)

# Save the results to an Excel file (optional)
output_file = 'C:/Users/Tim/Desktop/Gut_Matched results.xlsx'
results_df.to_excel(output_file, index=False)