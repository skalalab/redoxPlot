#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 24 09:56:42 2025

@author: allan
"""

#Base Libraries
import os
import shutil

#Core Functions
#%%
import cemia55s as cemia
rt = "/home/useradmin/Downloads/tmre_redox"
cell_lines = ["A549", "HeLa", "HFF", "HPDE", "IMR90", "MCF10A", "MCF7", "Panc1"]
treatment_dir = ["A1", "A2", "A3"]

#%%
# dfs = []
# for cell_line in cell_lines:
#     for treatment in treatment_dir:
#         treatment_name = "FCCP" if treatment == "A1" else "control" if treatment == "A2" else "Oligo"
#         treatment_path = os.path.join(rt, cell_line, treatment)
#         # Get the list of absolute paths that has _MitoMask_ in the treatment directory
#         cell_list = [
#             os.path.join(root, file)
#             for root, _, files in os.walk(treatment_path)
#             for file in files
#             if "_MitoMask_" in file
#         ]
#         print(cell_list[0])
#         df = cemia.measurement(treatment_path, cell_list, "")
#         df['cell_line'] = cell_line
#         df['treatment'] = treatment_name
#         dfs.append(df)

# # Concatenate all dataframes into one
# import pandas as pd
# mito_morph_df = pd.concat(dfs, ignore_index=True)

# # Save the dataframe to a CSV file
# output_path = os.path.join(rt, "mito_morph.csv")
# mito_morph_df.to_csv(output_path, index=False)

#%%
import pandas as pd
df = pd.read_csv(os.path.join(rt, "mito_morph.csv"))
df['cell_id'] = df.apply(
    lambda row: f"{row['cell_line']}_{row['treatment']}_{row['cell_name'].split('/')[-1].replace('MitoMask_', '').replace('.tiff', '')}", 
    axis=1
)

# Get the number of columns before dropping
original_shape = df.shape
print(f"Original dataframe shape: {original_shape}")

# Drop columns with 'std' or 'median' in their names
columns_to_drop = [col for col in df.columns if 'std' in col or 'median' in col]
print(f"Number of columns to drop: {len(columns_to_drop)}")
print(f"Columns being dropped: {columns_to_drop}")

# Drop the columns
df = df.drop(columns=columns_to_drop)

# Print the shape of the resulting dataframe
print(f"New dataframe shape: {df.shape}")
# %%
import pandas as pd
mito_df = pd.read_csv(os.path.join(rt, "tmre_analysis.csv"))
# Step 2: Create the cell_id column
def extract_file_suffix(file_name):
    """
    Extract the suffix from file name in format like "A1_N2_MitoI_12.tiff"
    to get "A1_N2_12"
    """
    parts = file_name.split('_')
    if len(parts) >= 4:
        prefix = parts[0]  # A1
        middle = parts[1]  # N2
        last_part = parts[-1]  # 12.tiff
        number = last_part.split('.')[0]  # 12
        return f"{prefix}_{middle}_{number}"
    return file_name

# Create the cell_id column combining cell_line, treatment, and file suffix
mito_df['cell_id'] = mito_df.apply(
    lambda row: f"{row['cell_line']}_{row['treatment']}_{extract_file_suffix(row['file_name'])}", 
    axis=1
)

# %%
merged_df = pd.merge(
    mito_df, 
    df.drop(['cell_line', 'treatment'], axis=1),  # Remove duplicate columns
    on='cell_id',
    how='inner'  # This will keep only matching rows
)

# remove columns with na values
merged_df = merged_df.dropna(axis=1, how='any')
#remove file_name column and cell_name column
merged_df = merged_df.drop(['file_name', 'cell_name'], axis=1)
print(merged_df.shape)
# Save the merged dataframe to a CSV file
output_path = os.path.join(rt, "merged_mito_df.csv")
merged_df.to_csv(output_path, index=False)
# %%
