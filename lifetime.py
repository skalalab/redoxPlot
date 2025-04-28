import pandas as pd
import numpy as np

# Define the user-specified columns to always include
cols_to_include = ["base_name", "cell_type", "cell_line", "treatment"]
csv_path = '/Users/allan/GitHub/redoxPlot/250120_cytoplasm_all-QC-PCA-Filter.csv'
filtered_selected_df = pd.DataFrame() # Initialize an empty DataFrame

print(f"Attempting to read: {csv_path}")

try:
    # Read the full CSV first to identify columns and filter
    # Using low_memory=False might help with mixed type warnings if they occur
    df_full = pd.read_csv(csv_path, low_memory=False)
    print(f"Successfully read {csv_path}. Shape: {df_full.shape}")

    # Filter rows where 'experiment' column is "TMRE"
    if 'experiment' in df_full.columns:
        df_filtered_rows = df_full[df_full['experiment'] == 'TMRE'].copy()
        if df_filtered_rows.empty:
             print(f"Warning: No rows found with experiment == 'TMRE' in {csv_path}.")
        else:
             print(f"Filtered rows by experiment == 'TMRE'. Shape: {df_filtered_rows.shape}")
    else:
        print(f"Warning: Column 'experiment' not found in {csv_path}. Cannot filter by experiment.")
        df_filtered_rows = df_full.copy() # Use the full dataframe if 'experiment' column is missing

    if not df_filtered_rows.empty:
        # Identify numerical columns from the filtered DataFrame
        numerical_cols = df_filtered_rows.select_dtypes(include=np.number).columns.tolist()
        print(f"Identified {len(numerical_cols)} numerical columns.")

        # Combine user-specified columns and numerical columns, ensuring uniqueness
        final_cols_to_keep = list(set(cols_to_include + numerical_cols))

        # Ensure the columns actually exist in the filtered DataFrame before selecting
        final_cols_present = [col for col in final_cols_to_keep if col in df_filtered_rows.columns]
        print(f"Will keep {len(final_cols_present)} columns.")

        # Report any requested columns that were not found
        missing_cols = list(set(final_cols_to_keep) - set(final_cols_present))
        if missing_cols:
            print(f"Warning: Columns specified or identified as numerical were not found in the filtered data and will be ignored: {missing_cols}")

        # Select the final columns
        filtered_selected_df = df_filtered_rows[final_cols_present]

        print("\nResulting DataFrame info:")
        print(f"Shape: {filtered_selected_df.shape}")
        print("\nFirst 5 rows:")
        print(filtered_selected_df.head())
        # print("\nColumn names:") # Uncomment to see all column names
        # print(filtered_selected_df.columns.tolist())

        # --- Modifications based on user request ---

        # 1. Remove 'Unnamed: 0' column if it exists
        if 'Unnamed: 0' in filtered_selected_df.columns:
            filtered_selected_df = filtered_selected_df.drop(columns=['Unnamed: 0'])
            print("Removed 'Unnamed: 0' column.")

        # --- Additional modification: Standardize treatment names ---
        if 'treatment' in filtered_selected_df.columns:
            print("Standardizing treatment names...")
            replacements = {
                '0-control': 'control',
                'cccp': 'FCCP',
                'oligomycin': 'Oligo'
            }
            # Check original values before replacement (optional, for verification)
            # print("Original unique treatments:", filtered_selected_df['treatment'].unique())
            filtered_selected_df['treatment'] = filtered_selected_df['treatment'].replace(replacements)
            # Check values after replacement (optional, for verification)
            # print("Updated unique treatments:", filtered_selected_df['treatment'].unique())
            print("Treatment names standardized.")
        else:
            print("Warning: 'treatment' column not found, skipping standardization.")

        # --- Additional modification: Standardize cell_line names ---
        if 'cell_line' in filtered_selected_df.columns:
            print("Standardizing cell_line names...")
            # Check original values before replacement (optional, for verification)
            # print("Original unique cell_lines:", filtered_selected_df['cell_line'].unique())
            filtered_selected_df['cell_line'] = filtered_selected_df['cell_line'].replace({'HPDE6': 'HPDE'})
            # Check values after replacement (optional, for verification)
            # print("Updated unique cell_lines:", filtered_selected_df['cell_line'].unique())
            print("Cell_line names standardized.")
        else:
            print("Warning: 'cell_line' column not found, skipping standardization.")


        # 2. Save the version including 'nint' and 'fint' (lifetime_df)
        lifetime_df = filtered_selected_df.copy()
        lifetime_df_path = 'lifetime_df.csv'
        lifetime_df.to_csv(lifetime_df_path, index=False)
        print(f"Saved DataFrame with 'nint'/'fint' (if present) to {lifetime_df_path}. Shape: {lifetime_df.shape}")

        # 3. Create and save the version excluding 'nint' and 'fint' (lifetime_noint_df)
        cols_to_drop_for_noint = [col for col in ['nint', 'fint'] if col in filtered_selected_df.columns]
        if cols_to_drop_for_noint:
            lifetime_noint_df = filtered_selected_df.drop(columns=cols_to_drop_for_noint)
            print(f"Removed columns {cols_to_drop_for_noint} for the 'noint' version.")
        else:
            lifetime_noint_df = filtered_selected_df.copy() # No columns to drop
            print("Columns 'nint' and 'fint' not found, 'noint' version is the same.")

        lifetime_noint_df_path = 'lifetime_noint_df.csv'
        lifetime_noint_df.to_csv(lifetime_noint_df_path, index=False)
        print(f"Saved DataFrame without 'nint'/'fint' to {lifetime_noint_df_path}. Shape: {lifetime_noint_df.shape}")


except FileNotFoundError:
    print(f"Error: File not found at {csv_path}")
except Exception as e:
    print(f"An error occurred while reading or processing the CSV: {e}")

# The final DataFrame 'filtered_selected_df' contains the desired data.
# You can now work with this DataFrame. For example:
# print(filtered_selected_df.describe())