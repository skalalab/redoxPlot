# partition a df 
import pandas as pd

cell_line = ['IMR90']
df = pd.read_csv("lifetime_noint_df.csv")
cell_lines_str = ', '.join(cell_line)
df = df[df['cell_line'].isin(cell_line)]
df.to_csv(f"{cell_lines_str}_lifetime_noint_df.csv", index=False)
df = pd.read_csv("tmre_df.csv")
df = df[df['cell_line'].isin(cell_line)]
df.to_csv(f"{cell_lines_str}_tmre_df.csv", index=False)