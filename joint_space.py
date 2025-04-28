import numpy as np
from sklearn.preprocessing import StandardScaler
from SCOT.src.scotv1 import * 
from SCOT.src.scotv2 import *

def create_latent_space(df1, df2, method):
    N_FACTORS = 20
    # assign the numeric columns to X1 and X2
    X1 = df1.select_dtypes(include=[np.number]).values.astype(np.float32)
    X2 = df2.select_dtypes(include=[np.number]).values.astype(np.float32)

    categorical_cols_df1 = df1.select_dtypes(include=['object'])
    categorical_cols_df2 = df2.select_dtypes(include=['object'])
    ids_df1 = df1["cell_id"].tolist()
    ids_df2 = df2["cell_id"].tolist()
    # scale each dataset separately (z‑score per feature)
    X1 = StandardScaler().fit_transform(X1)
    X2 = StandardScaler().fit_transform(X2)
    labels_stacked  = np.array(["Lifetime"] * len(X1) + ["TMRE"] * len(X2))
    ids_stacked    = np.array(ids_df1 + ids_df2)
    # stack the categorical columns
    categorical_cols_stacked = np.vstack((categorical_cols_df1.values, categorical_cols_df2.values))
    SCOT=SCOTv2([X1,X2])
    aligned_datasets=SCOT.align(normalize=True, k=50, eps=0.005, rho=0.1, projMethod="embedding",)
    print(len(aligned_datasets))
    print(aligned_datasets[0].shape)
    print(aligned_datasets[1].shape)
    pass