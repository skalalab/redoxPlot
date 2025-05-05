import numpy as np
from sklearn.preprocessing import StandardScaler
from SCOT.src.scotv1 import * 
from SCOT.src.scotv2 import *
import streamlit as st
import pandas as pd
@st.cache_data
def create_latent_space(df1, df2, method):
    N_FACTORS = 20
    # assign the numeric columns to X1 and X2
    X1 = df1.select_dtypes(include=[np.number]).values.astype(np.float32)
    X2 = df2.select_dtypes(include=[np.number]).values.astype(np.float32)

    categorical_cols_df1 = df1[["cell_id", "cell_type", "cell_line", "treatment"]].astype(str)
    categorical_cols_df2 = df2[["cell_id", "cell_type", "cell_line", "treatment"]].astype(str)

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
    #aligned_datasets=SCOT.align(normalize=True, k=50, eps=0.005, rho=0.1, projMethod="embedding",)
    #aligned_datasets=SCOT.align(normalize=True, k=100, eps=0.005, rho=0.1, projMethod="embedding",)
    aligned_datasets=SCOT.align(normalize=True, k=20, eps=0.005, rho=0.1, projMethod="embedding",)
    print(len(aligned_datasets))
    print(aligned_datasets[0].shape)
    print(aligned_datasets[1].shape)
    # combine the aligned datasets
    combined_data = np.vstack((aligned_datasets[0], aligned_datasets[1]))
    # based on the returned number of columns(.shape[1]), name the columns as Dim1, Dim2, Dim3, etc.
    col_names = [f"Dim{i+1}" for i in range(combined_data.shape[1])]
    # create a new DataFrame with the combined data and the column names
    combined_data = pd.DataFrame(combined_data, columns=col_names)
    # attach the categorical columns
    categorical_cols_stacked = pd.DataFrame(categorical_cols_stacked, columns=["cell_id", "cell_type", "cell_line", "treatment"])
    combined_data = pd.concat([combined_data, categorical_cols_stacked], axis=1)
    # attach the labels
    labels_stacked = labels_stacked.reshape(-1, 1)
    # give the labels a name
    labels_stacked = pd.DataFrame(labels_stacked, columns=["experiment"])
    combined_data = pd.concat([combined_data, labels_stacked], axis=1)
    # save the combined data df as a csv file
    
    combined_data.to_csv("IMR90_embed.csv", index=False)

    pass 