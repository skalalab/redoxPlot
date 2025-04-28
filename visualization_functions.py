import seaborn as sns
import matplotlib.pyplot as plt
from itertools import combinations
import streamlit as st
import plotly.graph_objects as go
import numpy as np
from sklearn.mixture import GaussianMixture

def glass_delta(group1, group2):
    mean_diff = np.mean(group1) - np.mean(group2)
    group2_sd = np.std(group2, ddof=1)  # Using Bessel's correction with ddof=1
    return mean_diff / group2_sd

def create_color_map(groups, overlap_point):
    # if points in the visulization is going to overlap, use a transparent color
    if overlap_point: 
        alpha = 0.6 if len(groups) > 1 else 1.0
    else:
        alpha = 1.0
    palette = sns.color_palette("tab10", n_colors=len(groups))
    color_sequence = [f"rgba({int(color[0]*255)}, {int(color[1]*255)}, {int(color[2]*255)}, {alpha})" for color in palette]
    color_map = {t: color_sequence[i] for i, t in enumerate(groups)}
    return color_map



def dimension_reduction_plot(df, method="UMAP", colored_by=[], exp_var=None):
    """create a plotly plot to visualize the dimension-reduced data
    """
    fig = go.Figure()
    if method == "Principal Component Analysis":
        axis_labels = ["PC1", "PC2"]
    elif method == "UMAP":
        axis_labels = ["UMAP1", "UMAP2"]
    else:
        axis_labels = ["dim1", "dim2"]

    # colored by unique combinations of the selected categorical columns
    df['unique_color_group'] = df[colored_by].agg('_'.join, axis=1)
    unique_color_groups = df['unique_color_group'].unique()

    color_map = create_color_map(unique_color_groups, overlap_point=True)

    # plot scatter plot iteratively, once for each color group
    for g in unique_color_groups:
        g_df =  df[df['unique_color_group'] == g]
        fig.add_trace(
            go.Scatter(
                x=g_df[axis_labels[0]],
                y=g_df[axis_labels[1]],
                mode='markers',
                name=f'{g}',
                text=g_df["cell_id"],
                customdata=g_df["image_name"],
                hovertemplate="<b>%{text}</b>",
                marker=dict(color=color_map[g])
            ),
    )
        
    fig.update_layout(
        hovermode='closest'
    )

    # Update axis labels to include explained variance
    if exp_var is not None: 
        fig.update_xaxes(title_text=f"{axis_labels[0]}({exp_var[0]:.2f}%)")
        fig.update_yaxes(title_text=f"{axis_labels[1]}({exp_var[1]:.2f}%)")
    else:
        fig.update_xaxes(title_text=f"{axis_labels[0]}")
        fig.update_yaxes(title_text=f"{axis_labels[1]}")
    # remove the column after plotting
    df.drop(columns=['unique_color_group'], inplace=True)
    return fig

