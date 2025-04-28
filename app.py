import pandas as pd
import scipy.stats as stats
import numpy as np
from plotly.subplots import make_subplots
from plotly import graph_objects as go
import seaborn as sns
from joint_space import create_latent_space
mito_df = pd.read_csv('tmre_df.csv')
mito_lifetime_df = pd.read_csv('tmre_lifetime.csv')

lifetime_features = ['na1', 'na2', 'nt1', 'nt2', 'ntm', 'fa1',
                     'fa2', 'ft1', 'ft2', 'ftm', 'nint', 'fint', 'normrr']
mito_features = ['mito_area', 'mito_intensity_sum', 'mito_intensity_avg', 'nuclei_area',
                'nuclei_intensity_sum', 'nuclei_intensity_avg', 'Delta Psi 1',
                'cell_area', 'mito_content']
mito_lifetime_features = lifetime_features + mito_features
from filter_widgets import filters_widget
import streamlit as st
st.set_page_config(layout="wide", initial_sidebar_state="collapsed")
col1, col2 = st.columns([1, 2.5])
with col1:
    method = st.selectbox(  
        "Select a visualization method", ["Mito Features Correlation", "Heatmap", "Lifetime-Mito Features Correlation", "Joint UMAP"])
    
    features = mito_lifetime_features if "Lifetime" in method else mito_features
    df = mito_df
    if method == "Mito Features Correlation" or method == "Lifetime-Mito Features Correlation":
        st.write("Select features to correlate")

        if method == "Lifetime-Mito Features Correlation":
            # filter out the lifetime features
            df = mito_lifetime_df
            features = mito_lifetime_features
    
        # create two single select widgets for feature selection: Select X and Select Y
        # default to be the first feature in the list
        selected_x = st.selectbox("Select X", features, index=0)
        # default to be the first feature not selected in X
        selected_y = st.selectbox("Select Y", [feature for feature in features if feature != selected_x], index=0)
    elif method == "Heatmap":
        pass
    elif method == "Joint UMAP":
        # upload two csv files
        tmre_df = st.file_uploader("Upload tmre feature dataset", type="csv")
        tmre_df = pd.read_csv(tmre_df) if tmre_df is not None else None
        lifetime_df = st.file_uploader("Upload lifetime feature dataset", type="csv")
        lifetime_df = pd.read_csv(lifetime_df) if lifetime_df is not None else None
    #     col1, col2 = st.columns(2)
    # # First number incrementor in the first column
    #     with col1:
        n_neighbors = st.number_input(
            "n_neighbors",
            value=15,  # Initial value
            step=5,             # Increment/Decrement step
            format="%d"            # Integer format
        )
        # Second number incrementor in the second column
        # with col2:
        min_dist = st.number_input(
            "min_dist",
            value=0.1,  # Initial value
            step=0.1,            
        )
        latent_space_method = st.selectbox(
            "Select a latent space method", [ "SCOT"])
with col2:
    if method == "Joint UMAP": 
        if tmre_df is not None and lifetime_df is not None:
            filtered_tmre_df, color_by_options = filters_widget(tmre_df, wildcard=True)
            available_cell_types = filtered_tmre_df["cell_type"].unique()
            available_cell_lines = filtered_tmre_df["cell_line"].unique()
            available_treatments = filtered_tmre_df["treatment"].unique()
            # use this to filter the lifetime_df
            filtered_lifetime_df = lifetime_df[
                (lifetime_df["cell_type"].isin(available_cell_types)) &
                (lifetime_df["cell_line"].isin(available_cell_lines)) &
                (lifetime_df["treatment"].isin(available_treatments))
            ]
            # create a joint latent space
            latent_space_matrix = create_latent_space(tmre_df, lifetime_df, method=latent_space_method)
    else: 
        filtered_df, color_by_options = filters_widget(df, wildcard=True)
    if method == "Mito Features Correlation":
        if not filtered_df.empty and selected_x and selected_y:
            # create a scatter plot of the selected features with correlation coefficient and p-value
            
            # calculate correlation coefficient and p-value
            corr_coef, p_value = stats.pearsonr(filtered_df[selected_x], filtered_df[selected_y])
            # create scatter plot with color_by_options
            filtered_df['unique_color_group'] = filtered_df[color_by_options].agg('_'.join, axis=1)
            unique_color_groups = filtered_df['unique_color_group'].unique()
            alpha = 0.6 if len(unique_color_groups) > 1 else 1.0
            palette = sns.color_palette("tab10", n_colors=len(unique_color_groups))
            color_sequence = [f"rgba({int(color[0]*255)}, {int(color[1]*255)}, {int(color[2]*255)}, {alpha})" for color in palette]
            color_map = {t: color_sequence[i] for i, t in enumerate(unique_color_groups)}

            fig = make_subplots(rows=1, cols=1)
            for group in unique_color_groups:
                group_df = filtered_df[filtered_df['unique_color_group'] == group]
                fig.add_trace(go.Scatter(
                x=group_df[selected_x],
                y=group_df[selected_y],
                mode='markers',
                marker=dict(
                    color=color_map[group],
                    size=5
                ),
                name=group
                ))
            # add correlation coefficient and p-value to the plot
            fig.add_annotation(text=f"Correlation coefficient: {corr_coef:.2f}<br>P-value: {p_value:.2e}", xref="paper", yref="paper", x=0.5, y=1.05, showarrow=False)
            # update layout
            fig.update_layout(title=f"Correlation between {selected_x} and {selected_y}", xaxis_title=selected_x, yaxis_title=selected_y)
            # show plot
            st.plotly_chart(fig, use_container_width=True)
    elif method == "Heatmap":
        if not filtered_df.empty:
            # create a correlation heatmap of all mito features
            import matplotlib.pyplot as plt

            # calculate correlation matrix
            corr_matrix = filtered_df[features].corr()

            # create a heatmap
            fig, ax = plt.subplots(figsize=(10, 8))
            sns.heatmap(corr_matrix, annot=True, fmt=".2f", cmap="coolwarm", cbar=True, ax=ax)

            # set titles and labels
            ax.set_title("Correlation Heatmap", fontsize=16)
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right')
            ax.set_yticklabels(ax.get_yticklabels(), rotation=0)

            # render the heatmap in Streamlit
            st.pyplot(fig)
            plt.close(fig)

    elif method == "Lifetime-Mito Features Correlation":
        if selected_x and selected_y:
            # plot correlation between selected_x and selected_y
            # calculate correlation coefficient and p-value
            # Extract feature names without the 'std_' prefix
            x_base = selected_x
            y_base = selected_y
            
            # Check if std columns exist
            # Add controls for error bars
            col_err1, col_err2 = st.columns(2)
            with col_err1:
                show_x_error = st.checkbox("Show X error bars", value=False)
            with col_err2:
                show_y_error = st.checkbox("Show Y error bars", value=True)
            
            # Check if std columns exist
            x_std_col = f"{selected_x}_std" if f"{selected_x}_std" in mito_lifetime_df.columns and show_x_error else None
            y_std_col = f"{selected_y}_std" if f"{selected_y}_std" in mito_lifetime_df.columns and show_y_error else None
            # add those std columns to filtered_df
            corr_coef, p_value = stats.pearsonr(mito_lifetime_df[selected_x], mito_lifetime_df[selected_y])
            # create scatter plot with color_by_options
            mito_lifetime_df['unique_color_group'] = mito_lifetime_df[["cell_line", "treatment"]].agg('_'.join, axis=1)
            unique_color_groups = mito_lifetime_df['unique_color_group'].unique()
            alpha = 1.0
            palette = sns.color_palette("tab20", n_colors=len(unique_color_groups))
            color_sequence = [f"rgba({int(color[0]*255)}, {int(color[1]*255)}, {int(color[2]*255)}, {alpha})" for color in palette]
            color_map = {t: color_sequence[i] for i, t in enumerate(unique_color_groups)}
            
            fig = make_subplots(rows=1, cols=1)
            for group in unique_color_groups:
                group_df = mito_lifetime_df[mito_lifetime_df['unique_color_group'] == group]
                
                # Add scatter plot
                fig.add_trace(go.Scatter(
                    x=group_df[selected_x],
                    y=group_df[selected_y],
                    mode='markers',
                    marker=dict(
                        color=color_map[group],
                        size=5
                    ),
                    error_x=dict(
                        type='data',
                        array=group_df[x_std_col] if x_std_col else None,
                        visible=True if x_std_col else False,
                        color=color_map[group],
                        thickness=1,
                        width=3
                    ),
                    error_y=dict(
                        type='data',
                        array=group_df[y_std_col] if y_std_col else None,
                        visible=True if y_std_col else False,
                        color=color_map[group],
                        thickness=1,
                        width=3
                    ),
                    name=group
                ))

            # add correlation coefficient and p-value to the plot
            fig.add_annotation(text=f"Correlation coefficient: {corr_coef:.2f}<br>P-value: {p_value:.2e}", xref="paper", yref="paper", x=0.5, y=1.05, showarrow=False)
            # update layout
            fig.update_layout(
                title=f"Correlation between {selected_x} and {selected_y}", 
                xaxis_title=selected_x, 
                yaxis_title=selected_y
            )
            # show plot
            st.plotly_chart(fig, use_container_width=True)
    elif method == "Joint UMAP":
        pass