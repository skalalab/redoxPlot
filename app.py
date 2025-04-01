import pandas as pd
import scipy.stats as stats
import numpy as np
from plotly.subplots import make_subplots
from plotly import graph_objects as go
import seaborn as sns
mito_df = pd.read_csv('merged_mito_df.csv')

mito_features = mito_df.columns.tolist()
# remove columns that are not numerical
mito_features = [col for col in mito_features if pd.api.types.is_numeric_dtype(mito_df[col])]

mito_lifetime_features = [] + mito_features
from filter_widgets import filters_widget
import streamlit as st
st.set_page_config(layout="wide", initial_sidebar_state="collapsed")
col1, col2 = st.columns([1, 2.5])
with col1:
    method = st.selectbox(  
        "Select a visualization method", ["Mito Features Correlation", "Heatmap"])
    
    features = mito_lifetime_features if "Lifetime" in method else mito_features
    df = mito_df
    if method == "Mito Features Correlation":
        st.write("Select features to correlate")
    
        # create two single select widgets for feature selection: Select X and Select Y
        # default to be the first feature in the list
        selected_x = st.selectbox("Select X", features, index=0)
        # default to be the first feature not selected in X
        selected_y = st.selectbox("Select Y", [feature for feature in features if feature != selected_x], index=0)
    elif method == "Heatmap":
        pass
with col2:
    filtered_df, color_by_options, cols = filters_widget(df, wildcard=True)
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
