from shiny import App, ui, render, reactive
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from htmltools import HTML
from pathlib import Path
from shinywidgets import output_widget, render_widget  
import plotly.graph_objects as go
import plotly.express as px
import umap
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
import seaborn as sns

df = pd.read_csv(Path(__file__).parent / "240710_cytoplasm_all-QC.csv")
measures = ["na1","na2","nt1","nt2","ntm","fa1","fa2","ft1","ft2","ftm","nint","fint","normrr"]
nadh_measures = ["na1", "na2","nt1","nt2","ntm"]
fadh_measures = ["fa1","fa2","ft1","ft2","ftm"]
grouped_df = pd.read_csv(Path(__file__).parent / "grouped_df.csv")
grouped_df = grouped_df.rename(columns={grouped_df.columns[0]: 'id'})
grouped_df_norm =  pd.read_csv(Path(__file__).parent / "grouped_norm.csv")
volcano_df = pd.read_csv(Path(__file__).parent / "volcano.csv")
column_names = grouped_df.columns.tolist()[1:8]
unique_values_dict = {col: grouped_df[col].unique().tolist() for col in column_names}
experiments_treatment_dict = {name: group["treatment"].unique().tolist() for name, group in df.groupby("experiment")}
experiments_cellline_dict = {name: group["cell_line"].unique().tolist() for name, group in df.groupby("experiment")}
experiments_cancer_dict = {name: group["cancer"].unique().tolist() for name, group in df.groupby("experiment")}
experiments_celltype_dict = {name: group["cell_type"].unique().tolist() for name, group in df.groupby("experiment")}
experiments_media_dict = {name: group["media"].unique().tolist() for name, group in df.groupby("experiment")}
vars_colnames = ["na1_mean", "na2_mean","nt1_mean","nt2_mean","ntm_mean", "fa1_mean","fa2_mean","ft1_mean","ft2_mean","ftm_mean", "nint_mean", "fint_mean", "normrr_mean"]
vars = [name.replace("_mean", "") for name in vars_colnames]
pca_filters = ["cell_line", "treatment", "cancer", "cell_type", "media"]

# generate colors for vars and cell lines
cm = plt.get_cmap('tab20')
var_colors = {var: cm(i) for i, var in enumerate(vars)}
cm = plt.get_cmap('tab10')
cell_line_colors = {cell_line: cm(i) for i, cell_line in enumerate(unique_values_dict[column_names[0]])}

# Define the UI layout
app_ui = ui.page_navbar(
    ui.nav_panel("1", 
        ui.navset_tab(  
            ui.nav_panel("Data Table", ui.layout_sidebar(
                ui.sidebar(
                    ui.layout_columns(ui.input_checkbox_group(column_names[0], "Cell Line",["All"], selected=["All"]),
                        ui.input_checkbox_group(column_names[5], "Cell Type", ["All"] + unique_values_dict[column_names[5]], selected=["All"]), 
                        ui.input_checkbox_group(column_names[6], "Tissue", ["All"] + unique_values_dict[column_names[6]], selected=["All"]),                       
                        col_widths=(4, 4, 4)),
                    ui.layout_columns(ui.input_checkbox_group(column_names[1], "Cancer", ["All"] + unique_values_dict[column_names[1]], selected=["All"]),
                        ui.input_checkbox_group(column_names[2], "Media", ["All"] + unique_values_dict[column_names[2]], selected=["All"]), 
                        col_widths=(5, 7)),
                    ui.layout_columns(
                        ui.input_checkbox_group(column_names[3], "Experiment",["All"] + unique_values_dict[column_names[3]], selected=["All"]),
                        ui.input_checkbox_group(column_names[4], "Treatment", []), 
                        col_widths=(6,6)), width=400, 
                ), HTML("<h2>Dataset Explorer</h2>"),
                   ui.output_data_frame("filtered_table"),
                )),
            ui.nav_panel("Plots", ui.layout_sidebar(
                ui.sidebar(   
                    ui.layout_columns(
                        ui.card(ui.input_selectize("select_experiment", "Select experiment", unique_values_dict[column_names[3]]), height="200px")
                        , col_widths=12),
                    ui.layout_columns(
                        ui.card(ui.input_selectize("select_cell_line", "Select cell line(s)", ["All cell lines"] + unique_values_dict[column_names[0]], multiple=True), height="200px"),
                        ui.output_ui("add_2nd_cellline"), col_widths=(6,6)),
                    ui.layout_columns(
                        ui.card(ui.input_selectize("select_var", "Select measure(s)", ["All measures"] + vars, multiple=True), height="200px"),
                        ui.output_ui("add_2nd_measure"), col_widths=(6,6)),
                        ui.output_ui("turn_on_2nd"), width=400,
                ), ui.card(ui.card_header(HTML("<h3>Plot percent change for treatments against control</h3>")),
                ui.output_plot("plotTreatmentsXcontrol", width="100%", height="60%"), height="600px"),
                )),          
            ui.nav_panel("Volcano Plot", ui.layout_sidebar(
                ui.sidebar(ui.layout_columns(
                        ui.card(ui.input_selectize("select_experiment_v", "Select experiment", unique_values_dict[column_names[3]], multiple=False), height="200px"),
                        ui.card(ui.input_selectize("select_var_v", "Select measure(s)", vars, multiple=False), height="150px"), 
                        ui.card(ui.input_numeric("log2fc_cutoff", "Specify the absolute log2 fold change cutoff", 0.5, min=0, max=5, step=0.1), height="150px"), 
                        col_widths=12),), 
                output_widget('volcano_plot')
                )),
            ui.nav_panel("PCA Plot", ui.layout_sidebar(
                ui.sidebar(ui.layout_columns(
                        ui.card(ui.input_selectize("select_experiment_pca", "Select experiment", unique_values_dict[column_names[3]], multiple=False), height="150px"),
                        ui.card(ui.input_selectize("select_filter_pca", "Filtered by", pca_filters + ["Nothing"], multiple=False), height="100px"),
                        ui.card(ui.input_selectize("filter_pca", None, [] , multiple=False), height="150px"),  
                        ui.card(ui.input_selectize("color_pca", "Colored by", pca_filters, multiple=False), height="150px"),  
                        ui.input_switch("nadh_vs_full", "All vs NADPH", False),
                        ui.input_switch("fad_vs_full", "All vs FAD", False),
                        ui.input_switch("loadings", "Show PC loadings?", False),
                        col_widths=12),), 
                ui.output_plot("plotPCA", width="100%", height="60%"), height="600px")
                ),
            ui.nav_panel("UMAP Plot", ui.layout_sidebar(
                ui.sidebar(ui.layout_columns(
                        ui.card(ui.input_selectize("select_experiment_umap", "Select experiment", unique_values_dict[column_names[3]], multiple=False), height="150px"),
                        ui.card(ui.input_selectize("select_filter_umap", "Filtered by", pca_filters+ ["Nothing"], multiple=False), height="100px"),
                        ui.card(ui.input_selectize("filter_umap", None, [] , multiple=False), height="150px"),  
                        ui.card(ui.input_selectize("color_umap", "Colored by", pca_filters, multiple=False), height="150px"),  
                        ui.layout_columns(ui.input_numeric("n_neighbors", "n_neighbors", 15, min=0, step=10),
                                          ui.input_numeric("min_dist", "min_dist", 0.5, min=0, step=0.1)
                        , col_widths=(6,6)),
                        ui.input_switch("nadh_vs_full_umap", "All vs NADPH", False),
                        ui.input_switch("fad_vs_full_umap", "All vs FAD", False),
                        col_widths=12),), 
                ui.output_plot("plotUMAP", width="100%", height="60%"), height="600px")
                ),
            ui.nav_panel("t-SNE Plot", ui.layout_sidebar(
                ui.sidebar(ui.layout_columns(
                        ui.card(ui.input_selectize("select_experiment_tsne", "Select experiment", unique_values_dict[column_names[3]], multiple=False), height="150px"),
                        ui.card(ui.input_selectize("select_filter_tsne", "Filtered by", pca_filters, multiple=False), height="100px"),
                        ui.card(ui.input_selectize("filter_tsne", None, [] , multiple=False), height="150px"),  
                        ui.card(ui.input_selectize("color_tsne", "Colored by", pca_filters, multiple=False), height="150px"),  
                        ui.input_numeric("perplexity", "Perplexity", 20, min=0, step=10),
                        col_widths=12),), 
                ui.output_plot("plotTSNE", width="100%", height="60%"), height="600px")
                ),
            id="tab",)
    ),     
    ui.nav_panel("2", "Page 2"),  
    ui.nav_panel("3", "Page 3"),  
    title="Redox Plots",  
    id="page",  
    fillable = True, 
)  


# Define the server logic
def server(input, output, session):
    previous_states = {"nadh_vs_full": False, "fad_vs_full": False, "loadings": False}

    @reactive.Effect
    def update_switches():
        nonlocal previous_states

        # Check if `nadh_vs_full` was just turned on
        if input.nadh_vs_full() and not previous_states["nadh_vs_full"]:
            session.send_input_message("fad_vs_full", {"value": False})
            session.send_input_message("loadings", {"value": False})

        # Check if `fad_vs_full` was just turned on
        elif input.fad_vs_full() and not previous_states["fad_vs_full"]:
            session.send_input_message("nadh_vs_full", {"value": False})
            session.send_input_message("loadings", {"value": False})

        # Check if `loadings` was just turned on
        elif input.loadings() and not previous_states["loadings"]:
            session.send_input_message("nadh_vs_full", {"value": False})
            session.send_input_message("fad_vs_full", {"value": False})

        # Update the previous states
        previous_states["nadh_vs_full"] = input.nadh_vs_full()
        previous_states["fad_vs_full"] = input.fad_vs_full()
        previous_states["loadings"] = input.loadings()

    previous_states_umap = {"nadh_vs_full_umap": False, "fad_vs_full_umap": False}
    @reactive.Effect
    def _():
        nonlocal previous_states_umap

        # Check if `nadh_vs_full` was just turned on
        if input.nadh_vs_full_umap() and not previous_states_umap["nadh_vs_full_umap"]:
            session.send_input_message("fad_vs_full_umap", {"value": False})
        # Check if `fad_vs_full` was just turned on
        elif input.fad_vs_full_umap() and not previous_states_umap["fad_vs_full_umap"]:
            session.send_input_message("nadh_vs_full_umap", {"value": False})
        # Update the previous states
        previous_states_umap["nadh_vs_full_umap"] = input.nadh_vs_full_umap()
        previous_states_umap["fad_vs_full_umap"] = input.fad_vs_full_umap()

    @render.plot
    def plotPCA(): 
        experiment = input["select_experiment_pca"]()
        filtered_by = input["select_filter_pca"]()
        filter = input["filter_pca"]()
        color_by = input["color_pca"]()
        if filter is None or experiment is None or color_by is None or filtered_by is None:
            return None
        show_loadings = input["loadings"]()
        show_nadh = input["nadh_vs_full"]()
        show_fadh = input["fad_vs_full"]()
        num_plots = 1
        if show_loadings or show_nadh or show_fadh:
            num_plots = 2
        fig, axes = plt.subplots(1, num_plots, figsize=(10/num_plots, 5))
        if filtered_by == "Nothing":
            subsetDF = df[df["experiment"] == experiment]
        else: 
            subsetDF = df[(df["experiment"] == experiment) & (df[filtered_by] == filter)]
        if subsetDF.empty:
            return None
        sd_df = StandardScaler().fit_transform(subsetDF[measures])
        pca = PCA(n_components=2)
        principalComponents = pca.fit_transform(sd_df)
        principalDF = pd.DataFrame(data=principalComponents, columns=['PC1', 'PC2'])

       
        if num_plots == 2:
            if show_loadings: 
                loadings = pca.components_.T  # Transposing so each column corresponds to a PC
                loadings_df = pd.DataFrame(loadings, columns=['PC1', 'PC2'], index=measures)
                sns.heatmap(loadings_df, annot=True, cmap='coolwarm', center=0)
                axes[1].set_title('Contribution of Original Variables to PCs', fontsize=10, fontweight='bold')
            else:
                partial_measure = nadh_measures if show_nadh else fadh_measures
                partial_sd_df = StandardScaler().fit_transform(subsetDF[partial_measure])
                partial_pca = PCA(n_components=2)
                partial_principalComponents = partial_pca.fit_transform(partial_sd_df)
                partial_principalDF = pd.DataFrame(data=partial_principalComponents, columns=['PC1', 'PC2'])
                subsetDF = subsetDF.reset_index(drop=True)
                partial_principalDF  = pd.concat([partial_principalDF, subsetDF[color_by]], axis=1)
                treatments = partial_principalDF[color_by].unique()
                for t in treatments:
                    t_df = partial_principalDF[partial_principalDF[color_by] == t]
                    axes[1].scatter(t_df["PC1"],t_df["PC2"], label=t)
                axes[1].legend()
            
                axes[1].set_xlabel(f'Principal Component 1: {round(partial_pca.explained_variance_ratio_[0]*100, 2)}%')
                axes[1].set_ylabel(f'Principal Component 2: {round(partial_pca.explained_variance_ratio_[1]*100, 2)}%')

            ax = axes[0]
        else:
            ax = axes

        if not show_nadh or not show_fadh:
            subsetDF = subsetDF.reset_index(drop=True)
        principalDF = pd.concat([principalDF, subsetDF[color_by]], axis=1)
        treatments = principalDF[color_by].unique()
        for t in treatments:
            t_df = principalDF[principalDF[color_by] == t]
            ax.scatter(t_df["PC1"],t_df["PC2"], label=t)
        ax.legend()
    
        ax.set_xlabel(f'Principal Component 1: {round(pca.explained_variance_ratio_[0]*100, 2)}%')
        ax.set_ylabel(f'Principal Component 2: {round(pca.explained_variance_ratio_[1]*100, 2)}%')
        return fig
    
    @render.plot
    def plotUMAP(): 
        experiment = input["select_experiment_umap"]()
        filtered_by = input["select_filter_umap"]()
        filter = input["filter_umap"]()
        color_by = input["color_umap"]()
        if filter is None or experiment is None or color_by is None or filtered_by is None:
            return None
        
        n_neighbors = input["n_neighbors"]()
        min_dist = input["min_dist"]()
        show_nadh = input["nadh_vs_full_umap"]()
        show_fadh = input["fad_vs_full_umap"]()
        num_plots = 1
        if show_nadh or show_fadh:
            num_plots = 2
        fig, axes = plt.subplots(1, num_plots, figsize=(10/num_plots, 5))

        if filtered_by == "Nothing":
            subsetDF = df[df["experiment"] == experiment]
        else: 
            subsetDF = df[(df["experiment"] == experiment) & (df[filtered_by] == filter)]
        if subsetDF.empty:
            return None
        sd_df = StandardScaler().fit_transform(subsetDF[measures])
        reducer = umap.UMAP(n_components=2, random_state=42, n_neighbors=n_neighbors,
        min_dist=min_dist)
        umapDF = pd.DataFrame(data=reducer.fit_transform(sd_df), columns=['umap1', 'umap2'])

        if num_plots == 2:
            partial_measure = nadh_measures if show_nadh else fadh_measures
            partial_sd_df = StandardScaler().fit_transform(subsetDF[partial_measure])
            partial_umap = umap.UMAP(n_components=2, random_state=42, n_neighbors=n_neighbors, min_dist=min_dist)
            partial_umapDF = pd.DataFrame(data=partial_umap.fit_transform(partial_sd_df), columns=['umap1', 'umap2'])
            subsetDF = subsetDF.reset_index(drop=True)
            partial_umapDF  = pd.concat([partial_umapDF, subsetDF[color_by]], axis=1)
            treatments = partial_umapDF[color_by].unique()
            for t in treatments:
                t_df = partial_umapDF[partial_umapDF[color_by] == t]
                axes[1].scatter(t_df['umap1'],t_df["umap2"], label=t)
            axes[1].legend()
        
            axes[1].set_xlabel(f'UMAP Component 1')
            axes[1].set_ylabel(f'UMAP Component 2')
            ax = axes[0]
        else:
            ax = axes

        if not show_nadh or not show_fadh:
            subsetDF = subsetDF.reset_index(drop=True)
        umapDF = pd.concat([umapDF, subsetDF[color_by]], axis=1)
        treatments = umapDF[color_by].unique()
        for t in treatments:
            t_df = umapDF[umapDF[color_by] == t]
            ax.scatter(t_df["umap1"],t_df["umap2"], label=t)
        ax.legend()
    
        ax.set_xlabel(f'UMAP Component 1')
        ax.set_ylabel(f'UMAP Component 2')
        return fig
    
    @render.plot
    def plotTSNE(): 
        experiment = input["select_experiment_tsne"]()
        filtered_by = input["select_filter_tsne"]()
        filter = input["filter_tsne"]()
        color_by = input["color_tsne"]()
        if filter is None or experiment is None or color_by is None or filtered_by is None:
            return None
        perplexity = input["perplexity"]()
        fig, ax = plt.subplots(1, 1, figsize=(10, 5))

        subsetDF = df[(df["experiment"] == experiment) & (df[filtered_by] == filter)]
        if subsetDF.empty:
            return None
        sd_df = StandardScaler().fit_transform(subsetDF[measures])
        tsne = TSNE(n_components=2, perplexity=perplexity)
        tsneDF = pd.DataFrame(data=tsne.fit_transform(sd_df), columns=['t-SNE1', 't-SNE2'])
        subsetDF = subsetDF.reset_index(drop=True)
        tsneDF = pd.concat([tsneDF, subsetDF[color_by]], axis=1)
        treatments = tsneDF[color_by].unique()
       

        for t in treatments:
            t_df = tsneDF[tsneDF[color_by] == t]
            ax.scatter(t_df["t-SNE1"],t_df["t-SNE2"], label=t)
        ax.legend()
    
        ax.set_xlabel(f't-SNE Component 1')
        ax.set_ylabel(f't-SNE Component 2')
        return fig
    
    # dynmaically update the filter options
    @reactive.effect
    def _():
        filter_by = input["select_filter_pca"]()
        experiment = input["select_experiment_pca"]()

        if filter_by == "cell_line":
            ui.update_selectize("filter_pca", choices=experiments_cellline_dict[experiment])
        elif filter_by == "treatment":
            ui.update_selectize("filter_pca", choices=experiments_treatment_dict[experiment])
        elif filter_by == "cancer":
            ui.update_selectize("filter_pca", choices=experiments_cancer_dict[experiment])
        elif filter_by == "cell_type":
            ui.update_selectize("filter_pca", choices=experiments_celltype_dict[experiment])
        elif filter_by == "media":
            ui.update_selectize("filter_pca", choices=experiments_media_dict[experiment])
        else:
            ui.update_selectize("filter_pca", choices=[])
        if filter_by != "Nothing":
            color_by = pca_filters.copy()
            color_by.remove(filter_by)
            ui.update_selectize("color_pca", choices=color_by)
        else: 
            ui.update_selectize("color_pca", choices=pca_filters)

    
    @reactive.effect  
    def _():
        filter_by = input["select_filter_umap"]()
        experiment = input["select_experiment_umap"]()

        if filter_by == "cell_line":
            ui.update_selectize("filter_umap", choices=experiments_cellline_dict[experiment])
        elif filter_by == "treatment":
            ui.update_selectize("filter_umap", choices=experiments_treatment_dict[experiment])
        elif filter_by == "cancer":
            ui.update_selectize("filter_umap", choices=experiments_cancer_dict[experiment])
        elif filter_by == "cell_type":
            ui.update_selectize("filter_umap", choices=experiments_celltype_dict[experiment])
        elif filter_by == "media":
            ui.update_selectize("filter_umap", choices=experiments_media_dict[experiment])
        else:
            ui.update_selectize("filter_umap", choices=[])
        if filter_by != "Nothing":
            color_by = pca_filters.copy()
            color_by.remove(filter_by)
            ui.update_selectize("color_umap", choices=color_by)
        else: 
            ui.update_selectize("color_umap", choices=pca_filters)
    # dynmaically update the filter options
    @reactive.effect
    def _():
        filter_by = input["select_filter_tsne"]()
        experiment = input["select_experiment_tsne"]()

        if filter_by == "cell_line":
            ui.update_selectize("filter_tsne", choices=experiments_cellline_dict[experiment])
        elif filter_by == "treatment":
            ui.update_selectize("filter_tsne", choices=experiments_treatment_dict[experiment])
        elif filter_by == "cancer":
            ui.update_selectize("filter_tsne", choices=experiments_cancer_dict[experiment])
        elif filter_by == "cell_type":
            ui.update_selectize("filter_tsne", choices=experiments_celltype_dict[experiment])
        elif filter_by == "media":
            ui.update_selectize("filter_tsne", choices=experiments_media_dict[experiment])
        else:
            ui.update_selectize("filter_tsne", choices=[])
        color_by = pca_filters.copy()
        color_by.remove(filter_by)
        ui.update_selectize("color_tsne", choices=color_by)

    @render_widget
    def volcano_plot():
        experiment = input["select_experiment_v"]()
        measure = input["select_var_v"]()
        log2fc_cutoff = input["log2fc_cutoff"]()
        plot_title = "Volcano Plot for " + measure + " in experiment " + experiment
        x_axis_title = "log2 fold change" 
        y_axis_title = "-log10 pvalue" 
        point_radius = 6
        fig = go.Figure()
        fig.update_layout(
            title=plot_title,
            xaxis_title= x_axis_title,
            yaxis_title=y_axis_title,
            paper_bgcolor= 'white',
            plot_bgcolor='white',
        )
        colors = []
        data = volcano_df[(volcano_df['experiment'] == experiment)]
        log2fc = data[measure + '_log2fc'].tolist()
        nlog10pv = data[measure + '_pvadj'].apply(lambda x: -np.log10(x)).tolist()
        # assign color to each point
        for i in range(0, len(log2fc)):
            if nlog10pv[i] > 2:
                if log2fc[i] > log2fc_cutoff:
                    colors.append('#db3232')
                elif log2fc[i] < log2fc_cutoff * -1:
                    colors.append('#3f65d4')
                else:
                    colors.append('rgba(150,150,150,0.5)')
            else:
                colors.append('rgba(150,150,150,0.5)')
        fig.add_trace(
            go.Scattergl(
                x = log2fc,
                y = nlog10pv,
                mode = 'markers',
                hovertemplate = data['cell_line'] + ', ' + data['treatment'] + '<extra></extra>',
                marker= {
                    'color':colors,
                    'size':point_radius,
                }
            )
        )
        return fig
    
    # Reactive expression to filter the DataFrame based on the inputs
    @render.data_frame
    def filtered_table():
        filtered_df = grouped_df
        for col in column_names[:7]:
            if "All" not in input[col]():
                filtered_df = filtered_df[filtered_df[col].isin(input[col]())]
        return filtered_df
        
    second_var_switch = reactive.value(False)
    second_cell_line_switch = reactive.value(False)

    @reactive.Effect
    @reactive.event(input["2nd_cellline_switch"])
    def _():
        second_cell_line_switch.set(input["2nd_cellline_switch"]())

    @reactive.Effect
    @reactive.event(input["2nd_measure_switch"])
    def _():
        second_var_switch.set(input["2nd_measure_switch"]())

    @render.ui
    def turn_on_2nd():
        selected_cell_line = input["select_cell_line"]()
        selected_measure = input["select_var"]()
        one_cell_line = (len(selected_cell_line) == 1 and "All cell lines" not in selected_cell_line)
        one_measure = (len(selected_measure) == 1 and "All measures" not in selected_measure)
        all_cell_lines = (len(selected_cell_line) == 1 and "All cell lines" in selected_cell_line)
        all_measures = (len(selected_measure) == 1 and "All measures" in selected_measure)
        more_than_one_cell_line = (len(selected_cell_line) > 1)
        more_than_one_measures = (len(selected_measure) > 1)
        if one_cell_line and (more_than_one_measures or all_measures):
            return ui.input_switch("2nd_cellline_switch", "Do you want to plot another cell line?", False)
        elif one_measure and (more_than_one_cell_line or all_cell_lines):
            return ui.input_switch("2nd_measure_switch", "Do you want to plot another measure?", False)
        else: 
            second_var.set("None")
            second_cell_line.set("None")
            second_var_switch.set(False)
            second_cell_line_switch.set(False)
            
        return None   
     
    @render.ui
    def add_2nd_cellline():      
        if input["2nd_cellline_switch"]() and second_cell_line_switch():           
            selected_cell_line = input["select_cell_line"]()
            if len(selected_cell_line) < 1:
                return None
            choices = unique_values_dict[column_names[0]].copy()
            choices.remove(selected_cell_line[0])
            return ui.card(ui.input_selectize("select_2nd_cell_line", "Select a 2nd cell line", choices, multiple=False), height="200px")      
        else:
            return None
        
    @render.ui
    def add_2nd_measure():
        if input["2nd_measure_switch"]() and second_var_switch():
            selected_measure = input["select_var"]()
            if len(selected_measure) < 1:
                return None
            choices = vars.copy()
            choices.remove(selected_measure[0])
            return ui.card(ui.input_selectize("select_2nd_var", "Select a 2nd measure", choices, multiple=False), height="200px")
        else:
            return None
    
    @reactive.Effect   
    def _():
        for col in column_names[:7]:
            if col == "treatment":
                continue
            i = input[col]()
            if "All" in i:
                ui.update_checkbox_group(col, choices=["All"], selected=["All"])
            else:
                ui.update_checkbox_group(col, choices=["All"] + unique_values_dict[col], selected=i)
    
    @reactive.Effect
    def _():
        def filterTreatmentbyExperiment(exp):
            treatments = []
            if "All" in exp:
                return ["All"]
            else:
                for i in exp:
                    treatments += experiments_treatment_dict[i]
            return treatments

        x = input[column_names[3]]()
        # Can also set the label and select items
        ui.update_checkbox_group(
            column_names[4],
            label="Treatment",
            choices=filterTreatmentbyExperiment(x),
            selected=filterTreatmentbyExperiment(x),
        )

    @reactive.Effect
    def _():
        selected_cell_line = input["select_cell_line"]()
        selected_measure = input["select_var"]()
        cell_line_empty = (selected_cell_line is None or len(selected_cell_line) == 0)
        measure_empty = (selected_measure is None or len(selected_measure) == 0)
        one_cell_line = (len(selected_cell_line) == 1 and "All cell lines" not in selected_cell_line)
        one_measure = (len(selected_measure) == 1 and "All measures" not in selected_measure)
        all_cell_lines = (len(selected_cell_line) == 1 and "All cell lines" in selected_cell_line)
        all_measures = (len(selected_measure) == 1 and "All measures" in selected_measure)
        more_than_one_cell_line = (len(selected_cell_line) > 1)
        more_than_one_measures = (len(selected_measure) > 1)
        
        if cell_line_empty and measure_empty:
            # print("all empty")
            ui.update_selectize(id="select_cell_line", label="Select cell line(s)", choices=["All cell lines"] + unique_values_dict[column_names[0]])
            ui.update_selectize(id="select_var", label="Select measure(s)", choices=["All measures"] + vars)
        elif cell_line_empty:
            if one_measure:
                # print("cell line empty, one measure")
                ui.update_selectize(id="select_cell_line", label="Select cell line(s)", choices=["All cell lines"] + unique_values_dict[column_names[0]])
                ui.update_selectize(id="select_var", label="Select measure(s)", choices=vars, selected=selected_measure)
            elif all_measures:
                # print("cell line empty, all measures")
                ui.update_selectize(id="select_var", label="You have chosen all measures", choices=selected_measure, selected=selected_measure)
                ui.update_selectize(id="select_cell_line", label="Select cell line(s)", choices=unique_values_dict[column_names[0]])                
            elif more_than_one_measures:
              #  print("cell line empty, more than one measures")            
                ui.update_selectize(id="select_cell_line", label="Select cell line(s).", choices=unique_values_dict[column_names[0]])
        elif measure_empty:
            if one_cell_line:
                # print("measure empty, one cell line")
                ui.update_selectize(id="select_var", label="Select measure(s)", choices=["All measures"] + vars)
                ui.update_selectize(id="select_cell_line", label="Select cell line(s)", choices=unique_values_dict[column_names[0]], selected=selected_cell_line)
            elif all_cell_lines:
                # print("measure empty, all cell lines")
                ui.update_selectize(id="select_cell_line", label="You have chosen all cell lines", choices=selected_cell_line, selected=selected_cell_line)
                ui.update_selectize(id="select_var", label="Select measure(s)", choices=vars)
            elif more_than_one_cell_line:
                # print("measure empty, more than one cell line")
                ui.update_selectize(id="select_var", label="Select measure(s)", choices=vars)
        elif one_cell_line and one_measure:
            # print("one cell line, one measure")
            ui.update_selectize(id="select_cell_line", label="Select cell line(s)", choices=unique_values_dict[column_names[0]], selected=selected_cell_line)
            ui.update_selectize(id="select_var", label="Select measure(s)", choices=vars, selected=selected_measure)
        elif one_cell_line and (all_measures or more_than_one_measures):
            # print("one cell line, more than one measures")
            ui.update_selectize(id="select_cell_line", label="If you want > 1 cell lines, please select only 1 measure or add another cell line below.", choices=selected_cell_line, selected=selected_cell_line)
            if all_measures:
                ui.update_selectize(id="select_var", label="You have chosen all measures.", choices=selected_measure, selected=selected_measure)
        elif one_measure and (all_cell_lines or more_than_one_cell_line):
            # print("one measure, more than one cell lines")
            ui.update_selectize(id="select_var", label="If you want > 1 measures, please select only 1 cell line or add another measure below.", choices=selected_measure, selected=selected_measure)
            if all_cell_lines:
                ui.update_selectize(id="select_cell_line", label="You have chosen all cell lines.", choices=selected_cell_line, selected=selected_cell_line)
        else: 
            print("illegal combination")

    second_var = reactive.value("None")
    second_cell_line = reactive.value("None")
   
    @reactive.Effect
    @reactive.event(input["2nd_measure_switch"], input["select_2nd_var"])
    def _():
        second_var.set(input["select_2nd_var"] if input["2nd_measure_switch"]() else "None")
    
    @reactive.Effect
    @reactive.event(input["2nd_cellline_switch"], input["select_2nd_cell_line"])
    def _():
        second_cell_line.set(input["select_2nd_cell_line"] if input["2nd_cellline_switch"]() else "None")

    @render.plot
    def plotTreatmentsXcontrol():
        experiment = input["select_experiment"]()
        cell_lines = input["select_cell_line"]()
        measures = input["select_var"]() 
        var2 = second_var()
        cell_line2 = second_cell_line()
        
        if experiment:
            treatments = experiments_treatment_dict[experiment]
            treatments = [t for t in treatments if "control" not in t]
        if "All cell lines" in cell_lines:
            cell_lines = unique_values_dict[column_names[0]]       
        if "All measures" in measures:
            measures = vars
        if (len(measures) > 1 and len(cell_lines) > 1) or len(measures) < 1 or len(cell_lines) < 1:
            return None

        num_plots = 1 if (var2 == "None" and cell_line2 == "None") else 2
        fig, axes = plt.subplots(1, num_plots, figsize=(10/num_plots, 5))

        for i in range(0, num_plots):
            ax = axes if num_plots == 1 else axes[i] 
           
            if i == 1: # plot the 2nd graph:
                g2 = var2() if var2 != "None" else cell_line2()
                print(g2)

            if len(measures) > 1:
                ys = measures
                isY = "measures"
            else:
                ys = cell_lines
                isY = "cell_lines"
            y_dict = {}
            for y in ys:
                y_dict[y] = ([], []) # the first list is the measure, the second list is the error (se)
                for treatment in treatments:
                    if isY == "measures":
                        cellline = cell_lines[0] if i == 0 else g2
                        row = grouped_df_norm[(grouped_df_norm["experiment"] == experiment) & (grouped_df_norm["treatment"] == treatment) & (grouped_df_norm["cell_line"] == cellline)]
                        if len(row) > 1 or len(row) == 0:
                            value = np.nan
                            error = np.nan
                        else:
                            value = float(row[vars_colnames[vars.index(y)]].iloc[0])
                            error = float(row[vars_colnames[vars.index(y)].replace("_mean", "_se")].iloc[0])
                        y_dict[y][0].append(value)
                        y_dict[y][1].append(error)
                    else:
                        row = grouped_df_norm[(grouped_df_norm["experiment"] == experiment) & (grouped_df_norm["treatment"] == treatment) & (grouped_df_norm["cell_line"] == y)]
                        if len(row) > 1 or len(row) == 0:
                            value = np.nan
                            error = np.nan
                        else:
                            measure = measures[0] if i == 0 else g2
                            value = float(row[vars_colnames[vars.index(measure)]].iloc[0])
                            error = float(row[vars_colnames[vars.index(measure)].replace("_mean", "_se")].iloc[0])
                        y_dict[y][0].append(value)
                        y_dict[y][1].append(error)
            
            # plot the data
            for y_name, (values, errors) in y_dict.items():
                # Plot each measure with error bars
                if isY == "measures":
                    ax.errorbar(treatments, values, yerr=errors, fmt='o', label=y_name, color=var_colors[y_name])
                else:
                    ax.errorbar(treatments, values, yerr=errors, fmt='o', label=y_name, color=cell_line_colors[y_name])
                # ax.errorbar(treatments, values, yerr=errors, fmt='o', label=y_name)
            
            ax.axhline(y=100, color='gray', linestyle='-', label='control')
        
            if i == 0:
                ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=10)
           
            # Customization
            ax.set_xlabel('Treatments', fontsize=14, fontweight='bold')
            ax.tick_params(axis='x', labelsize=8)  
            ax.tick_params(axis='y', labelsize=10)
            if (i == 0):
                if isY == "measures":
                    ax.set_title(f"{cell_lines[0]}", fontsize=14, fontweight='bold')
                else:
                    ax.set_title(f"{measures[0]}", fontsize=14, fontweight='bold')
            else:
                ax.set_title(f"{g2}", fontsize=14, fontweight='bold')

         # Set common y-axis label
        fig.text(0, 0.5, '% change VS control', va='center', rotation='vertical', fontsize=12, fontweight='bold')
       
        fig.subplots_adjust(right=0.9, wspace=0.1)
        return fig

# Create the Shiny app
app = App(app_ui, server)

if __name__ == "__main__":
    app.run()