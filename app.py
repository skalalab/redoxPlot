from shiny import App, ui, render, reactive
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from htmltools import HTML
from pathlib import Path

grouped_df = pd.read_csv(Path(__file__).parent / "grouped_df.csv")
grouped_df = grouped_df.rename(columns={grouped_df.columns[0]: 'id'})
column_names = grouped_df.columns.tolist()[1:8]
unique_values_dict = {col: grouped_df[col].unique().tolist() for col in column_names}
experiments_treatment_dict = {name: group["treatment"].tolist() for name, group in grouped_df.groupby("experiment")}
vars_colnames = ["na1_mean", "na2_mean","nt1_mean","nt2_mean","ntm_mean", "fa1_mean","fa2_mean","ft1_mean","ft2_mean","ftm_mean", "nint_mean", "fint_mean", "normrr_mean"]
vars = [name.replace("_mean", "") for name in vars_colnames[:-1]] + ["rr"]

# generate colors for vars and cell lines
cm = plt.get_cmap('tab20')
var_colors = {var: cm(i) for i, var in enumerate(vars)}
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
                ui.panel_sidebar(   
                    ui.layout_columns(
                        ui.card(ui.input_selectize("select_experiment", "Select experiment", unique_values_dict[column_names[3]]), height="250px"),
                        ui.card(ui.input_selectize("select_cell_line", "Select cell line(s)", ["All cell lines"] + unique_values_dict[column_names[0]], multiple=True), height="180px"),
                        ui.card(ui.input_selectize("select_var", "Select measure(s)", ["All measures"] + vars, multiple=True), height="250px"),
                        col_widths=12
                    ),
                ), ui.panel_main(
                    HTML("<h2>Plot percent change for treatments against control</h2>"),
                    ui.output_plot("plotTreatmentsXcontrol"),

                ))),          
            ui.nav_panel("Volcano Plot?", "To be continued..."),
            id="tab",)
    ),     
    ui.nav_panel("2", "Page B content"),  
    ui.nav_panel("3", "Page C content"),  
    title="Redox Plots",  
    id="page",  
    fillable = True, 
)  


# Define the server logic
def server(input, output, session):
    # Reactive expression to filter the DataFrame based on the inputs
    @render.data_frame
    def filtered_table():
        filtered_df = grouped_df
        for col in column_names[:7]:
            if "All" not in input[col]():
                filtered_df = filtered_df[filtered_df[col].isin(input[col]())]
        return filtered_df
    
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
                # print("cell line empty, more than one measures")
                ui.update_selectize(id="select_cell_line", label="You can only choose 1 cell line now. If you want other cell lines, please select only 1 measure.", choices=unique_values_dict[column_names[0]])
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
                ui.update_selectize(id="select_var", label="You can only choose 1 measure now. If you want other measures, please select only 1 cell line.", choices=vars)
        elif one_cell_line and one_measure:
            # print("one cell line, one measure")
            ui.update_selectize(id="select_cell_line", label="Select cell line(s)", choices=unique_values_dict[column_names[0]], selected=selected_cell_line)
            ui.update_selectize(id="select_var", label="Select measure(s)", choices=vars, selected=selected_measure)
        elif one_cell_line and (all_measures or more_than_one_measures):
            # print("one cell line, more than one measures")
            ui.update_selectize(id="select_cell_line", label="You have chosen 1 cell line. If you want other cell lines, please select only 1 measure.", choices=selected_cell_line, selected=selected_cell_line)
            if all_measures:
                ui.update_selectize(id="select_var", label="You have chosen all measures.", choices=selected_measure, selected=selected_measure)
        elif one_measure and (all_cell_lines or more_than_one_cell_line):
            # print("one measure, more than one cell lines")
            ui.update_selectize(id="select_var", label="You have chosen 1 measure. If you want other measures, please select only 1 cell line.", choices=selected_measure, selected=selected_measure)
            if all_cell_lines:
                ui.update_selectize(id="select_cell_line", label="You have chosen all cell lines.", choices=selected_cell_line, selected=selected_cell_line)
        else: 
            print("illegal combination")

    @render.plot
    def plotTreatmentsXcontrol():
        experiment = input["select_experiment"]()
        cell_lines = input["select_cell_line"]()
        measures = input["select_var"]()    
        if experiment:
            treatments = experiments_treatment_dict[experiment]
            treatments = [t for t in treatments if "control" not in t]
        if "All cell lines" in cell_lines:
            cell_lines = unique_values_dict[column_names[0]]       
        if "All measures" in measures:
            measures = vars
        if (len(measures) > 1 and len(cell_lines) > 1) or len(measures) < 1 or len(cell_lines) < 1:
            return None
              
        fig, ax = plt.subplots(1,1)
        ax.set_title("% Change vs Control in " + r"$\bf{" +experiment+ "}$" + " across treatments", fontsize=15)
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
                    row = grouped_df[(grouped_df["experiment"] == experiment) & (grouped_df["treatment"] == treatment) & (grouped_df["cell_line"] == cell_lines[0])]
                    if len(row) > 1 or len(row) == 0:
                        value = np.nan
                        error = np.nan
                    else:
                        value = float(row[vars_colnames[vars.index(y)]].iloc[0])
                        error = float(row[vars_colnames[vars.index(y)].replace("_mean", "_se")].iloc[0])
                    y_dict[y][0].append(value)
                    y_dict[y][1].append(error)
                else:
                    row = grouped_df[(grouped_df["experiment"] == experiment) & (grouped_df["treatment"] == treatment) & (grouped_df["cell_line"] == y)]
                    if len(row) > 1 or len(row) == 0:
                        value = np.nan
                        error = np.nan
                    else:
                        value = float(row[vars_colnames[vars.index(measures[0])]].iloc[0])
                        error = float(row[vars_colnames[vars.index(measures[0])].replace("_mean", "_se")].iloc[0])
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
        ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=12)
        plt.subplots_adjust(right=0.8)
        # Customization
        ax.set_xlabel('Treatments', fontsize=14, fontweight='bold')
        ax.set_ylabel('% change VS control', fontsize=14, fontweight='bold')
        ax.tick_params(axis='x', labelsize=10)  
        ax.tick_params(axis='y', labelsize=12)  

        return fig

# Create the Shiny app
app = App(app_ui, server)

if __name__ == "__main__":
    app.run()