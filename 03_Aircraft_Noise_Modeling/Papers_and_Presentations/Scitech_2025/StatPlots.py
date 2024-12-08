import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np 
import os
import matplotlib.cm as cm
from RCAIDE.Framework.Core import  Data
import scipy.stats as stats

from matplotlib.patches import Rectangle
from matplotlib.collections import PolyCollection
from matplotlib.colors import to_rgb
from matplotlib.legend_handler import HandlerTuple

from scipy.stats import f_oneway
from statsmodels.stats.multicomp import pairwise_tukeyhsd

from scipy.stats import f_oneway, kruskal, shapiro, levene

import scikit_posthocs as sp
from sklearn.datasets import load_iris






race_dictionary = {
    "Total": "B03002001",
    "White": "B03002003",
    "Black": "B03002004",
    "American Indian": "B03002005",
    "Asian": "B03002006",
    "Pacific Islander": "B03002007",
    "Hispanic": "B03002012"
}

income_dictionary = {
    "Total": "B19001001",
    "<$10k": "B19001002",
    "$10k to $15k": "B19001003",
    "$15k to $20k": "B19001004",
    "$20k to $25k": "B19001005",
    "$25k to $30k": "B19001006",
    "$30k to $35k": "B19001007",
    "$35k to $40k": "B19001008",
    "$40k to $45k": "B19001009",
    "$45k to 50k": "B19001010",
    "$50k to $60k": "B19001011",
    "$60k to $75k": "B19001012",
    "$75k to $100k": "B19001013",
    "$100k to $125k": "B19001014",
    "$125k to $150": "B19001015",
    "$150k to $200k": "B19001016",
    "$200k+ ": "B19001017"
}


income_aggregation = {
    '<50K': [
        "<$10k",
        "$10k to $15k",
        "$15k to $20k",
        "$20k to $25k",
        "$25k to $30k",
        "$30k to $35k",
        "$35k to $40k",
        "$40k to $45k",
        "$45k to $50k",
    ],
    '50k to 100k': [
        "$50k to $60k",
        "$60k to $75k",
        "$75k to $100k",
    ],
    '100k to 150k': [
        "$100k to $125k",
        "$125k to $150",
    ],
    '150k to 200k': [
        "$150k to $200k",
    ],
    '200k+': [
        "$200k+ ",
    ],
}





# data_file_tract = pd.read_csv('HC_Data_Tract_CA_ONLY.csv')
# data_file_tract = pd.read_csv('HC_Data_Tract_LogPOP_Only.csv')

data_file_tract = pd.read_csv('TR_Data_Tract.csv')


# Define the relevant columns
race_columns = ['Asian','Black','Hispanic','Native American','Pacific Islander', 'White']

income_columns =['<50K','50k to 100k','100k to 150k','150k to 200k','200k+']
income_raw  = list(income_dictionary.keys())[1:]
struct_columns = ['Churches','Hospitals/Medical Centers','Schools']

for income in income_columns:
    data_file_tract[income]= data_file_tract[income_aggregation[income]].sum(axis=1)


noise_column = 'L_dn'

def stackedbplot_race(data_frame, column_names, noise_type, save_name,col_names):
    folder_path = 'Bar_Plots'

    # Check if the folder exists, and if not, create it
    if not os.path.isdir(folder_path):
        os.mkdir(folder_path)

    noise_threshold_1 = 45
    noise_threshold_2 = 60

    data_frame_thresh45 = data_frame[data_frame[noise_type] >= noise_threshold_1]
    data_frame_thresh60 = data_frame[data_frame[noise_type] >= noise_threshold_2]

    values45 = data_frame_thresh45[column_names].sum() / data_frame[column_names].sum() * 100
    values2_45 = abs(data_frame[column_names].sum() - data_frame_thresh45[column_names].sum()) / data_frame[column_names].sum() * 100

    values60 = data_frame_thresh60[column_names].sum() / data_frame[column_names].sum() * 100
    values2_60 = abs(data_frame[column_names].sum() - data_frame_thresh60[column_names].sum()) / data_frame[column_names].sum() * 100

    # X positions for bars
    x = np.arange(len(column_names))  # number of columns
    bar_width = .25 # Width of each bar

    # pp = define_plot_parameters()

    base_colors = cm.tab20(np.linspace(0, 1, 20))

    plt.figure(figsize=(4,4))
    plt.rcParams["font.family"] = "Times New Roman"
    parameters = {'axes.labelsize': 11,
        'legend.fontsize': 9,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10}
    plt.rcParams.update(parameters)


    # Plot bars for each column with custom colors
    for idx, column in enumerate(column_names):
        base_color = base_colors[idx*2]
        lighter_color = [min(1, c + 0.3) for c in base_color]  

        # Plot bars for 45 dB threshold
        plt.bar(x[idx], values45[column], width=bar_width, color=base_color, label=None if idx > 0 else r"> 45 dBA $L_{dn}$")
        plt.bar(x[idx], values2_45[column], width=bar_width, bottom=values45[column], color=lighter_color, label=None if idx > 0 else   r"< 45 dBA $L_{dn}$")

        # Plot bars for 60 dB threshold with hatching
        plt.bar(x[idx] + bar_width, values60[column], width=bar_width, color=base_color, hatch='//', label=None if idx > 0 else  r"> 60 dBA $L_{dn}$")
        plt.bar(x[idx] + bar_width, values2_60[column], width=bar_width, bottom=values60[column], color=lighter_color, hatch='//', label=None if idx > 0 else  r" <60 dBA $L_{dn}$")

    # Add labels, title, and legend
    plt.xticks(x+bar_width/2, col_names)  # Set x-tick labels and positions
    plt.ylabel('Population %')
    plt.ylim(0,120)
    plt.legend(loc='upper center',ncol=2)

    # Save plot
    plot_file_name = os.path.join(folder_path, f'TR_{save_name}{noise_threshold_1}{noise_threshold_2}dB_Plot.png')
    plt.tight_layout()
    plt.savefig(plot_file_name, dpi=1100)
    plt.close()



def split_violin_plot(data_tract, data_columns, noise, x_axis,y_axis,col_names,y_ax):
    folder_path = 'Violin_Plots_V2'

    # Check if the folder exists, and if not, create it
    if not os.path.isdir(folder_path):
        os.mkdir(folder_path)

    noise_thresholds = [45,60]
    data=[]
    for threshold in noise_thresholds:
        filtered_tract = data_tract[data_tract[noise] >= threshold]

        val =filtered_tract[data_columns]
        val = val.reindex(sorted(val.columns), axis=1)
        val = val.set_axis(col_names, axis=1)
        # val['CA'] = filtered_tract['CA']
        val["L_dn"] = filtered_tract["L_dn"]
        val['threshold'] = threshold

        data.append(val)
    data = pd.concat(data)
    path_file = 'Data_sets'
    # Check if the folder exists, and if not, create it
    if not os.path.isdir(path_file):
        os.mkdir(path_file)
    file_path_name = os.path.join(path_file, f'{x_axis}{y_axis}.csv')
    data = data.apply(pd.to_numeric, errors='coerce') 
    data.to_csv(file_path_name,index=False)

    df = pd.DataFrame(data).apply(pd.to_numeric, errors='coerce') 
    

    # Reshape the data using melt()
    df_melted = df.melt(
        id_vars=["threshold",'L_dn'],
        value_vars=col_names, 
        var_name="Race", 
        value_name="vals"
    )
    
    df_melted = df_melted.dropna()


    base_colors = cm.tab20(np.linspace(0, 1, 20))
    for idx, column in enumerate(col_names):
        base_color = base_colors[idx*2]
        lighter_color = [min(1, c + 0.3) for c in base_color]  


    # Create a custom color palette
    plt.figure(figsize=(4,4))
    plt.rcParams["font.family"] = "Times New Roman"
    parameters = {'axes.labelsize': 11,
        'legend.fontsize': 9,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10}
    plt.rcParams.update(parameters)
    # sns.violinplot(data=df_melted,
    #     x="Race", y="vals", inner = 'quart', hue="threshold", gap=.5,
    #     split=True
    # )
    # # plt.xticks(col_names) 
    # plt.ylabel(f'{y_axis}')
    # plt.xlabel(f'{x_axis}')
    # # Save the plot

    # Generate the base color palette
    base_colors = cm.tab20(np.linspace(0, 1, 20))

    # Define the column names or categories you're working with
    col_names = df_melted['Race'].unique()

    # Create alternating colors (dark, light for each pair)
    custom_palette = []
    for idx, column in enumerate(col_names):
        base_color = base_colors[idx * 2]         # Darker color
        lighter_color = [min(1, c + 0.3) for c in base_color]  # Lighter version of the dark color
        custom_palette.extend([base_color, lighter_color])    # Append both colors

    # Set plot parameters
    plt.figure(figsize=(4, 4))
    plt.rcParams["font.family"] = "Times New Roman"
    parameters = {
        'axes.labelsize': 11,
        'legend.fontsize': 9,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10
    }
    plt.rcParams.update(parameters)

    # Create the violin plot using the custom palette
    ax = sns.violinplot(
        data=df_melted,
        x="Race", y="vals", inner='quart', hue="threshold", gap=.1,
        split=True, palette=custom_palette
    )

    plt.ylabel(f'{y_ax}')
    plt.xlabel(f'{x_axis}')

    # plt.ylabel(r'$L_{dn}$')
    # plt.xlabel(r'$L_{dn}$')

    handles = []
    for ind, violin in enumerate(ax.findobj(PolyCollection)):
        rgb = to_rgb(custom_palette[ind])
        if ind % 2 != 0:  # Lighten for the second split
            rgb = 0.5 + 0.5 * np.array(rgb)
        violin.set_facecolor(rgb)
        handles.append(Rectangle((0, 0), 0, 0, facecolor=rgb, edgecolor='black'))

    ax.legend(
        # handles=[tuple(handles[::2]), tuple(handles[1::2])],  # Group dark and light shades
        labels=[r"> 45 dBA $L_{dn}$",r"> 60 dBA $L_{dn}$"],
        handlelength=2,
        handler_map={tuple: HandlerTuple(ndivide=None, pad=0)},
        loc='upper left'
    )
    # ax.legend(handles=[tuple(handles)], labels=[r"> 45 dBA $L_{dn}$"], loc='upper left',handler_map={tuple: HandlerTuple(ndivide=None, pad=0)} )

    plot_file_name = os.path.join(folder_path, f'TRviolin_plot_{x_axis}{y_axis}.png')
    plt.tight_layout()
    plt.savefig(plot_file_name,dpi=800)


    # df_melted = df_melted.dropna()

    # groups = [group['vals'] for _, group in df_melted.groupby(['Race', 'threshold'])]
    # # f_statistic, p_value = f_oneway(*groups)
    # stat, p_value = kruskal(*groups)
    # # stat, p_value = f_oneway(*groups)



    # groups = []
    # group_labels = []

    # for (race, threshold), group in df_melted.groupby(['Race', 'threshold']):
    #     groups.append(group['vals'].values)
    #     group_labels.append(f"Race {race}, threshold {threshold}")

    # # Directory to save the plots
    # output_dir = "plots"
    # os.makedirs(output_dir, exist_ok=True)


    # df_melted['group'] = df_melted['Race'] + df_melted['threshold'].astype(str)

    # anova_results = {
    # 'KWF-statistic': [stat],
    # 'p-value': [p_value]
    # }

    # print(x_axis,anova_results)
    # # Convert to a DataFrame for easier export
    # anova_df = pd.DataFrame(anova_results)


    # df_melted['group'] = df_melted['Race'] + df_melted['threshold'].astype(str)
    # dunn_p_values = sp.posthoc_dunn(df_melted, val_col='vals', group_col='group', p_adjust='holm')

    # dunn_p_values_df = pd.DataFrame(dunn_p_values)

    # # Save to a CSV file
    # # dunn_p_values_df.to_csv(output_file, index=True)

    # anova_df = pd.DataFrame(anova_results)
    # # anova_df.to_csv(anova_file, index=False)

    # # Combine Dunn's and ANOVA results into one file (optional)
    # combined_file = os.path.join(path_file,f"HC_{x_axis}combined_stat_results.csv")

    # with open(combined_file, 'w') as f:
    #     f.write("KWF Results\n")
    #     anova_df.to_csv(f, index=False)
    #     f.write("\nDunn's Test Results\n")
    #     dunn_p_values_df.to_csv(f, index=True)

def stats_analysis(data_tract, data_columns, noise, x_axis,y_axis,col_names,ac_type):

    noise_thresholds = [45]
    data=[]
    for threshold in noise_thresholds:
        filtered_tract = data_tract[data_tract[noise] >= threshold]

        val =filtered_tract[data_columns]
        val = val.reindex(sorted(val.columns), axis=1)
        val = val.set_axis(col_names, axis=1)
        # val['CA'] = filtered_tract['CA']
        val["L_dn"] = filtered_tract["L_dn"]
        val['threshold'] = threshold

        data.append(val)
    data = pd.concat(data)
    path_file = 'Data_sets'
    # Check if the folder exists, and if not, create it
    if not os.path.isdir(path_file):
        os.mkdir(path_file)
    file_path_name = os.path.join(path_file, f'{x_axis}{y_axis}.csv')
    data = data.apply(pd.to_numeric, errors='coerce') 
    data.to_csv(file_path_name,index=False)

    df = pd.DataFrame(data).apply(pd.to_numeric, errors='coerce') 
    

    # Reshape the data using melt()
    df_melted = df.melt(
        id_vars=["threshold",'L_dn'],
        value_vars=col_names, 
        var_name="Race", 
        value_name="vals"
    )
    
    df_melted = df_melted.dropna()

    groups = [group['vals'] for _, group in df_melted.groupby(['Race', 'threshold'])]
    # f_statistic, p_value = f_oneway(*groups)
    stat, p_value = kruskal(*groups)
    # stat, p_value = f_oneway(*groups)



    groups = []
    group_labels = []

    for (race, threshold), group in df_melted.groupby(['Race', 'threshold']):
        groups.append(group['vals'].values)
        group_labels.append(f"Race {race}, threshold {threshold}")

    # Directory to save the plots
    output_dir = "plots"
    os.makedirs(output_dir, exist_ok=True)



    stat_results = {
    'KWF-statistic': [stat],
    'p-value': [p_value]
    }

    print(x_axis,stat_results)
    # Convert to a DataFrame for easier export
    anova_df = pd.DataFrame(stat_results)


    df_melted['group'] = df_melted['Race'] + df_melted['threshold'].astype(str)
    df_melted.to_csv(f"{ac_type}_{x_axis}df_melted.csv",index=False)
    dunn_p_values = sp.posthoc_dunn(df_melted, val_col='vals', group_col='group', p_adjust='holm')

    dunn_p_values_df = pd.DataFrame(dunn_p_values)


    # Combine Dunn's and ANOVA results into one file (optional)
    combined_file = os.path.join(path_file,f"{ac_type}_{x_axis}combined_stat_results.csv")

    with open(combined_file, 'w') as f:
        f.write("KWF Results\n")
        anova_df.to_csv(f, index=False)
        f.write("\nDunn's Test Results\n")
        dunn_p_values_df.to_csv(f, index=True)


#------------------------#
col_names_race = ['Asian','Black','Hisp.','Nat.\nAmer.','Pac.\nIsld.', 'White']
# col_names_race = ['Asian','Black','Hisp.','NatAmer','PacIsld', 'White']

col_names_inc = ['<50K','50k-100k','100k-150k','150k-200k','200k+']
col_names_struct = ['Churches','Medical\nCenters','Schools']

# stackedbplot_race(data_file_tract,race_columns,'L_dn','Race',col_names_race)
# stackedbplot_race(data_file_tract,income_columns,'L_dn','Income',col_names_inc)
# stackedbplot_race(data_file_tract,struct_columns,'L_dn','Structures',col_names_struct)



# race_columns = ['Asian_CA','Black_CA','Hispanic_CA','Native American_CA','Pacific Islander_CA', 'White_CA']
# income_columns =['<50K_CA','50k to 100k_CA','100k to 150k_CA','150k to 200k_CA','200k+_CA']

race_columns = ['Asian_LogPOP','Black_LogPOP','Hispanic_LogPOP','Native American_LogPOP','Pacific Islander_LogPOP', 'White_LogPOP']
income_columns =['<50K_LogPOP','50k to 100k_LogPOP','100k to 150k_LogPOP','150k to 200k_LogPOP','200k+_LogPOP']


# files= ['HC_Raw_Data_Tract.csv','TR_Raw_Data_Tract.csv']
files= ['TR_Raw_Data_Tract.csv']

# y_axis_list = ['CA']
# y_ax_label = ['CA']

y_axis_list = ['Log(PD)L_dn']
y_ax_label = [r'$L_{dn}\cdot log_{10}(PD)$']
ac_list= ['HC','TR']
i=0
for file in files:
    data_file_tract = pd.read_csv(file)
    # stats_analysis(data_file_tract,race_columns,'L_dn','Race',y_axis_list[0],col_names_race,ac_list[i])
    

    split_violin_plot(data_file_tract,race_columns,'L_dn','Race',y_axis_list[i],col_names_race,y_ax_label[i])
    split_violin_plot(data_file_tract,income_columns,'L_dn','Income ($)',y_axis_list[i],col_names_inc,y_ax_label[i])
    # split_violin_plot(data_file_tract,race_columns,'L_dn','L_dn',['L_dn'],col_names_race,'L_dn')

    i+=1

#show the entire table
#have a table statically singifcant
#True false - 

# % Please add the following required packages to your document preamble:
# % \usepackage{multirow}
# % \usepackage[table,xcdraw]{xcolor}
# % Beamer presentation requires \usepackage{colortbl} instead of \usepackage[table,xcdraw]{xcolor}
# \begin{table}[]
# \begin{tabular}{|c|c|cccc|cccc|}
# \hline
#  &  & \multicolumn{4}{c|}{\textbf{\textgreater{}45 dBA}} & \multicolumn{4}{c|}{\textbf{\textgreater{}60dBA}} \\ \cline{3-10} 
# \multirow{-2}{*}{\textbf{Aircraft}} & \multirow{-2}{*}{\textbf{Race}} & \multicolumn{1}{c|}{Asian} & \multicolumn{1}{c|}{Black} & \multicolumn{1}{c|}{Hispanic} & White & \multicolumn{1}{c|}{Asian} & \multicolumn{1}{c|}{Black} & \multicolumn{1}{c|}{Hispanic} & White \\ \hline
#  & Asian & \multicolumn{1}{c|}{} & \multicolumn{1}{c|}{} & \multicolumn{1}{c|}{} &  & \multicolumn{1}{c|}{} & \multicolumn{1}{c|}{} & \multicolumn{1}{c|}{} &  \\ \cline{2-10} 
#  & Black & \multicolumn{1}{c|}{} & \multicolumn{1}{c|}{} & \multicolumn{1}{c|}{} &  & \multicolumn{1}{c|}{} & \multicolumn{1}{c|}{} & \multicolumn{1}{c|}{} &  \\ \cline{2-10} 
#  & Hispanic & \multicolumn{1}{c|}{} & \multicolumn{1}{c|}{} & \multicolumn{1}{c|}{} &  & \multicolumn{1}{c|}{} & \multicolumn{1}{c|}{} & \multicolumn{1}{c|}{} &  \\ \cline{2-10} 
# \multirow{-4}{*}{Stop Tilt Rotor} & White & \multicolumn{1}{c|}{} & \multicolumn{1}{c|}{} & \multicolumn{1}{c|}{} &  & \multicolumn{1}{c|}{} & \multicolumn{1}{c|}{} & \multicolumn{1}{c|}{} &  \\ \hline
#  & Asian & \multicolumn{1}{c|}{\cellcolor[HTML]{FFCCC9}FALSE} & \multicolumn{1}{c|}{\cellcolor[HTML]{FFCCC9}FALSE} & \multicolumn{1}{c|}{\cellcolor[HTML]{FFCCC9}FALSE} & \cellcolor[HTML]{FFCCC9}FALSE & \multicolumn{1}{c|}{N/A} & \multicolumn{1}{c|}{N/A} & \multicolumn{1}{c|}{N/A} & N/A \\ \cline{2-10} 
#  & Black & \multicolumn{1}{c|}{\cellcolor[HTML]{FFCCC9}FALSE} & \multicolumn{1}{c|}{\cellcolor[HTML]{FFCCC9}FALSE} & \multicolumn{1}{c|}{\cellcolor[HTML]{FFCCC9}FALSE} & \cellcolor[HTML]{FFCCC9}FALSE & \multicolumn{1}{c|}{N/A} & \multicolumn{1}{c|}{N/A} & \multicolumn{1}{c|}{N/A} & N/A \\ \cline{2-10} 
#  & Hispanic & \multicolumn{1}{c|}{\cellcolor[HTML]{FFCCC9}FALSE} & \multicolumn{1}{c|}{\cellcolor[HTML]{FFCCC9}FALSE} & \multicolumn{1}{c|}{\cellcolor[HTML]{FFCCC9}FALSE} & \cellcolor[HTML]{FFCCC9}FALSE & \multicolumn{1}{c|}{N/A} & \multicolumn{1}{c|}{N/A} & \multicolumn{1}{c|}{N/A} & N/A \\ \cline{2-10} 
# \multirow{-4}{*}{Tilt Rotor} & White & \multicolumn{1}{c|}{\cellcolor[HTML]{FFCCC9}FALSE} & \multicolumn{1}{c|}{\cellcolor[HTML]{FFCCC9}FALSE} & \multicolumn{1}{c|}{\cellcolor[HTML]{FFCCC9}FALSE} & \cellcolor[HTML]{FFCCC9}FALSE & \multicolumn{1}{c|}{N/A} & \multicolumn{1}{c|}{N/A} & \multicolumn{1}{c|}{N/A} & N/A \\ \hline
#  & Asian & \multicolumn{1}{c|}{\cellcolor[HTML]{FFCCC9}FALSE} & \multicolumn{1}{c|}{\cellcolor[HTML]{FFCCC9}FALSE} & \multicolumn{1}{c|}{\cellcolor[HTML]{9AFF99}TRUE} & \cellcolor[HTML]{FFCCC9}FALSE & \multicolumn{1}{c|}{\cellcolor[HTML]{FFCCC9}FALSE} & \multicolumn{1}{c|}{\cellcolor[HTML]{FFCCC9}FALSE} & \multicolumn{1}{c|}{\cellcolor[HTML]{9AFF99}TRUE} & \cellcolor[HTML]{FFCCC9}FALSE \\ \cline{2-10} 
#  & Black & \multicolumn{1}{c|}{\cellcolor[HTML]{FFCCC9}FALSE} & \multicolumn{1}{c|}{\cellcolor[HTML]{FFCCC9}FALSE} & \multicolumn{1}{c|}{\cellcolor[HTML]{9AFF99}TRUE} & \cellcolor[HTML]{FFCCC9}FALSE & \multicolumn{1}{c|}{\cellcolor[HTML]{FFCCC9}FALSE} & \multicolumn{1}{c|}{\cellcolor[HTML]{FFCCC9}FALSE} & \multicolumn{1}{c|}{\cellcolor[HTML]{9AFF99}TRUE} & \cellcolor[HTML]{FFCCC9}FALSE \\ \cline{2-10} 
#  & Hispanic & \multicolumn{1}{c|}{\cellcolor[HTML]{9AFF99}TRUE} & \multicolumn{1}{c|}{\cellcolor[HTML]{9AFF99}TRUE} & \multicolumn{1}{c|}{\cellcolor[HTML]{FFCCC9}FALSE} & \cellcolor[HTML]{FFCCC9}FALSE & \multicolumn{1}{c|}{\cellcolor[HTML]{9AFF99}TRUE} & \multicolumn{1}{c|}{\cellcolor[HTML]{9AFF99}TRUE} & \multicolumn{1}{c|}{\cellcolor[HTML]{FFCCC9}FALSE} & \cellcolor[HTML]{FFCCC9}FALSE \\ \cline{2-10} 
# \multirow{-4}{*}{Hexacopter} & White & \multicolumn{1}{c|}{\cellcolor[HTML]{FFCCC9}FALSE} & \multicolumn{1}{c|}{\cellcolor[HTML]{FFCCC9}FALSE} & \multicolumn{1}{c|}{\cellcolor[HTML]{FFCCC9}FALSE} & \cellcolor[HTML]{FFCCC9}FALSE & \multicolumn{1}{c|}{\cellcolor[HTML]{FFCCC9}FALSE} & \multicolumn{1}{c|}{\cellcolor[HTML]{FFCCC9}FALSE} & \multicolumn{1}{c|}{\cellcolor[HTML]{FFCCC9}FALSE} & \cellcolor[HTML]{FFCCC9}FALSE \\ \hline
# \end{tabular}
# \end{table}