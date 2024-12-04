import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np 
import os
import matplotlib.cm as cm
from RCAIDE.Framework.Core import  Data


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
    "$45k to $45k": "B19001010",
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



data_file_tract = pd.read_csv('HC_Data_Tract_CA.csv')

# Define the relevant columns
race_columns = ['Asian','Black','Hispanic','American Indian','Pacific Islander', 'White']

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
    noise_threshold_2 = 55

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
        plt.bar(x[idx], values45[column], width=bar_width, color=base_color, label=None if idx > 0 else  r"> 45 dbA $L_{dn}$")
        plt.bar(x[idx], values2_45[column], width=bar_width, bottom=values45[column], color=lighter_color, label=None if idx > 0 else   r"< 45 dbA $L_{dn}$")

        # Plot bars for 60 dB threshold with hatching
        plt.bar(x[idx] + bar_width, values60[column], width=bar_width, color=base_color, hatch='//', label=None if idx > 0 else  r"> 65 dbA $L_{dn}$")
        plt.bar(x[idx] + bar_width, values2_60[column], width=bar_width, bottom=values60[column], color=lighter_color, hatch='//', label=None if idx > 0 else  r" <65 dbA $L_{dn}$")

    # Add labels, title, and legend
    plt.xticks(x+bar_width/2, col_names)  # Set x-tick labels and positions
    plt.ylabel('Population %')
    plt.ylim(0,120)
    plt.legend(loc='upper center',ncol=2)

    # Save plot
    plot_file_name = os.path.join(folder_path, f'{save_name}{noise_threshold_1}{noise_threshold_2}dB_Plot.png')
    plt.tight_layout()
    plt.savefig(plot_file_name, dpi=1100)
    plt.close()



def split_violin_plot(data_tract, data_columns, noise, x_axis):
    folder_path = 'Violin_Plots'

    # Check if the folder exists, and if not, create it
    if not os.path.isdir(folder_path):
        os.mkdir(folder_path)

    noise_thresholds = [45, 55]
    data=[]
    for threshold in noise_thresholds:
        filtered_tract = data_tract[data_tract[noise] >= threshold]
        val = filtered_tract[data_columns]/data_tract[data_columns]
        val['CA'] = filtered_tract['CA']
        val["L_dn"] = filtered_tract["L_dn"]
        val["L_dnlog(PD)"] = filtered_tract["L_dnlog(PD)"]
        val['threshold'] = threshold
        data.append(val)
    data = pd.concat(data)
    data.to_csv('data.csv',index=False)
    print(data)
    # Create the DataFrame
    df = pd.DataFrame(data)

    # Reshape the data using melt()
    df_melted = df.melt(
        id_vars=["threshold", "L_dn","L_dnlog(PD)", "CA"], 
        value_vars=data_columns, 
        var_name="Race", 
        value_name="Population"
    )
    y_val = ["L_dn","L_dnlog(PD)","CA"]
    for value in y_val:
        plt.figure(figsize=(12, 8))
        sns.violinplot(
            x="Race", y=value, hue="threshold", 
            data=df_melted, split=True, inner="quartile", palette="Set2"
        )
        plt.ylabel("CA")
        plt.tight_layout()

        # Save the plot
        plot_file_name = os.path.join(folder_path, f'violin_plot_{x_axis}{value}.png')
        plt.tight_layout()
        plt.savefig(plot_file_name)
    


    # for threshold in noise_thresholds:
    #     filtered_bg = data_tract[data_tract[noise] >= threshold]

    #     # merged.to_csv('debug.csv', index=False)
    #     percentages = filtered_bg['CA']
    #     percentages['threshold'] = threshold
    #     results.append(percentages)
    #     filtered_bg = 0
    
    # results = pd.concat(results)
    # results.to_csv('violinplot_results.csv', index=False)
    # # Reshape the data for plotting (melt into long format)
    # results_melted = results.melt(id_vars=['threshold'], value_vars=[f"{col}" for col in data_columns], 
    #                                var_name="variable", value_name="Percentage")
    # # Plot the split violin plot using seaborn
    # plt.figure(figsize=(10, 6))
    # sns.violinplot(x='variable', y='Percentage', hue='threshold', data=results_melted,split=True, inner="quart", palette="Set2")
    # plt.title("Race Percentage by Noise Threshold")
    # plt.ylabel("Percentage of Population")
    # plot_file_name = os.path.join(folder_path, f'pls_vplot_Plot.png')
    # plt.tight_layout()
    # plt.savefig(plot_file_name)

#------------------------#
col_names_race = ['Asian','Black','Hisp.','Nat.\nAmer.','Pac.\nIsld.', 'White']
col_names_inc = ['<50K','50k-100k','100k-150k','150k-200k','200k+']
col_names_struct = ['Churches','Medical\nCenters','Schools']

# stackedbplot_race(data_file_tract,race_columns,'L_dn','Race',col_names_race)
# stackedbplot_race(data_file_tract,income_columns,'L_dn','Income',col_names_inc)
# stackedbplot_race(data_file_tract,struct_columns,'L_dn','Structures',col_names_struct)

split_violin_plot(data_file_tract,race_columns,'L_dn','Race')
split_violin_plot(data_file_tract,income_columns,'L_dn','Income')
split_violin_plot(data_file_tract,struct_columns,'L_dn','Structures')


# folder_path = 'Bar_Plots'

# # Check if the folder exists, and if not, create it
# if not os.path.isdir(folder_path):
#     os.mkdir(folder_path)

# noise_thresholds = [45, 65]
# data=[]
# for threshold in noise_thresholds:
#     filtered_tract = data_file_tract[data_file_tract['L_dn'] >= threshold]
#     filtered_tract['threshold'] = threshold
#     data.append(filtered_tract)
# data = pd.concat(data)
# # Create the DataFrame
# df = pd.DataFrame(data)

# # Reshape the data using melt()
# df_melted = df.melt(
#     id_vars=["threshold", "CA"], 
#     value_vars=["Asian", "Black", "Hispanic", "American Indian", "Pacific Islander", "White"], 
#     var_name="Race", 
#     value_name="Population"
# )
# print(df_melted)
# # Plot the split violin plot with CA as the y-axis
# plt.figure(figsize=(12, 8))
# sns.violinplot(
#     x="Race", y="CA", hue="threshold", 
#     data=df_melted, split=True, inner="quartile", palette="Set2"
# )
# plt.ylabel("CA")
# plt.xlabel("Race")
# plt.title("Split Violin Plot: CA by Race and Noise Thresholds")
# plt.tight_layout()

# # Show and save the plot
# plt.savefig("split_violin_plot_ca.png")