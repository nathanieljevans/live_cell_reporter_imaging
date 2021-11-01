import logging
import os
import math
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.colors import ListedColormap
import cell_reporter_list as cr_list
from Set_parameters import source_path

# ---------------------------------------------- SET UP / OPTIONS-------------------------------------------------------
# size of table to display
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

# input cell & reporter
cell_reporter_list = cr_list.cell_reporter
list_index = 0
print("""
cell_reporter_list:""")
for cr in cr_list.cell_reporter:
    print(f"     {list_index} = {cr}")
    list_index += 1

cell_reporter_input = input("input code for cell_reporter(s): ")
if int(cell_reporter_input) > len(cell_reporter_list) - 1:
    print("Warning: wrong cell_reporter input")
    input("please RE-RUN")

# input channel to lead sorting
hmap_df_c = pd.read_csv(f"{source_path}/output_abs/clover/clover_all_cell.csv", index_col="track_index")
hmap_df_m = pd.read_csv(f"{source_path}/output_abs/mscarlet/mscarlet_all_cell.csv", index_col="track_index")
reporter_channel_list = ["clover", "mscarlet"]
print(
    """
    reporter_channel_list:
    0 = clover
    1 = mscarlet
    """
)
lead_ch = input("reporter_channel to lead sorting: ")
lead_df = ""
follow_df = ""
lead_channel_folder = ""
follow_channel_folder = ""
if int(lead_ch) > 1:
    print("wrong color input")
    input("please RE-RUN")
elif int(lead_ch) == 0:
    lead_df = hmap_df_c
    follow_df = hmap_df_m
    lead_channel_folder = "clover"
    follow_channel_folder = "mscarlet"
elif int(lead_ch) == 1:
    lead_df = hmap_df_m
    follow_df = hmap_df_c
    lead_channel_folder = "mscarlet"
    follow_channel_folder = "clover"

cell_reporter = cell_reporter_list[int(cell_reporter_input)]

# Create output folder and subfolders
output_path = f"{source_path}/output_abs_sb_{lead_channel_folder}"
qc_directory = "quality_control"
qc_path = os.path.join(output_path, qc_directory)
os.makedirs(output_path, exist_ok=True)
os.makedirs(qc_path, exist_ok=True)
logging.basicConfig(filename=f'{output_path}/log_file.log', level=logging.INFO)
for rc in reporter_channel_list:
    rc_path = os.path.join(output_path, rc)
    os.makedirs(rc_path, exist_ok=True)
    heatmap_path = os.path.join(rc_path, "heatmap")
    os.makedirs(heatmap_path, exist_ok=True)


# merge cytoplasm and nucleus dataframe based on cn_lookup
def clover_mscarlet_merge(dataframe_1, dataframe_2):
    cm_merge = pd.merge(dataframe_1, dataframe_2, left_index=True, right_index=True, how="inner", suffixes=("_a", "_b"), validate="1:1")
    return cm_merge


def split_dataframe_lead(dataframe):
    dataframe1 = dataframe.iloc[:, 0:len(hmap_df_c.columns)].copy()
    dataframe1.columns = dataframe1.columns.str.replace("_a", "")
    return dataframe1


def split_dataframe_follow(dataframe):
    dataframe1 = dataframe.iloc[:, len(hmap_df_c.columns):len(hmap_df_cm_merge.columns)].copy()
    dataframe1.columns = dataframe1.columns.str.replace("_b", "")
    return dataframe1


def groups(dataframe):
    df_groups = dataframe.groupby("cell__treatment", sort=False)
    return df_groups


# generate heatmap
# commented code to generate one figure with multiple heatmap subplots (problems: not able to arrange the order of subplots)
def heatmap(df_groups, channel, df_open, df_end, colors, minimum, maximum, data_identifier, destination):
    num_of_subplots = len(df_groups)
    n_rows = math.ceil(num_of_subplots / 5)
    fig_width = ""
    fig_height = ""
    n_columns = ""
    if num_of_subplots < 5:
        fig_width = num_of_subplots * 4
        fig_height = 5
        n_columns = num_of_subplots
    elif num_of_subplots > 5:
        fig_width = 20
        fig_height = n_rows * 5
        n_columns = 5
    elif num_of_subplots == 5:
        fig_width = 20
        fig_height = n_rows * 5
        n_columns = 5
    i = 1
    fig = plt.figure(figsize=(fig_width, fig_height), tight_layout=True)  # width x height
    for name, group in df_groups:
        # fig, ax = plt.subplots(figsize=(20,20))# width x height
        # fig = plt.figure(figsize=(15, 15))
        dataframe = group.iloc[:, df_open:group.shape[1] - df_end]
        current_palette = sns.color_palette(colors, n_colors=1000)
        cmap = ListedColormap(sns.color_palette(current_palette).as_hex())
        sns.set_context("paper", font_scale=2)  # rc={"font.size":2,"axes.labelsize":2})
        # fig.subplots_adjust(hspace=0.1, wspace=0.1)
        ax = fig.add_subplot(n_rows, n_columns, i)  # row, column, position
        # sns.heatmap(dataframe, ax=ax, cmap="RdBu_r", center=0, vmin=-1, vmax=1, yticklabels=False)
        title = name.replace(f"{cell_reporter}__", "")
        total_cell_number = group.shape[0]
        sns.heatmap(dataframe, cmap=cmap, vmin=minimum, vmax=maximum, yticklabels=False)
        ax.set_title(f"{title}\nn={total_cell_number}", fontsize="small")
        i += 1
    fig.savefig(f"{output_path}/{channel}/{destination}/{data_identifier}_heatmap.pdf", dpi=300)
    plt.close(fig)
    return


# export divided dataframe in groups
def export_gps(df_groups, channel, data_identifier):
    for name, group in df_groups:
        group.to_csv(f"{output_path}/{channel}/{name}_{data_identifier}.csv", index=True)
    return


# export complete hmap_dataframe
def export_csv(dataframe, channel):
    dataframe.to_csv(f"{output_path}/{channel}/{channel}_all_cell.csv", index=True)
    return


hmap_df_cm_merge = clover_mscarlet_merge(lead_df, follow_df)
hmap_df_lead = split_dataframe_lead(hmap_df_cm_merge)
hmap_df_lead_groups = groups(hmap_df_lead)
hmap_df_follow = split_dataframe_follow(hmap_df_cm_merge)
hmap_df_follow_groups = groups(hmap_df_follow)
heatmap(hmap_df_lead_groups, lead_channel_folder, df_open=0, df_end=3, colors="RdBu_r", data_identifier="reporter", minimum=0, maximum=1, destination="heatmap")
heatmap(hmap_df_follow_groups, follow_channel_folder, df_open=0, df_end=3, colors="RdBu_r", data_identifier="reporter", minimum=0, maximum=1, destination="heatmap")
export_gps(hmap_df_lead_groups, lead_channel_folder, "reporter")
export_gps(hmap_df_follow_groups, follow_channel_folder, "reporter")
export_csv(hmap_df_lead, lead_channel_folder)
export_csv(hmap_df_follow, follow_channel_folder)
