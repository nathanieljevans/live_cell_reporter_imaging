import logging
import os
import numpy as np
from matplotlib.colors import ListedColormap
from sklearn.cluster import KMeans
import math
# from sklearn.manifold import TSNE
# from sklearn.decomposition import PCA
# from sklearn.cluster import DBSCAN
# from sklearn.mixture import GaussianMixture
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import fcluster
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import cell_reporter_list as cr_list
from Set_parameters import treatment_order
from Set_parameters import source_path

print(treatment_order)
# ---------------------------------------------- SET UP / OPTIONS-------------------------------------------------------
# Size of table to display
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

# Input cell & reporter
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

cell_reporter = cell_reporter_list[int(cell_reporter_input)]

clover_reporter = input("""
Signaling pathway of clover reporter: """)
mscarlet_reporter = input("""
Signaling pathway of mScarlet reporter: """)

print("""
INPUT DATA TYPE:
    raw data (abs) = 0
    normalized data (t0norm) = 1""")
data_input = int(input("""Input code for input data type: """))
data_type = ""
folder_name = ""
if data_input == 0:
    data_type = "abs"
    folder_name = "abs"
elif data_input == 1:
    data_type = "t0norm"
    folder_name = "t0norm"

# Create output folder and subfolders
parent_directory = f"{source_path}/output_path_coord_{folder_name}"
plot_directory = "plots"
plot_path = os.path.join(parent_directory, plot_directory)
quadrant_dir = "quadrant"
quandrant_plot_output_path = os.path.join(plot_path, quadrant_dir)
kmeans_dir = "k_mean_clustering"
k_mean_plot_output_path = os.path.join(plot_path, kmeans_dir)
h_cluster_directory = "h_cluster"
h_cluster_output_path = os.path.join(plot_path, h_cluster_directory)
pearson_corr_directory = "pearson_corr"
pearson_corr_output_path = os.path.join(plot_path, pearson_corr_directory)
# tsne_dynamics_dir = "tsne_dynamics"
# tsne_dynamics_output_path = os.path.join(plot_path, tsne_dynamics_dir)
# tsne_dir = "tsne"
# tsne_output_path = os.path.join(plot_path, tsne_dir)
# pca_dynamics_dir = "pca_dynamics"
# pca_dynamics_output_path = os.path.join(plot_path, pca_dynamics_dir)
# pca_dir = "pca"
# pca_output_path = os.path.join(plot_path, pca_dir)
# dbscan_dir = "dbscan"
# dbscan_output_path = os.path.join(plot_path, dbscan_dir)
# gaussian_mix_dir = "gaussian_mix"
# gaussian_mix_output_path = os.path.join(plot_path, gaussian_mix_dir)
# umap_dir = "umap"
# umap_output_path = os.path.join(plot_path, umap_dir)
os.makedirs(plot_path, exist_ok=True)
os.makedirs(quandrant_plot_output_path, exist_ok=True)
os.makedirs(k_mean_plot_output_path, exist_ok=True)
os.makedirs(h_cluster_output_path, exist_ok=True)
os.makedirs(pearson_corr_output_path, exist_ok=True)
# os.makedirs(tsne_dynamics_output_path, exist_ok=True)
# os.makedirs(tsne_output_path, exist_ok=True)
# os.makedirs(pca_dynamics_output_path, exist_ok=True)
# os.makedirs(pca_output_path, exist_ok=True)
# os.makedirs(dbscan_output_path, exist_ok=True)
# os.makedirs(gaussian_mix_output_path, exist_ok=True)
# os.makedirs(umap_output_path, exist_ok=True)

# Read clover and mScarlet reporter activity data from abs folder
hmap_df_c = pd.read_csv(f"{source_path}/output_{data_type}/clover/clover_all_cell.csv", index_col="track_index")
hmap_df_m = pd.read_csv(f"{source_path}/output_{data_type}/mscarlet/mscarlet_all_cell.csv", index_col="track_index")


# Melt pivot table to put all the reporter activity value in one single column
def unpivot_df(dataframe, val_name):
    dataframe2 = dataframe.iloc[:, 0: len(dataframe.columns) - 3]
    dataframe2["track_index"] = dataframe2.index
    dataframe3 = pd.melt(dataframe2, id_vars=["track_index"], var_name=["time"], value_name=val_name)
    return dataframe3


# Merge Clover and mScarlet reporter activity together using unique identifier "track_index_time"
# Drop the row with missing reporter data
# Divide experiment time into 16 slices
def merge_df(dataframe1, dataframe2):
    dataframe1["track_index_time"] = dataframe1["track_index"] + "_" + dataframe1["time"]
    dataframe2["track_index_time"] = dataframe2["track_index"] + "_" + dataframe2["time"]
    dataframe3 = pd.merge(dataframe1, dataframe2.iloc[:, 2:4], how="inner", on="track_index_time", validate="1:1")
    dataframe4 = dataframe3.drop(["track_index_time"], axis=1)
    before_drop = dataframe4.shape[0]
    dataframe4.dropna(axis=0, how="any", subset=["clover", "mscarlet"], inplace=True)
    dataframe4.rename(columns={"clover": clover_reporter, "mscarlet": mscarlet_reporter}, inplace=True)
    after_drop = dataframe4.shape[0]
    num_of_drops = before_drop - after_drop
    percentage_of_drops = round(((num_of_drops / before_drop) * 100), 2)
    logging.info(f"clover_mScarlet_merged_number of drops: {num_of_drops}")
    logging.info(f"clover_mScarlet_merged_percentage of drops: {percentage_of_drops}%")
    dataframe4[["site", "track", "cell__treatment"]] = dataframe4["track_index"].str.split(pat="_", n=2, expand=True)
    dataframe4["time"] = dataframe4["time"].astype("float64")
    # Divide time into 16 time slices
    time_increments = dataframe4["time"].max(axis=0) / 16
    time_list = list(range(17))
    time_slice_list = [i * time_increments for i in time_list]
    dataframe4["time_range"] = pd.cut(dataframe4["time"], time_slice_list, include_lowest=True)
    # Standardize scaling of data for k-means, tsne and pca clustering
    scaled_data = StandardScaler().fit_transform(dataframe4.loc[:, [clover_reporter, mscarlet_reporter, "time"]])
    dataframe5 = pd.DataFrame(scaled_data, index=dataframe4.index,
                              columns=[f"{clover_reporter}_scaled", f"{mscarlet_reporter}_scaled", "time_scaled"])
    dataframe6 = dataframe4.merge(dataframe5, how="inner", left_index=True, right_index=True, validate="1:1")
    dataframe6[["cell", "treatment"]] = dataframe6["cell__treatment"].str.split(pat="__", n=1, expand=True)
    # drop_list = ["0.625nm_dasatinib", "1.25nm_dasatinib", "2.5nm_dasatinib", "5nm_dasatinib", "100nm_dasatinib", "1000nm_dasatinib"]
    # for drop in drop_list:
    # dataframe6.drop(dataframe6.loc[dataframe6["treatment"] == drop].index, inplace=True)
    dataframe6.reset_index(drop=True, inplace=True)
    dataframe6['treatment'] = pd.Categorical(dataframe6['treatment'], categories=treatment_order, ordered=True)
    export_csv(dataframe6, "clover_mscarlet_merged", False)
    return dataframe6


# Divide dataframe into groups based on "cell__treatment" info
def groups(dataframe):
    df_groups = dataframe.groupby("cell__treatment", sort=False)
    return df_groups


# Export divided dataframe in groups
def export_gps(df_groups):
    for name, group in df_groups:
        group.to_csv(f"{parent_directory}/{name}.csv", index=False)
    return


# Export dataframe
def export_csv(dataframe, filename, index_export):
    dataframe.to_csv(f"{parent_directory}/{filename}.csv", index=index_export)
    return


################### QUADRANT #########################################################################################
# Generate threshold values to divide scatterplot into 4 parts based on "untreated" median
def quadrant_lines_threshold(dataframe):
    # dataframe[["cell", "treatment"]] = dataframe["cell__treatment"].str.split(pat="__", n=1, expand=True)
    # Use "contain "untreated"" to include multiple cell lines (e.g. OVAR__heya8_untreated, OVAR__ovcar8_untreated)
    untreated_dataframe = dataframe[dataframe["treatment"].str.contains("untreated")]
    ## ------------------------GENERATE TIME-DEPENDENT THRESHOLD VALUES -----------------------------------
    # clover_mean = untreated_dataframe.groupby("time_range")[clover_reporter].mean()
    # mscarlet_mean = untreated_dataframe.groupby("time_range")[mscarlet_reporter].mean()
    # clover_std = untreated_dataframe.groupby("time_range")[clover_reporter].std()
    # mscarlet_std = untreated_dataframe.groupby("time_range")[mscarlet_reporter].std()
    # merged_mean = pd.merge(clover_mean, mscarlet_mean, how="inner", on="time_range", validate="1:1")
    # merged_mean.rename(columns={clover_reporter: f"{clover_reporter}_mean", mscarlet_reporter: f"{mscarlet_reporter}_mean"}, inplace=True)
    # merged_std = pd.merge(clover_std, mscarlet_std, how="inner", on="time_range", validate="1:1")
    # merged_std.rename(columns={clover_reporter: f"{clover_reporter}_std", mscarlet_reporter: f"{mscarlet_reporter}_std"}, inplace=True)
    # merged_threshold = pd.merge(merged_mean, merged_std, how="inner", on="time_range", validate="1:1")
    # merged_threshold[f"{clover_reporter}_threshold"] = merged_threshold[f"{clover_reporter}_mean"] + merged_threshold[f"{clover_reporter}_std"]*0
    # merged_threshold[f"{mscarlet_reporter}_threshold"] = merged_threshold[f"{mscarlet_reporter}_mean"] + merged_threshold[f"{mscarlet_reporter}_std"]*0
    return untreated_dataframe


# Generate a dataframe of cell distribution per time index based on reporter coordination
def quadrant_stats_df(df_groups, time_subgroup_index):
    rows = []
    for name, group in df_groups:
        dataframe = group
        time_subgroup = dataframe.groupby(time_subgroup_index)
        for time_name, time_group in time_subgroup:
            dataframe1 = time_group
            total_count = len(dataframe1.index)
            # Output # of rows that matches the rules
            upper_left_count = len(dataframe1[(dataframe1[f"{clover_reporter}"] < clover_threshold) & (
                    dataframe1[f"{mscarlet_reporter}"] > mscarlet_threshold)].index)
            upper_right_count = len(dataframe1[(dataframe1[f"{clover_reporter}"] > clover_threshold) & (
                    dataframe1[f"{mscarlet_reporter}"] > mscarlet_threshold)].index)
            bottom_left_count = len(dataframe1[(dataframe1[f"{clover_reporter}"] < clover_threshold) & (
                    dataframe1[f"{mscarlet_reporter}"] < mscarlet_threshold)].index)
            bottom_right_count = len(dataframe1[(dataframe1[f"{clover_reporter}"] > clover_threshold) & (
                    dataframe1[f"{mscarlet_reporter}"] < mscarlet_threshold)].index)
            upper_left_percentage = round(((upper_left_count / total_count) * 100), 2)
            upper_right_percentage = round(((upper_right_count / total_count) * 100), 2)
            bottom_left_percentage = round(((bottom_left_count / total_count) * 100), 2)
            bottom_right_percentage = round(((bottom_right_count / total_count) * 100), 2)
            # Add a row of assigned values
            rows.append([name, time_name, total_count, upper_left_count, upper_right_count, bottom_left_count,
                         bottom_right_count, upper_left_percentage, upper_right_percentage, bottom_left_percentage,
                         bottom_right_percentage])
    df = pd.DataFrame(rows, columns=["cell__treatment", "time", "total_count",
                                     f"{clover_reporter}_low__{mscarlet_reporter}_high_count",
                                     f"{clover_reporter}_high__{mscarlet_reporter}_high_count",
                                     f"{clover_reporter}_low__{mscarlet_reporter}_low_count",
                                     f"{clover_reporter}_high__{mscarlet_reporter}_low_count",
                                     f"{clover_reporter}_low__{mscarlet_reporter}_high_percentage",
                                     f"{clover_reporter}_high__{mscarlet_reporter}_high_percentage",
                                     f"{clover_reporter}_low__{mscarlet_reporter}_low_percentage",
                                     f"{clover_reporter}_high__{mscarlet_reporter}_low_percentage"])
    return df


# Assign quadrant cluster id to each cell 0) upper left, 1) upper right, 2) lower left, 3) lower right
def quadrant_membership(dataframe):
    dataframe1 = dataframe.copy()
    dataframe1["cluster_id"] = ""
    # upper left
    dataframe1.loc[(dataframe1[f"{clover_reporter}"] < clover_threshold) & (
            dataframe1[f"{mscarlet_reporter}"] > mscarlet_threshold), "cluster_id"] = 0
    # upper right
    dataframe1.loc[(dataframe1[f"{clover_reporter}"] > clover_threshold) & (
            dataframe1[f"{mscarlet_reporter}"] > mscarlet_threshold), "cluster_id"] = 1
    # bottom left
    dataframe1.loc[(dataframe1[f"{clover_reporter}"] < clover_threshold) & (
            dataframe1[f"{mscarlet_reporter}"] < mscarlet_threshold), "cluster_id"] = 2
    # bottom right
    dataframe1.loc[(dataframe1[f"{clover_reporter}"] > clover_threshold) & (
            dataframe1[f"{mscarlet_reporter}"] < mscarlet_threshold), "cluster_id"] = 3
    dataframe1.loc[(dataframe1[f"{clover_reporter}"] == clover_threshold) & (
            dataframe1[f"{mscarlet_reporter}"] > mscarlet_threshold), "cluster_id"] = 0
    dataframe1.loc[(dataframe1[f"{clover_reporter}"] == clover_threshold) & (
            dataframe1[f"{mscarlet_reporter}"] < mscarlet_threshold), "cluster_id"] = 3
    dataframe1.loc[(dataframe1[f"{clover_reporter}"] < clover_threshold) & (
            dataframe1[f"{mscarlet_reporter}"] == mscarlet_threshold), "cluster_id"] = 2
    dataframe1.loc[(dataframe1[f"{clover_reporter}"] > clover_threshold) & (
            dataframe1[f"{mscarlet_reporter}"] == mscarlet_threshold), "cluster_id"] = 1
    return dataframe1


# Create 16 scatter plots in one figure with quadrant lines
def quad_path_coordination_scatterplot(df_groups):
    x_min = df_groups[clover_reporter].min(axis=1).min()
    x_max = df_groups[clover_reporter].max(axis=1).min()
    y_min = df_groups[mscarlet_reporter].min(axis=1).min()
    y_max = df_groups[mscarlet_reporter].max(axis=1).min()
    for name, group in df_groups:
        dataframe = group
        sns.set_context("talk")
        ax = sns.relplot(data=dataframe, x=clover_reporter, y=mscarlet_reporter, col="time_range", col_wrap=4,
                         legend="full", s=50, alpha=0.5)
        # To allow iteration of subplots, use "axes.flatten()"
        axes = ax.axes.flatten()
        # Add threshold line and percentage of cell distribution in each quarter for all subplots
        quadrant_group_slice = quadrant_stats_time_range_results.loc[
            quadrant_stats_time_range_results['cell__treatment'] == name]
        time_range_count = 0
        for i in axes:
            # Add threshold line
            i.plot([x_min + 0.1, x_max - 0.1], [mscarlet_threshold, mscarlet_threshold], ls='--', linewidth=2,
                   color='red')
            i.plot([clover_threshold, clover_threshold], [y_min + 0.1, y_max - 0.1], ls='--', linewidth=2, color='red')
            # Add percentage of cell distribution
            up_left = quadrant_group_slice.iloc[time_range_count][
                f"{clover_reporter}_low__{mscarlet_reporter}_high_percentage"]
            up_right = quadrant_group_slice.iloc[time_range_count][
                f"{clover_reporter}_high__{mscarlet_reporter}_high_percentage"]
            bot_left = quadrant_group_slice.iloc[time_range_count][
                f"{clover_reporter}_low__{mscarlet_reporter}_low_percentage"]
            bot_right = quadrant_group_slice.iloc[time_range_count][
                f"{clover_reporter}_high__{mscarlet_reporter}_low_percentage"]
            i.text(x_min + 0.1, y_max - 0.1, f"{up_left}%", fontsize=15, color="red", weight='semibold')
            i.text(x_max - 0.15, y_max - 0.1, f"{up_right}%", fontsize=15, color="red", weight='semibold')
            i.text(x_min + 0.1, y_min + 0.1, f"{bot_left}%", fontsize=15, color="red", weight='semibold')
            i.text(x_max - 0.15, y_min + 0.1, f"{bot_right}%", fontsize=15, color="red", weight='semibold')
            time_range_count += 1
        title = name.replace(f"{cell_reporter}__", "")
        plt.savefig(f"{quandrant_plot_output_path}/quad_{title}.pdf", dpi=300)
        plt.close()
    return


# Generate a line graph to show the trend of the 4 categories of reporter coordination
# Generate a line graph to show the trend of the 4 categories of reporter coordination after normalized to t0
# The time axis (x-axis) of the line graph use all the image timepoints instead of time-range created
def quadrant_lineplot(dataframe):
    dataframe1 = pd.melt(dataframe, id_vars=["cell__treatment", "time"],
                         value_vars=[f"{clover_reporter}_low__{mscarlet_reporter}_high_percentage",
                                     f"{clover_reporter}_high__{mscarlet_reporter}_high_percentage",
                                     f"{clover_reporter}_low__{mscarlet_reporter}_low_percentage",
                                     f"{clover_reporter}_high__{mscarlet_reporter}_low_percentage"],
                         var_name="quadrant_class", value_name="quadrant_percentage")
    dataframe1[["cell", "treatment"]] = dataframe1["cell__treatment"].str.split(pat="__", n=1, expand=True)
    dataframe1['treatment'] = pd.Categorical(dataframe1['treatment'], categories=treatment_order,
                                             ordered=True)
    dataframe1_t0 = dataframe1.loc[dataframe1["time"] == 0]
    dataframe1_t0_2 = dataframe1_t0.rename(columns={"quadrant_percentage": "quadrant_percentage_t0"})
    dataframe1_add_t0 = pd.merge(dataframe1,
                                 dataframe1_t0_2[["cell__treatment", "quadrant_class", "quadrant_percentage_t0"]],
                                 on=["cell__treatment", "quadrant_class"], how="left", validate="m:1")
    dataframe1_add_t0["t0_normalized_percentage"] = dataframe1_add_t0["quadrant_percentage"] - dataframe1_add_t0[
        "quadrant_percentage_t0"]
    treatment_n = dataframe1["treatment"].nunique()
    if treatment_n < 5:
        col_n = treatment_n
    else:
        col_n = 5
    sns.set_context("talk")
    sns.relplot(data=dataframe1, x="time", y="quadrant_percentage", hue="quadrant_class", kind="line",
                col="treatment", col_wrap=col_n, legend="full", markers=True)
    plt.savefig(f"{quandrant_plot_output_path}/quadrant_line.pdf", dpi=600)
    plt.close()
    sns.relplot(data=dataframe1_add_t0, x="time", y="t0_normalized_percentage", hue="quadrant_class", kind="line",
                col="treatment", col_wrap=col_n, legend="full", markers=True)
    plt.savefig(f"{quandrant_plot_output_path}/normalized_quadrant_line.pdf", dpi=600)
    plt.close()
    return


################################ K-MEAN CLUSTERING ####################################################################
# Create elbow plot per treatment to define treatment-based cluster number "k"
# Perform k-mean clustering for all samples disregarding treatment groups
# For k-means clustering with raw data, post-tsne and post-pca
# Feed the algorithm with scaled data
# Plot the clustering with scaled data
########################################################################################################################

def clustering_on_df(dataframe, data_origin):
    x_value = ""
    y_value = ""
    plot_output_pathway = ""
    sub_filename = ""
    dataframe_selected = ""
    pathway = ""
    if data_origin == "k_mean":
        dataframe_selected = dataframe
        x_value = f"{clover_reporter}_scaled"
        y_value = f"{mscarlet_reporter}_scaled"
        plot_output_pathway = k_mean_plot_output_path
        sub_filename = "k_mean"
        pathway = ""
    # elif data_origin == "tsne_dynamics":
    #     x_value = "tsne_x"
    #     y_value = "tsne_y"
    #     plot_output_pathway = tsne_dynamics_output_path
    #     sub_filename = "tsne_dynamics"
    #     dataframe_selected = dataframe.loc[dataframe["perplexity_value"] == dynamics_best_perplexity]
    #     pathway = tsne_dynamics_path
    # elif data_origin == "pca_dynamics":
    #     x_value = k_pca_dy_x
    #     y_value = k_pca_dy_y
    #     plot_output_pathway = pca_dynamics_output_path
    #     sub_filename = "pca_dynamics"
    #     dataframe_selected = dataframe
    #     pathway = pca_dynamics_path
    # elif data_origin == "pca":
    #     x_value = k_pca_x
    #     y_value = k_pca_y
    #     plot_output_pathway = pca_output_path
    #     sub_filename = "pca"
    #     dataframe_selected = dataframe
    #     pathway = ""
    else:
        print("Define data origin")
    # create empty dictionary to keep all the created dataframe with cluster id-labelled objects
    dataframe1 = dataframe_selected.loc[:, [x_value, y_value]]
    # Elbow plot to determine the number of cluster for each treatment
    # Fitting multiple k-means algorithms and storing the values in an empty list
    sum_of_squared_distance = []
    for cluster in range(1, 20):
        kmeans = KMeans(n_clusters=cluster, init='k-means++')
        kmeans.fit(dataframe1)
        sum_of_squared_distance.append(kmeans.inertia_)
    # Converting the results into a dataframe and plotting them
    elbow_plot_dataframe = pd.DataFrame({'Cluster': range(1, 20), 'sum_of_squared_distance': sum_of_squared_distance})
    plt.figure(figsize=(12, 6))
    plt.plot(elbow_plot_dataframe['Cluster'], elbow_plot_dataframe['sum_of_squared_distance'], marker='o')
    plt.xlabel('Number of clusters')
    plt.ylabel('Inertia')
    working_fig = plt.gcf()
    plt.show()
    working_fig.savefig(f"{plot_output_pathway}/{pathway}_{sub_filename}_elbow.pdf", dpi=300)
    # Input the number of cluster for each treatment
    k = int(input(f"number of clusters of {sub_filename} based on {pathway}: "))
    plt.close()
    # Run k mean clustering algorithm
    kmeans = KMeans(n_clusters=k, init="k-means++", max_iter=1000, n_init=10, random_state=0)
    y_kmeans = kmeans.fit_predict(dataframe1)
    dataframe1["cluster_id"] = y_kmeans
    dataframe_center = pd.DataFrame(kmeans.cluster_centers_)
    dataframe_center.columns = [x_value, y_value]
    # Visualising the clusters
    color_list = sns.color_palette(cluster_color_palette, n_colors=k)
    sns.set_context("paper")
    ax = sns.scatterplot(data=dataframe1, x=x_value, y=y_value, s=30, legend="full", hue="cluster_id", alpha=0.3,
                         palette=color_list)
    ax.legend(loc='upper center', bbox_to_anchor=(-0.5, 1.2, 2, 0.1), ncol=k + 1, fontsize="small")
    # Plot the centroid. This time we're going to use the cluster centres
    # Attribute that returns here the coordinates of the centroid.
    sns.scatterplot(data=dataframe_center, x=x_value, y=y_value, s=50, color='yellow', legend=False)
    plt.savefig(f"{plot_output_pathway}/{sub_filename}_{pathway}_cluster.pdf", dpi=300)
    plt.close()
    # Add "cluster_id" to each cell at different time_point
    dataframe2 = pd.merge(dataframe_selected, dataframe1["cluster_id"], how="inner", left_index=True, right_index=True,
                          validate="1:1")
    return dataframe2


# Add a "cluster_id_t0" column to assign the starting (t0) cluster id to each cell at different time
def cell_cluster_id_t0_df(dataframe, data_origin):
    sub_filename = ""
    if data_origin == "k_mean":
        sub_filename = "k_mean_all_cell_cluster_id"
    elif data_origin == "tsne":
        sub_filename = "tsne_all_cell_cluster_id"
    elif data_origin == "dbscan":
        sub_filename = "dbscan_all_cell_cluster_id"
    elif data_origin == "quadrant":
        sub_filename = "quad_all_cell_cluster_id"
    elif data_origin == "gaussian_mix":
        sub_filename = "gaussian_mix_all_cell_cluster_id"
    else:
        print("Define data origin")
    dataframe_t0 = dataframe.loc[dataframe["time"] == 0]
    dataframe_t0_2 = dataframe_t0.rename(columns={"cluster_id": "cluster_id_t0"})
    dataframe_cluster_id_t0 = pd.merge(dataframe, dataframe_t0_2[["track_index", "cluster_id_t0"]], on="track_index",
                                       how="left", validate="m:1")
    export_csv(dataframe_cluster_id_t0, sub_filename, False)
    return dataframe_cluster_id_t0


# Create 16 scatter plots in one figure WITHOUT quadrant lines
# Use raw (non-scaled) data to plot the scatterplots
def path_coordination_scatterplot(dataframe, palette, data_origin):
    x_value = ""
    y_value = ""
    plot_output_pathway = ""
    sub_filename = ""
    hue = ""
    time_div = ""
    if data_origin == "raw":
        x_value = clover_reporter
        y_value = mscarlet_reporter
        plot_output_pathway = plot_path
        hue = None
        sub_filename = ""
        time_div = "time_range"
    elif data_origin == "k_mean":
        x_value = clover_reporter
        y_value = mscarlet_reporter
        plot_output_pathway = k_mean_plot_output_path
        hue = "cluster_id_t0"
        sub_filename = "k_mean_cluster_id"
        time_div = "time_range"
    # elif data_origin == "tsne":
    #     x_value = clover_reporter
    #     y_value = mscarlet_reporter
    #     plot_output_pathway = tsne_dynamics_output_path
    #     hue = "cluster_id_t0"
    #     sub_filename = "tsne_cluster_id"
    #     time_div = "time_range"
    # elif data_origin == "dbscan":
    # x_value = clover_reporter
    # y_value = mscarlet_reporter
    # plot_output_pathway = dbscan_output_path
    # hue = "cluster_id_t0"
    # sub_filename = "dbscan_cluster_id"
    # time_div = "time_slice"
    # elif data_origin == "gaussian_mix":
    # x_value = clover_reporter
    # y_value = mscarlet_reporter
    # plot_output_pathway = gaussian_mix_output_path
    # hue = "cluster_id_t0"
    # sub_filename = "gau_mix_cluster_id"
    # time_div = "time_range"
    else:
        print("Define data origin")
    df_groups = groups(dataframe)
    for name, group in df_groups:
        if "cluster_id" in dataframe.columns:
            color_number = group[hue].nunique()
        else:
            color_number = 1
        color_list = sns.color_palette(palette, n_colors=color_number)
        sns.set_context("talk")
        sns.relplot(data=group, x=x_value, y=y_value, col=time_div, col_wrap=4,
                    legend="full", s=50, alpha=0.8, hue=hue, palette=color_list)
        title = name.replace(f"{cell_reporter}__", "")
        plt.savefig(f"{plot_output_pathway}/{title}_{sub_filename}_scatter.pdf", dpi=300)
        plt.close()
    return


# Switch clustered dataframe structure to pivot table
def cluster_switch_melt_to_pivot_df(dataframe, data_origin):
    df_id_sort = dataframe.sort_values(by="cluster_id", ascending=True)
    cluster_label = df_id_sort["cluster_id"].unique().tolist()
    sub_filename = ""
    if data_origin == "k_mean":
        sub_filename = "k_mean_"
    elif data_origin == "tsne":
        sub_filename = "tsne_"
    elif data_origin == "dbscan":
        sub_filename = "dbscan_"
    elif data_origin == "quadrant":
        sub_filename = "quad_"
    elif data_origin == "gaussian_mix":
        sub_filename = "gaussian_mix_"
    else:
        print("Define data origin")
    dataframe_matrix = dataframe.pivot(index="track_index", columns="time", values="cluster_id")
    dataframe_matrix["track_index"] = dataframe_matrix.index
    dataframe_matrix[["site", "track", "cell__treatment"]] = dataframe_matrix["track_index"].str.split(pat="_", n=2,
                                                                                                       expand=True)
    dataframe_matrix.drop(["site", "track"], axis=1, inplace=True)
    dataframe_matrix[["cell", "treatment"]] = dataframe_matrix["cell__treatment"].str.split(pat="__", n=1, expand=True)
    dataframe_matrix_2 = dataframe_matrix.iloc[:, 0:len(dataframe_matrix.columns) - 2]
    # Report the most frequent cluster_id for each cell
    dataframe_matrix["mode"] = dataframe_matrix_2.mode(axis=1).iloc[:, 0]
    dataframe_matrix['treatment'] = pd.Categorical(dataframe_matrix['treatment'], categories=treatment_order,
                                                   ordered=True)
    # Sort cell tracks based on 1) treatment, 2) cluster_id at t0, 3) most frequent cluster_id afterwards
    dataframe_matrix.sort_values(by=["treatment", 0.0, "mode"], ascending=True, inplace=True)
    export_csv(dataframe_matrix, f"{sub_filename}pivot_all_cell_cluster_id", True)
    return dataframe_matrix, cluster_label


# Create barplot to show cluster switch frequency for each starting cluster
def cluster_switch_frequency(dataframe, data_origin):
    sub_filename = ""
    plot_output_path = ""
    # cluster_id_t0_dataframe = ""
    if data_origin == "kmeans":
        sub_filename = "kmeans"
        plot_output_path = k_mean_plot_output_path
        # cluster_id_t0_dataframe = k_cluster_id_t0_dataframe
    # elif data_origin == "dbscan":
    # sub_filename = "dbscan"
    # plot_output_path = dbscan_output_path
    # cluster_id_t0_dataframe = db_cluster_id_t0_dataframe
    elif data_origin == "quadrant":
        sub_filename = "quad"
        plot_output_path = quandrant_plot_output_path
        # cluster_id_t0_dataframe = q_cluster_id_t0_dataframe
    # elif data_origin == "gaussian_mix":
    # sub_filename = "gaussian_mix"
    # plot_output_path = gaussian_mix_output_path
    # cluster_id_t0_dataframe = gau_cluster_id_t0_dataframe
    else:
        print("Define data origin")
    # Forward fill empty cells so they won't falsely regard as cluster switch
    dataframe1 = dataframe.fillna(method="ffill", axis=1, inplace=False)
    # Make a copy so the changes downwards won't affect the heatmap plots
    dataframe2 = dataframe1.copy()
    working_df = dataframe2.iloc[:, 0:len(dataframe2.columns) - 5]
    # Count total number of cluster switch per cell track
    df_diff = working_df.diff(axis=1) != 0
    df_diff.drop(0.0, axis=1, inplace=True)
    df_diff = df_diff.astype(int)  # change True & False to 1 & 0
    switch_counts = df_diff.apply(pd.Series.value_counts, axis=1, sort=False, normalize=False)
    switch_counts.reset_index(drop=False, inplace=True)
    switch_counts[["site", "track", "cell__treatment"]] = switch_counts["track_index"].str.split(pat="_", n=2,
                                                                                                 expand=True)
    switch_counts[["cell", "treatment"]] = switch_counts["cell__treatment"].str.split(pat="__", n=1, expand=True)
    # dataframe_id = pd.merge(switch_counts, cluster_id_t0_dataframe.loc[:, ["track_index", "cluster_id_t0"]], how="inner", on="track_index", validate="1:m")
    dataframe_id = pd.merge(switch_counts, dataframe.loc[:, ["track_index", "mode"]],
                            how="inner", left_on="track_index", right_index=True, validate="1:1")
    # dataframe_id.drop_duplicates(keep='first', inplace=True)
    # Reset the wrong index after drop duplicates
    dataframe_id.reset_index(drop=True, inplace=True)
    dataframe_id.rename(columns={0: "unchange_frequency_per_cell", 1: "switching_frequency_per_cell"}, inplace=True)
    treatment_n = dataframe_id["treatment"].nunique()
    if treatment_n < 5:
        col_n = treatment_n
    else:
        col_n = 5
    sns.set_context("talk")
    sns.catplot(data=dataframe_id, x="mode", y="switching_frequency_per_cell", kind="box", col="treatment",
                col_wrap=col_n)
    plt.savefig(f"{plot_output_path}/{sub_filename}_cluster_switch_frequency_box.pdf", dpi=300)
    plt.close()
    sns.catplot(data=dataframe_id, x="mode", y="switching_frequency_per_cell", kind="swarm", col="treatment",
                col_wrap=col_n)
    plt.savefig(f"{plot_output_path}/{sub_filename}_cluster_switch_frequency_swarm.pdf", dpi=300)
    plt.close()
    export_csv(dataframe_id, f"{sub_filename}_cluster_switch_frequency_dataframe", False)
    weight_list = range(0, len(df_diff.columns))
    df_diff_weighted = df_diff.copy()
    df_diff_weighted *= np.array(weight_list)
    df_diff_weighted["sum"] = df_diff_weighted.sum(axis=1)
    df_diff["sum"] = df_diff_weighted["sum"]
    # df_diff["sum"] = df_diff.sum(axis=1)
    df_diff.reset_index(drop=False, inplace=True)
    df_diff[["site", "track", "cell__treatment"]] = df_diff["track_index"].str.split(pat="_", n=2,
                                                                                     expand=True)
    df_diff[["cell", "treatment"]] = df_diff["cell__treatment"].str.split(pat="__", n=1, expand=True)
    df_diff['treatment'] = pd.Categorical(df_diff['treatment'], categories=treatment_order,
                                          ordered=True)
    df_diff.sort_values(by=["treatment", "sum"], inplace=True)
    df_groups = groups(df_diff)
    # fig = plt.figure(figsize=(50, 25))
    fig = plt.figure(figsize=(20, 10), tight_layout=True)
    i = 1
    # for c_treatment in cell_treatment_list:
    for name, group in df_groups:
        # work_df = dataframe.loc[dataframe["cell__treatment"] == c_treatment]
        # df = work_df.iloc[:, 0:len(work_df.columns) - 5]
        df = group.iloc[:, 1:len(group.columns) - 6]
        sns.set_context("talk", font_scale=1)
        ax = fig.add_subplot(2, 5, i)  # row, column, position
        title = name.replace(f"{cell_reporter}__", "")
        total_cell_number = group.shape[0]
        sns.heatmap(data=df, ax=ax, yticklabels=False)
        # colorbar = ax.collections[0].colorbar
        # r = colorbar.vmax - colorbar.vmin
        # colorbar.set_ticks([colorbar.vmin + r / cluster_number * (0.5 + i) for i in range(cluster_number)])
        # colorbar.set_ticklabels(cluster_label)
        ax.set_title(f"{title}\nn={total_cell_number}", fontsize="small")
        i += 1
    fig.savefig(f"{plot_output_path}/{sub_filename}_cluster_switch_time_heatmap.pdf", dpi=300)
    plt.close(fig)
    return


# Create cluster switching trajectory heatmap for individual cell
def cluster_switch_trajectory_heatmap(dataframe, data_origin):
    plot_output_path = ""
    cluster_label = ""
    if data_origin == "k_mean":
        plot_output_path = k_mean_plot_output_path
        cluster_label = k_cluster_label
    # elif data_origin == "dbscan":
    # plot_output_path = dbscan_output_path
    # cluster_label = db_cluster_label
    elif data_origin == "quadrant":
        plot_output_path = quandrant_plot_output_path
        cluster_label = q_cluster_label
    # elif data_origin == "gaussian_mix":
    # plot_output_path = gaussian_mix_output_path
    # cluster_label = gau_cluster_label
    else:
        print("Define data origin")
    # cell_treatment_list = dataframe["cell__treatment"].unique().tolist()
    df_groups = groups(dataframe)
    # fig = plt.figure(figsize=(50, 25))
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
    fig = plt.figure(figsize=(fig_width, fig_height), tight_layout=True)
    i = 1
    # for c_treatment in cell_treatment_list:
    for name, group in df_groups:
        # work_df = dataframe.loc[dataframe["cell__treatment"] == c_treatment]
        # df = work_df.iloc[:, 0:len(work_df.columns) - 5]
        df = group.iloc[:, 0:len(group.columns) - 5]
        cluster_number = len(cluster_label)
        current_palette = sns.color_palette(cluster_color_palette, n_colors=cluster_number)
        sns.set_context("talk", font_scale=1)
        ax = fig.add_subplot(n_rows, n_columns, i)  # row, column, position
        title = name.replace(f"{cell_reporter}__", "")
        total_cell_number = group.shape[0]
        sns.heatmap(data=df, ax=ax, cmap=current_palette, yticklabels=False)
        colorbar = ax.collections[0].colorbar
        r = colorbar.vmax - colorbar.vmin
        colorbar.set_ticks([colorbar.vmin + r / cluster_number * (0.5 + i) for i in range(cluster_number)])
        colorbar.set_ticklabels(cluster_label)
        ax.set_title(f"{title}\nn={total_cell_number}", fontsize="small")
        i += 1
    fig.savefig(f"{plot_output_path}/cluster_switch_heatmap.pdf", dpi=300)
    plt.close(fig)
    return


# Find out the population percentage of each cluster per timepoint per treatment group
# Create dataframe for the lineplot and create lineplot for cluster populations
# Normalized to t0 percentage for better dynamic comparison
# Create lineplot using the t0 normalized data
def cluster_population_lineplot(dataframe_gps, data_origin):
    plot_output_path = ""
    sub_filename = ""
    if data_origin == "k_mean":
        plot_output_path = k_mean_plot_output_path
        sub_filename = "k_mean"
    # elif data_origin == "tsne":
    #     plot_output_path = tsne_dynamics_output_path
    #     sub_filename = "tsne"
    # elif data_origin == "dbscan":
    # plot_output_path = dbscan_output_path
    # sub_filename = "dbscan"
    # elif data_origin == "gaussian_mix":
    # plot_output_path = gaussian_mix_output_path
    # sub_filename = "gaussian_mix"
    else:
        print("Define data origin")
    rows = []
    for name, group in dataframe_gps:
        dataframe = group
        time_list = dataframe["time"].unique().tolist()
        dataframe_t0 = dataframe.loc[dataframe["time"] == 0]
        for time_id in time_list:
            working_dataframe = dataframe.loc[dataframe["time"] == time_id]
            total = len(working_dataframe.index)
            cluster_id_list = working_dataframe["cluster_id"].unique().tolist()
            for working_cluster_id in cluster_id_list:
                working_cluster_dataframe_slice = working_dataframe.loc[
                    working_dataframe["cluster_id"] == working_cluster_id]
                working_cluster_counts = len(working_cluster_dataframe_slice.index)
                working_cluster_percentage = working_cluster_counts / total * 100
                working_cluster_percentage_t0 = len(
                    dataframe_t0.loc[dataframe_t0["cluster_id"] == working_cluster_id]) / total * 100
                rows.append([working_cluster_id, working_cluster_counts, total, working_cluster_percentage,
                             working_cluster_percentage_t0, time_id, name])
    df = pd.DataFrame(rows, columns=["cluster_id", "cluster_id_counts", "total_counts", "cluster_percentage",
                                     "cluster_percentage_t0", "time", "cell__treatment"])
    df["cluster_percentage_normalized"] = df["cluster_percentage"] - df["cluster_percentage_t0"]
    df[["cell", "treatment"]] = df["cell__treatment"].str.split(pat="__", n=1, expand=True)
    df['treatment'] = pd.Categorical(df['treatment'], categories=treatment_order, ordered=True)
    export_csv(df, f"{sub_filename}_cluster_switch_lineplot_data", False)
    cluster_number = df["cluster_id"].nunique()
    color_list = sns.color_palette("muted", n_colors=cluster_number)
    treatment_n = df["treatment"].nunique()
    if treatment_n < 5:
        col_n = treatment_n
    else:
        col_n = 5
    sns.set_context("talk")
    # plot lineplot using abs values
    sns.relplot(data=df, x="time", y="cluster_percentage", hue="cluster_id", kind="line", palette=color_list,
                col="treatment", col_wrap=col_n, legend="full", markers=True)
    plt.savefig(f"{plot_output_path}/{sub_filename}_cluster_population_lineplot.pdf", dpi=600)
    plt.close()
    # plot lineplot using t0 normalized values
    sns.relplot(data=df, x="time", y="cluster_percentage_normalized", hue="cluster_id", kind="line",
                palette=color_list, col="treatment", col_wrap=col_n, legend="full", markers=True)
    plt.savefig(f"{plot_output_path}/{sub_filename}_cluster_population_lineplot_normalized.pdf", dpi=300)
    plt.close()
    return df


############ Hierarchical clustering ###################################################################################
def fancy_dendrogram(*args, **kwargs):
    max_d = kwargs.pop('max_d', None)
    if max_d and 'color_threshold' not in kwargs:
        kwargs['color_threshold'] = max_d
    annotate_above = kwargs.pop('annotate_above', 0)

    ddata = dendrogram(*args, **kwargs)

    if not kwargs.get('no_plot', False):
        plt.title('Hierarchical Clustering Dendrogram (truncated)')
        plt.xlabel('sample index or (cluster size)')
        plt.ylabel('distance')
        for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            if y > annotate_above:
                plt.plot(x, y, 'o', c=c)
                plt.annotate("%.3g" % y, (x, y), xytext=(0, -5),
                             textcoords='offset points',
                             va='top', ha='center')
        if max_d:
            plt.axhline(y=max_d, c='k')
    return ddata


def hierarchical_clustering(dataframe, reporter2):
    lead_reporter = ""
    sec_reporter = ""
    lead_color = ""
    sec_color = ""
    if reporter2 == clover_reporter:
        lead_reporter = clover_reporter
        sec_reporter = mscarlet_reporter
        lead_color = "Greens"
        sec_color = "Reds"
    elif reporter2 == mscarlet_reporter:
        lead_reporter = mscarlet_reporter
        sec_reporter = clover_reporter
        lead_color = "Reds"
        sec_color = "Greens"
    else:
        print("PLEASE define reporter!")
    dataframe.sort_values(by="treatment", inplace=True)
    lead_reporter_pivot = dataframe.pivot(index="track_index", columns="time", values=lead_reporter)
    z = linkage(lead_reporter_pivot, "ward")
    c, coph_dists = cophenet(z, pdist(lead_reporter_pivot))
    print(c)
    plt.figure(figsize=(25, 10), tight_layout=True)
    sns.set_context("paper")
    plt.xlabel("track_index")
    plt.ylabel("distance")
    fancy_dendrogram(z, leaf_rotation=90, leaf_font_size=6, labels=lead_reporter_pivot.index,
                     above_threshold_color="#AAAAAA", truncate_mode='lastp', p=30, show_contracted=True,
                     annotate_above=6)
    plt.show()
    plt.close()
    max_d = float(input(f"{lead_reporter}_cluster cutoff: "))
    fancy_dendrogram(z, leaf_rotation=90, leaf_font_size=6, labels=lead_reporter_pivot.index,
                     above_threshold_color="#AAAAAA", truncate_mode='lastp', p=15, show_contracted=True,
                     annotate_above=6, max_d=max_d)
    plt.savefig(f"{h_cluster_output_path}/{lead_reporter}_dendrogram.pdf", dpi=300)
    plt.close()
    clusters = fcluster(z, max_d, criterion='distance')
    # k = 2
    # fcluster(z, k, criterion='maxclust')
    cluster_df = pd.DataFrame({"track_index": lead_reporter_pivot.index, "cluster_id": clusters})
    lead_pivot_cluster = lead_reporter_pivot.merge(cluster_df, how="inner", left_index=True, right_on="track_index",
                                                   validate="1:1")
    lead_pivot_cluster[["site", "track", "cell__treatment"]] = lead_pivot_cluster["track_index"].str.split(pat="_",
                                                                                                           n=2,
                                                                                                           expand=True)
    lead_pivot_cluster[["cell", "treatment"]] = lead_pivot_cluster["cell__treatment"].str.split(pat="__", n=1,
                                                                                                expand=True)
    lead_melt_cluster = pd.melt(lead_pivot_cluster,
                                id_vars=["track_index", "cluster_id", "site", "track", "cell", "treatment",
                                         "cell__treatment"],
                                var_name="time", value_name=lead_reporter)
    lead_melt_cluster['treatment'] = pd.Categorical(lead_melt_cluster['treatment'], categories=treatment_order,
                                                    ordered=True)
    melt_cluster_merged = lead_melt_cluster.merge(dataframe[["track_index", "time", sec_reporter]], how="inner",
                                                  on=["track_index", "time"], validate="1:1")
    # Sort cell tracks based on 1) treatment, 2) cluster_id at t0, 3) most frequent cluster_id afterwards
    lead_melt_cluster.sort_values(by=["treatment", "cluster_id"], ascending=True, inplace=True)
    cluster_n = melt_cluster_merged["cluster_id"].nunique()
    if cluster_n < 5:
        col_n = cluster_n
    else:
        col_n = 5
    sns.set_context("paper")
    sns.relplot(data=melt_cluster_merged, x="time", y=lead_reporter, kind="line", col="cluster_id", hue="track_index",
                col_wrap=col_n, legend=None, ci=None, palette=lead_color)
    plt.savefig(f"{h_cluster_output_path}/{lead_reporter}_cluster_patterns.pdf", dpi=300)
    plt.close()
    sns.relplot(data=melt_cluster_merged, x="time", y=sec_reporter, kind="line", col="cluster_id",
                hue="track_index",
                col_wrap=col_n, legend=None, ci=None, palette=sec_color)
    plt.savefig(f"{h_cluster_output_path}/{lead_reporter}_clustered_{sec_reporter}_patterns.pdf", dpi=300)
    plt.close()
    sns.set_context("talk")
    sns.relplot(data=melt_cluster_merged, x="time", y=lead_reporter, kind="line", col="cluster_id", hue="treatment",
                col_wrap=col_n, legend="full", ci=None, palette=lead_color)
    plt.savefig(f"{h_cluster_output_path}/represented_{lead_reporter}_cluster_patterns.pdf", dpi=300)
    plt.close()
    sns.relplot(data=melt_cluster_merged, x="time", y=sec_reporter, kind="line", col="cluster_id", hue="treatment",
                col_wrap=col_n, legend="full", ci=None, palette=sec_color)
    plt.savefig(f"{h_cluster_output_path}/represented_{lead_reporter}_clustered_{sec_reporter}_patterns.pdf", dpi=300)
    plt.close()
    cluster_counts = (lead_melt_cluster.groupby(["treatment"], sort=False)['cluster_id']
                      .value_counts(normalize=True)
                      .rename('percentage')
                      .mul(100)
                      .reset_index())
    # sns.set(font_scale=0.6)
    # p = sns.barplot(data=cluster_counts, x="treatment", y="percentage", hue="cluster_id")
    # _ = plt.setp(p.get_xticklabels(), rotation=90)  # Rotate labels
    # p.legend(loc='upper center', bbox_to_anchor=(-0.5, 1.01, 2, 0.1), ncol=5, fontsize="small")
    treatment_n = cluster_counts["treatment"].nunique()
    if treatment_n < 5:
        col_n2 = treatment_n
    else:
        col_n2 = 5
    sns.catplot(data=cluster_counts, x="cluster_id", y="percentage", kind="bar", col="treatment", col_wrap=col_n2,
                legend="full")
    # ax = sns.countplot(data=erk_melt_cluster, x="cluster_id", hue="treatment")
    # ax.legend(loc='upper center', bbox_to_anchor=(-0.5, 1, 2, 0.1), ncol=4, fontsize="small")
    # g.set_xticklabels(rotation=90)
    plt.savefig(f"{h_cluster_output_path}/{lead_reporter}_cluster_population_per_treatment.pdf", dpi=300)
    plt.close()
    export_csv(lead_pivot_cluster, f"{lead_reporter}_hierarchical_cluster_df", False)
    return lead_pivot_cluster


def h_cluster_heatmap(dataframe, path_reporter):
    num_of_subplots = dataframe["cluster_id"].nunique()
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
    dataframe.sort_values(by="cluster_id", ascending=True, inplace=True)
    df_group = dataframe.groupby("cluster_id", sort=False)
    i = 1
    fig = plt.figure(figsize=(fig_width, fig_height), tight_layout=True)  # width x height
    for name, group in df_group:
        work_df = group.iloc[:, :len(dataframe.columns) - 7]
        current_palette = sns.color_palette("RdBu_r", n_colors=1000)
        cmap = ListedColormap(sns.color_palette(current_palette).as_hex())
        sns.set_context("paper", font_scale=2)  # rc={"font.size":2,"axes.labelsize":2})
        fig.add_subplot(n_rows, n_columns, i)  # row, column, position
        # sns.heatmap(dataframe, ax=ax, cmap="RdBu_r", center=0, vmin=-1, vmax=1, yticklabels=False)
        # title = name.replace(f"{cell_reporter}__", "")
        title = f"cluster_{name}"
        total_cell_number = len(work_df.index)
        ax = sns.heatmap(work_df, cmap=cmap, vmin=0, vmax=1, yticklabels=False)
        ax.set_title(f"{title}\nn={total_cell_number}", fontsize="small")
        i += 1
    fig.savefig(f"{h_cluster_output_path}/clustered_{path_reporter}_heatmap.pdf", dpi=300)
    plt.close(fig)
    return


########### Pearson correlation ########################################################################################
def overall_pearson_corr(dataframe):
    clover_pivot = dataframe.pivot(index="track_index", columns="time", values=clover_reporter)
    mscarlet_pivot = dataframe.pivot(index="track_index", columns="time", values=mscarlet_reporter)
    overall_corr = clover_pivot.corrwith(mscarlet_pivot, axis=1, method="pearson")
    overall_corr_df = pd.DataFrame(overall_corr, index=clover_pivot.index, columns=["overall_pearson_r"])
    overall_corr_df.reset_index(drop=False, inplace=True)
    overall_corr_df[["site", "track", "cell__treatment"]] = overall_corr_df["track_index"].str.split(pat="_", n=2,
                                                                                                     expand=True)
    overall_corr_df[["cell", "treatment"]] = overall_corr_df["cell__treatment"].str.split(pat="__", n=1, expand=True)
    overall_corr_df['treatment'] = pd.Categorical(overall_corr_df['treatment'], categories=treatment_order,
                                                  ordered=True)
    overall_corr_df.sort_values(by=["treatment"], ascending=True, inplace=True)
    export_csv(overall_corr_df, "overall_pearson_corr_df", False)
    sns.set_context("paper")
    sns.boxplot(data=overall_corr_df, x="treatment", y="overall_pearson_r")
    sns.swarmplot(data=overall_corr_df, x="treatment", y="overall_pearson_r", color=".25", size=3)
    plt.savefig(f"{pearson_corr_output_path}/overall_pearson_corr_vs_treatment.pdf", dpi=300)
    plt.close()
    # Merge overall corr with h_clusters of the interested reporter (clover or mScarlet)
    for r in reporter_list:
        corr_cluster_df = overall_corr_df.merge(pivot_cluster_df_all[r][["track_index", "cluster_id"]],
                                                on="track_index", how="inner",
                                                validate="1:1")
        sns.set_context("paper")
        sns.boxplot(data=corr_cluster_df, x="cluster_id", y="overall_pearson_r")
        sns.swarmplot(data=corr_cluster_df, x="cluster_id", y="overall_pearson_r", color=".25", size=3)
        plt.savefig(f"{pearson_corr_output_path}/overall_pearson_corr_vs_{r}_h_cluster.pdf", dpi=300)
        plt.close()
    return


def fractioned_pearson_corr(dataframe, win_size):
    # Set window size to compute moving window synchrony.
    # Sort by track_index and time for correct rolling
    dataframe.sort_values(by=["track_index", "time"], ascending=True, inplace=True)
    rolling_df = {}
    # Set rolling window size ~ 3h
    r_window_size = win_size
    df_groups = dataframe.groupby("track_index", sort=False)
    for name, group in df_groups:
        rolling_r = group[clover_reporter].rolling(window=r_window_size).corr(
            group[mscarlet_reporter])  # , center=True, min_periods=1
        rolling_df_part = pd.DataFrame(rolling_r, columns=["pearson_corr_r"])
        rolling_df_part["track_index"] = name
        rolling_df_part[["time", "treatment"]] = group[["time", "treatment"]]
        rolling_df[name] = rolling_df_part
    rolling_df_all = pd.concat(rolling_df, join="outer", ignore_index=True)
    rolling_df_all.dropna(axis=0, inplace=True)
    export_csv(rolling_df_all, f"winsize_{r_window_size}_pearson_df_all", False)
    treatment_n = rolling_df_all["treatment"].nunique()
    if treatment_n < 5:
        col_n = treatment_n
    else:
        col_n = 5
    sns.set_context("talk")
    sns.relplot(data=rolling_df_all, x="time", y="pearson_corr_r", col="treatment", col_wrap=col_n, kind="line", ci="sd", estimator="median")
    plt.savefig(f"{pearson_corr_output_path}/winsize_{r_window_size}_pearson_corr_lineplot.pdf", dpi=300)
    plt.close()
    pearson_r_pivot_df = rolling_df_all.pivot(index="track_index", columns="time", values="pearson_corr_r")
    # f, ax = plt.subplots(2, 1, figsize=(14, 6), sharex="all")
    # df.rolling(window=30, center=True).median().plot(ax=ax[0])
    # ax[0].set(xlabel='Frame', ylabel='Smiling Evidence')
    # rolling_r.plot(ax=ax[1])
    # ax[1].set(xlabel='Frame', ylabel='Pearson r')
    # plt.suptitle("Smiling data and rolling window correlation")
    z = linkage(pearson_r_pivot_df, "ward")
    c, coph_dists = cophenet(z, pdist(pearson_r_pivot_df))
    print(c)
    plt.figure(figsize=(25, 10), tight_layout=True)
    sns.set_context("paper")
    plt.xlabel("track_index")
    plt.ylabel("distance")
    fancy_dendrogram(z, leaf_rotation=90, leaf_font_size=6, labels=pearson_r_pivot_df.index,
                     above_threshold_color="#AAAAAA", truncate_mode='lastp', p=30, show_contracted=True,
                     annotate_above=6)
    plt.show()
    plt.close()
    max_d = float(input(f"Winsize_{r_window_size}_pearson_cluster cutoff: "))
    fancy_dendrogram(z, leaf_rotation=90, leaf_font_size=6, labels=pearson_r_pivot_df.index,
                     above_threshold_color="#AAAAAA", truncate_mode='lastp', p=15, show_contracted=True,
                     annotate_above=6, max_d=max_d)
    plt.savefig(f"{pearson_corr_output_path}/winsize_{r_window_size}_pearson_corr_dendrogram.pdf", dpi=300)
    plt.close()
    clusters = fcluster(z, max_d, criterion='distance')
    pearson_r_pivot_df["cluster_id"] = clusters
    pearson_r_pivot_df.reset_index(drop=False, inplace=True)
    pearson_cluster_melt_df = pd.melt(pearson_r_pivot_df,
                                      id_vars=["track_index", "cluster_id"],
                                      var_name="time", value_name="pearson_corr_r")
    merge_df_all = dataframe.merge(pearson_cluster_melt_df, how="left", on=["track_index", "time"], validate="1:1")
    merge_df_all.sort_values(by=["track_index", "time"], ascending=True, inplace=True)
    merge_df_all["cluster_id"].fillna(method="bfill", axis=0, inplace=True)
    merge_df_all['treatment'] = pd.Categorical(merge_df_all['treatment'], categories=treatment_order,
                                               ordered=True)
    merge_df_all.sort_values(by=["treatment", "cluster_id"], ascending=True, inplace=True)
    cluster_n = merge_df_all["treatment"].nunique()
    if cluster_n < 5:
        col_n2 = cluster_n
    else:
        col_n2 = 5
    sns.set_context("talk")
    sns.relplot(data=merge_df_all, x="time", y="pearson_corr_r", kind="line", col="cluster_id", col_wrap=col_n2,
                legend="full", ci="sd")
    plt.savefig(f"{pearson_corr_output_path}/represented_pearson_corr_clusters.pdf", dpi=300)
    plt.close()
    sns.set_context("paper")
    sns.relplot(data=merge_df_all, x="time", y="pearson_corr_r", kind="line", hue="track_index", palette="Blues",
                col="cluster_id", col_wrap=col_n2, legend=None, ci=None)
    plt.savefig(f"{pearson_corr_output_path}/pearson_corr_clusters.pdf", dpi=300)
    plt.close()
    sns.set_context("paper")
    cluster_counts = (merge_df_all.groupby(["treatment"], sort=False)['cluster_id']
                      .value_counts(normalize=True)
                      .rename('percentage')
                      .mul(100)
                      .reset_index())
    sns.catplot(data=cluster_counts, x="cluster_id", y="percentage", kind="bar", col="treatment", col_wrap=col_n)
    plt.savefig(f"{pearson_corr_output_path}/pearson_corr_cluster_population.pdf", dpi=300)
    plt.close()
    export_csv(merge_df_all, f"winsize{r_window_size}_pearson_corr_cluster_dataframe", False)
    return merge_df_all


def pearson_reporter_lineplot(dataframe):
    reporter_melt_df = pd.melt(dataframe,
                               id_vars=["track_index", "time", "site", "track", "cell__treatment", "time_range",
                                        f"{clover_reporter}_scaled", f"{mscarlet_reporter}_scaled", "cell", "treatment",
                                        "cluster_id", "pearson_corr_r", "time_scaled"], var_name="reporter",
                               value_name="reporter_activity")
    cluster_n = reporter_melt_df["treatment"].nunique()
    if cluster_n < 5:
        col_n = cluster_n
    else:
        col_n = 5
    sns.set_context("talk")
    sns.relplot(data=reporter_melt_df, x="time", y="reporter_activity", col="cluster_id", col_wrap=col_n, kind="line",
                hue="reporter", legend="full", ci=None, palette=["green", "red"], estimator="median")
    plt.savefig(f"{pearson_corr_output_path}/pearson_corr_cluster_reporter_dynamics.pdf", dpi=300)
    plt.close()
    return


def crosscorr(datax, datay, lag=0, wrap=False):
    """ Lag-N cross correlation.
    Shifted data filled with NaNs

    :param
    ----------
    lag : int, default 0
    datax, datay : pandas.Series objects of equal length
    Returns
    ----------
    crosscorr : float
    """
    if wrap:
        shiftedy = datay.shift(lag)
        shiftedy.iloc[:lag] = datay.iloc[-lag:].values
        return datax.corr(shiftedy)
    else:
        return datax.corr(datay.shift(lag))


def overall_crosscorr(dataframe, corr_direction, f):
    rows = []
    df_groups = dataframe.groupby("track_index", sort=False)
    for name, group in df_groups:
        d1 = group[clover_reporter]
        d2 = group[mscarlet_reporter]
        time_unit = 1
        frames = f
        rs = [crosscorr(d1, d2, lag) for lag in range(-int(time_unit * frames), int(time_unit * frames))]
        offset = ""
        if corr_direction == "pos":
            offset = np.ceil(len(rs) / 2) - np.argmax(rs)
        elif corr_direction == "neg":
            offset = np.ceil(len(rs) / 2) - np.argmin(rs)
        else:
            print("WRONG corr_direction input!")
        rows.append([name, offset])
        # Plot individual time lagged correlation lineplot
        # end = (time_unit * frames * 2)
        # increment = end / 6
        # f, ax = plt.subplots(figsize=(14, 3))
        # ax.plot(rs)
        # ax.axvline(np.ceil(len(rs) / 2), color='k', linestyle='--', label='Center')
        # ax.axvline(np.argmax(rs), color='r', linestyle='--', label='Peak synchrony')
        # ax.set(title=f'Offset = {offset} frames\nS1 leads <> S2 leads', ylim=[-1.1, 1.1], xlim=[0, end], xlabel='Offset',
        #        ylabel='Pearson r')
        # ax.set_xticks([0, increment, increment * 2, (increment * 3), (increment * 4), (increment * 5), end])
        # ax.set_xticklabels([-increment * 3, -increment * 2, -increment, 0, increment, increment * 2, increment * 3])
        # plt.legend()
        # plt.show()
    df = pd.DataFrame(rows, columns=["track_index", f"{clover_reporter}_lead_<_0_<_{mscarlet_reporter}_lead"])
    merged_df = dataframe.merge(df, on="track_index", how="inner", validate="m:1")
    merged_df['treatment'] = pd.Categorical(merged_df['treatment'], categories=treatment_order, ordered=True)
    merged_df.sort_values(by="treatment", inplace=True)
    offset_counts = (
        merged_df.groupby(["treatment"], sort=False)[f"{clover_reporter}_lead_<_0_<_{mscarlet_reporter}_lead"].value_counts(normalize=True).rename('percentage').mul(100).reset_index())
    if offset_counts["treatment"].nunique() < 5:
        col_wrap_value = offset_counts["treatment"].nunique()
    else:
        col_wrap_value = 5
    offset_counts[f"{clover_reporter}_lead_<_0_<_{mscarlet_reporter}_lead"] = offset_counts[
        f"{clover_reporter}_lead_<_0_<_{mscarlet_reporter}_lead"].astype(int)
    sns.set_context("talk")
    sns.catplot(data=offset_counts, x=f"{clover_reporter}_lead_<_0_<_{mscarlet_reporter}_lead", y="percentage",
                kind="bar",
                col="treatment", col_wrap=col_wrap_value)
    plt.savefig(f"{pearson_corr_output_path}/overall_time_lagged_corr_vs_treatment.pdf", dpi=300)
    export_csv(merged_df, "overall_time_lagged_correlation", False)
    plt.close()
    return


cluster_color_palette = "muted"
reporter_list = [clover_reporter, mscarlet_reporter]
unpivot_clover = unpivot_df(hmap_df_c, "clover")
unpivot_mscarlet = unpivot_df(hmap_df_m, "mscarlet")
unpivot_merge = merge_df(unpivot_clover, unpivot_mscarlet)
unpivot_gps = groups(unpivot_merge)
path_coordination_scatterplot(unpivot_merge, None, "raw")

################# create QUANDRANT plots and data ######################################################################
threshold_database = quadrant_lines_threshold(unpivot_merge)
clover_threshold = threshold_database[clover_reporter].median()
mscarlet_threshold = threshold_database[mscarlet_reporter].median()
quadrant_stats_time_data = quadrant_stats_df(unpivot_gps, "time")
quadrant_stats_time_range_results = quadrant_stats_df(unpivot_gps, "time_range")
export_csv(quadrant_stats_time_data, "quadrant_line_data", False)
export_csv(quadrant_stats_time_range_results, "quadrant_statistics", False)
quadrant_lineplot(quadrant_stats_time_data)
quad_path_coordination_scatterplot(unpivot_gps)
q_cluster_dataframe = quadrant_membership(unpivot_merge)
q_cluster_id_t0_dataframe = cell_cluster_id_t0_df(q_cluster_dataframe, "quadrant")
q_cluster_switch_pivot_dataframe, q_cluster_label = cluster_switch_melt_to_pivot_df(q_cluster_dataframe, "quadrant")
cluster_switch_frequency(q_cluster_switch_pivot_dataframe, "quadrant")
cluster_switch_trajectory_heatmap(q_cluster_switch_pivot_dataframe, "quadrant")

################# create K-mean cluster plots and data #################################################################
##### USE ALL CELLS #####
k_cluster_dataframe = clustering_on_df(unpivot_merge, "k_mean")
k_cluster_id_t0_dataframe = cell_cluster_id_t0_df(k_cluster_dataframe, "k_mean")
k_gp_cluster_id_t0_dataframe = groups(k_cluster_id_t0_dataframe)
path_coordination_scatterplot(k_cluster_id_t0_dataframe, cluster_color_palette, "k_mean")
cluster_population_lineplot(k_gp_cluster_id_t0_dataframe, "k_mean")
k_cluster_switch_pivot_dataframe, k_cluster_label = cluster_switch_melt_to_pivot_df(k_cluster_dataframe, "k_mean")
cluster_switch_frequency(k_cluster_switch_pivot_dataframe, "kmeans")
cluster_switch_trajectory_heatmap(k_cluster_switch_pivot_dataframe, "k_mean")

################# hierarchical clustering ##############################################################################
pivot_cluster_df_all = {}
for reporter in reporter_list:
    pivot_cluster = hierarchical_clustering(unpivot_merge, reporter)
    pivot_cluster_df_all[reporter] = pivot_cluster
    h_cluster_heatmap(pivot_cluster, reporter)

################# pearson correlation ##################################################################################
overall_pearson_corr(unpivot_merge)
pearson_corr_dataframe = fractioned_pearson_corr(unpivot_merge, 8)  # Argument 2) input value for the window size (number of frames)
pearson_reporter_lineplot(pearson_corr_dataframe)
overall_crosscorr(unpivot_merge, "neg", 30)  # Argument 2) "pos" or "neg" for second argument, Argument 3) input value for testing range of frameshift
