import logging
import glob
import os
import gc
import shutil
import math
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.colors import ListedColormap
import cell_reporter_list as cr_list
from Set_parameters import true_signal_highest
from Set_parameters import true_signal_lowest
from Set_parameters import signal_max
from Set_parameters import signal_min
from Set_parameters import treatment_order
from Set_parameters import source_path

# ---------------------------------------------- SET UP / OPTIONS-------------------------------------------------------
# size of table to display
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

# input , time interval, reporter color # & reporter translocation
rc1_input = ""
rc2_input = ""
rc1 = ""
rc2 = ""
rc_list = ""
reporter_translocation_1 = ""
reporter_translocation_2 = ""
reporter_number = input("number of reporters: ")
ch_counter = len(glob.glob1(f"{source_path}/Batch1", 'Cell_Cytoplasm_Intensity_Mean_*'))
reporter_channel_list = ["clover", "mscarlet"]
print(
    """reporter_channel_list:
    0 = clover
    1 = mscarlet
    """
)
if int(reporter_number) > 2:
    print("WARNING: more reporters than expected")
    input("please RE-RUN")
# case of 1 reporter
elif int(reporter_number) == 1:
    if int(reporter_number) > ch_counter - 1:
        print("WARNING: more reporters than expected")
        input("please RE-RUN")
    elif int(reporter_number) < ch_counter - 1:
        print("WARNING: expect more reporters")
        input("please RE-RUN")
    else:
        rc1_input = input("reporter 1 color: ")
    if int(rc1_input) > 1:
        print("wrong color input")
        input("please RE-RUN")
    else:
        rc1 = reporter_channel_list[int(rc1_input)]
        rc_list = [rc1]
        reporter_translocation_1 = input(
            "translocation of reporter 1 when activates (0. to cytoplasm   1. to nucleus): ")
        if int(reporter_translocation_1) > 1:
            print("WARNING: wrong translocation input")
            input("please RE-RUN")
# case of 2 reporters
elif int(reporter_number) == 2:
    if int(reporter_number) > ch_counter - 3:
        print("WARNING: more reporters than expected")
        input("please RE-RUN")
    elif int(reporter_number) < ch_counter - 3:
        print("WARNING: expect more reporters")
        input("please RE-RUN")
    else:
        rc1_input = input("reporter 1 color: ")
    if int(rc1_input) > 1:
        print("wrong color input")
        input("please RE-RUN")
    else:
        rc1 = reporter_channel_list[int(rc1_input)]
        reporter_translocation_1 = input(
            "translocation of reporter 1 when activates (0. to cytoplasm   1. to nucleus): ")
        if int(reporter_translocation_1) > 1:
            print("WARNING: wrong translocation input")
            input("please RE-RUN")
        rc2_input = input("reporter 2 color: ")
        if int(rc2_input) > 1:
            print("wrong color input")
            input("please RE-RUN")
        elif rc2_input == rc1_input and rc1_input != "2":
            print("Warning: two reporters with the same color assigned")
            input("please RE-RUN")
        elif rc1_input == "2" and rc2_input == "2":
            print("Warning: no reporter assigned")
            input("please RE-RUN")
        else:
            rc2 = reporter_channel_list[int(rc2_input)]
            rc_list = [rc1, rc2]
            reporter_translocation_2 = input(
                "translocation of reporter 2 when activates (0. to cytoplasm   1. to nucleus): ")
            if int(reporter_translocation_2) > 1:
                print("WARNING: wrong translocation input")
                input("please RE-RUN")

time_interval = int(input("""
time between frames (min): """))

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

cell_reporter = cell_reporter_list[int(cell_reporter_input)]

# Create output folder and subfolders
output_path = f"{source_path}/output_abs"
qc_directory = "quality_control"
qc_path = os.path.join(output_path, qc_directory)
os.makedirs(output_path, exist_ok=True)
os.makedirs(qc_path, exist_ok=True)
logging.basicConfig(filename=f'{output_path}/log_file.log', level=logging.INFO)

# data processing for all reporters
numerator = ''
denominator = ''
ch = ""
color = ""
palette = ""
folder_list = ["heatmap", "boxplot", "lineplot", "catplot"]
reporter_translocation = ""
for rc in rc_list:
    rc_directory = rc
    rc_path = os.path.join(output_path, rc_directory)
    os.makedirs(rc_path, exist_ok=True)
    for folder in folder_list:
        folder_directory = folder
        folder_path = os.path.join(rc_path, folder_directory)
        os.makedirs(folder_path, exist_ok=True)
    sorted_by_div_directory = "sorted_by_div"
    sorted_by_div_path = os.path.join(f"{output_path}/{rc}/heatmap", sorted_by_div_directory)
    os.makedirs(sorted_by_div_path, exist_ok=True)

    # link each reporter to its translocation when activates
    # add comments here
    if rc == rc1:
        reporter_translocation = reporter_translocation_1
    elif rc == rc2:
        reporter_translocation = reporter_translocation_2

    # assign numerator and denominator identity for c/n or n/c ratio calculation
    if reporter_translocation == "0":
        numerator = "cell_cytoplasm_intensity_mean"
        denominator = "nucleus_intensity_mean"
    elif reporter_translocation == "1":
        numerator = "nucleus_intensity_mean"
        denominator = "cell_cytoplasm_intensity_mean"
    else:
        print("wrong input on reporter translocation")
        break

    # assign ch#, color and palette based on reporter number and identity
    if reporter_number == "1":
        ch = 2
        if rc == "clover":
            color = "green"
            palette = "Greens"
        elif rc == "mscarlet":
            color = "red"
            palette = "Reds"
    elif reporter_number == "2":
        if rc == "mscarlet":
            ch = 2
            color = "red"
            palette = "Reds"
        elif rc == "clover":
            ch = 3
            color = "green"
            palette = "Greens"
    else:
        print("wrong reporter color input")
        break

    # ---------------------------------------------GET ALL DATA FILES---------------------------------------------------

    # combine files from different batches
    def file_combine(folder2, source, filename):
        batch_counter = len(glob.glob1(source_path, 'Batch*'))
        extension = '.csv'
        count = 1
        while count < batch_counter + 1:
            shutil.copy(f"{source_path}/Batch{count}/{source}{extension}",
                        f"{output_path}/quality_control/{filename}_{count}{extension}")
            count += 1
        all_filenames = [i for i in glob.glob(f"{output_path}/quality_control/{filename}*{extension}")]

        # combine all files in the list
        combined_csv = pd.concat([pd.read_csv(f, header=2) for f in all_filenames])
        # export to csv
        combined_csv.to_csv(f"{output_path}/{folder2}/combined_{filename}{extension}", index=False,
                            encoding="utf-8")
        return f"{output_path}/{folder2}/combined_{filename}{extension}"


    # -----------------------------------------------DATA FORMATTING--------------------------------------------------------

    # cleanup column header
    def clean_header(dataframe):
        dataframe2 = dataframe
        dataframe2.columns = dataframe2.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('(',
                                                                                                          '').str.replace(
            ')', '')
        return dataframe2


    # replace values in column "cellid" to column"id" & rename "cellid" as "id"
    # for nucleus dataframe only
    def cellid_to_id(dataframe):
        dataframe2 = dataframe
        dataframe2["id"] = dataframe2["cellid"].astype(str)
        dataframe2.drop(["cellid"], axis=1, inplace=True)
        dataframe2["trackid"] = dataframe2['id'].astype(str)
        return dataframe2


    # add "site" & "cn_lookup" columns
    def add_site_cnlookup(dataframe, key_column):
        ### ims. file coverted using imaris arena
        dataframe[["temp1", "temp2"]] = dataframe.original_image_id.str.split("_s", expand=True)
        dataframe[["site", "temp3"]] = dataframe.temp2.str.split(".ome", expand=True)
        ### split into .ome and use imaris converter
        # dataframe[["temp1", "temp2"]] = dataframe.original_image_id.str.split("_s", expand=True)
        # dataframe[["site", "temp3"]] = dataframe.temp2.str.split(".ome", expand=True)
        ### ims. file converted using imaris converter
        # dataframe[["temp1", "temp2"]] = dataframe.original_image_id.str.split("scene-", expand=True)
        # dataframe[["site", "temp3"]] = dataframe.temp2.str.split("-p", expand=True)
        # add line here to search the column names for "unnamed", then remove just this column
        # because IMARIS prints out column with "unnamed" in prefix and this may change in future
        dataframe = dataframe.loc[:, [key_column, "category", "birth_[s]", "death_[s]", "trackid", "id", "site"]]
        col = ["id", "trackid", "site"]
        dataframe[col] = dataframe[col].astype(str)
        dataframe["cn_lookup"] = (dataframe["site"] + "_" + dataframe["id"]).astype(str)
        return dataframe


    # merge cytoplasm and nucleus dataframe based on cn_lookup
    def cn_merge(dataframe_c, dataframe_n):
        return (pd.merge(dataframe_c, dataframe_n[["nucleus_intensity_mean", "cn_lookup"]],
                         on="cn_lookup", how="inner", suffixes=("_c", "_n")).astype({"trackid": "str"}))


    # add "cell__treatment" column to dataframe
    def df_pmap_merge(dataframe, p_map):
        dataframe_merge = pd.merge(dataframe, p_map[["site", "cell__treatment"]], on="site", how="left")
        return dataframe_merge


    ## DATA CLEAN UP AND PROCESSING
    ## NEED WRITTEN VERSION FOR METHODS/PROTOCOL IN PAPERS

    # 1) delete rows with cytoplasm & nucleus mean intensity <=700
    # 2) add "cn_ratio" ,
    # 3) delete rows with "cn_ratio" < 1 percentile and > 95 percentile
    def true_signal(dataframe):
        dataframe2 = dataframe.copy()
        index_names = dataframe2[
            (dataframe2["cell_cytoplasm_intensity_mean"] < 550) & (
                    dataframe2["nucleus_intensity_mean"] < 550)].index
        dataframe2.drop(index_names, axis=0, inplace=True)
        dataframe2["cn_ratio"] = (dataframe2[numerator]) / (dataframe2[denominator])
        index_name2 = dataframe2[(dataframe2["cn_ratio"] > dataframe2["cn_ratio"].quantile(true_signal_highest))].index
        index_name3 = dataframe2[(dataframe2["cn_ratio"] < dataframe2["cn_ratio"].quantile(true_signal_lowest))].index
        dataframe2.drop(index_name2, axis=0, inplace=True)
        dataframe2.drop(index_name3, axis=0, inplace=True)
        return dataframe2


    # add "time", "track_index" & "track_time_index" columns
    # delete rows with duplicated "track_time_index" values -- keep only one cell from 2 daughter cells

    def add_time(dataframe):
        dataframe2 = dataframe.copy()
        dataframe2["time_h"] = (dataframe2['birth_[s]'] * time_interval / 60).round(1)
        dataframe2["track_index"] = dataframe2["site"] + "_" + dataframe2["trackid"] + "_" + dataframe2[
            "cell__treatment"]
        dataframe2["track_time_index"] = dataframe2['track_index'] + "_" + dataframe2["time_h"].astype(str)
        dataframe2.drop_duplicates("track_time_index", inplace=True)
        return dataframe2


    # quality check -- display row with "nan" value
    def display_na(dataframe):
        dataframe_na = dataframe[dataframe.isna().any(axis=1)]
        return dataframe_na


    # -----------------------------------PREP DATA FOR HEATMAP VISUALIZATION------------------------------------------------

    # 1) generate matrix dataframe for reporter and n_size
    # 2) fill maximum ~10% missing values using linear method & then delete rows which still contain nan
    # 3) "cn_ratio" is normalized by abs. *max (98 percentile of all values)
    # 4) add "site" & "cell__treatment" column
    def df_to_hmap_nuc_only(dataframe, values, data_identity):
        dataframe_matrix = dataframe.pivot(index="track_index", columns="time_h", values=values)
        fill_limits = int(len(dataframe_matrix.columns) * 0.1)
        before_drop = dataframe_matrix.shape[0]
        dataframe_matrix.dropna(subset=[0.0], axis=0, inplace=True)
        dataframe_matrix.interpolate(method="linear", axis=1, limit=fill_limits, limit_direction="backward", inplace=True)
        dataframe_matrix.dropna(axis=0, inplace=True)
        # check in case  we miss many per site ---> check by full plate - 240 sites
        after_drop = dataframe_matrix.shape[0]
        num_of_drops = before_drop - after_drop
        percentage_of_drops = round(((num_of_drops / before_drop) * 100), 2)
        logging.info(f"{rc}_{data_identity}_number of drops: {num_of_drops}")
        logging.info(f"{rc}_{data_identity}_percentage of drops: {percentage_of_drops}%")
        # dataframe_matrix = dataframe_matrix.subtract(dataframe_matrix[0.0], axis=0)
        # dataframe_matrix = dataframe_matrix.div(dataframe_matrix[0.0], axis=0)
        # we are using MAX value per plate here to scale to
        maximum = dataframe_matrix.quantile(signal_max, axis=1).quantile(signal_max)
        minimum = dataframe_matrix.quantile(signal_min, axis=1).quantile(signal_min)
        dataframe_matrix = (dataframe_matrix - minimum) / (maximum - minimum)
        logging.info(f"{rc}_{data_identity}_maximum: {maximum}")
        logging.info(f"{rc}_{data_identity}_minimum: {minimum}\n")
        if data_identity == "reporter":
            dataframe_matrix["temp"] = dataframe_matrix.index
            dataframe_matrix[["site", "temp2", "cell__treatment"]] = dataframe_matrix["temp"].str.split(pat="_", n=2,
                                                                                                        expand=True)
            dataframe_matrix["site"] = dataframe_matrix["site"].astype(str)
            dataframe_matrix.drop(["temp", "temp2"], axis=1, inplace=True)
        elif data_identity == "n_size":
            dataframe_matrix = pd.merge(left=hmap_rc_sorted["cell__treatment"], right=dataframe_matrix, left_index=True,
                                        right_index=True, how="left")
        return dataframe_matrix


    # add columns to sort for reporter activity: "median"
    # threshold calculated based on "untreated" median
    # sort based on "above thres" + "median" or "median" only
    def hmap_sort(dataframe):
        dataframe2 = dataframe.copy()
        # ref_group = dataframe2.groupby("cell__treatment").get_group(f"{cell_reporter}__untreated")
        # threshold = ref_group.median(axis=1, skipna=True).median(axis=0, skipna=True)
        dataframe_num = dataframe2.iloc[:, 0:dataframe2.shape[1] - 2]
        # dataframe2["above_thres"] = dataframe_num[dataframe_num > threshold].count(axis=1)
        dataframe2["median"] = dataframe_num.median(skipna=True, axis=1)
        dataframe2[["cell", "treatment"]] = dataframe2["cell__treatment"].str.split(pat="__", n=1, expand=True)
        dataframe2['treatment'] = pd.Categorical(dataframe2['treatment'], categories=treatment_order, ordered=True)
        dataframe2.sort_values(by=["treatment", "median"], ascending=[True, False], inplace=True)
        dataframe2.drop(["cell", "treatment"], axis=1, inplace=True)
        # dataframe2.sort_values(by=["above_thres", "median"], ascending=False, inplace=True)
        return dataframe2


    # divide dataframe into groups based on "cell__treatment" info
    def groups(dataframe):
        df_groups = dataframe.groupby("cell__treatment", sort=False)
        return df_groups


    # export divided dataframe in groups
    def export_gps(df_groups, channel, data_identifier):
        for name, group in df_groups:
            group.to_csv(f"{output_path}/{channel}/{name}_{data_identifier}.csv", index=True)
        return


    # export complete hmap_dataframe
    def export_csv(dataframe, channel):
        dataframe.to_csv(f"{output_path}/{channel}/{channel}_all_cell.csv", index=True)
        return


    # export cell count summary over time - gives us # cells/time frame in total
    def export_count(df_count, name):
        df_count.to_csv(f"{output_path}/quality_control/{name}.csv", index=False)
        return


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


    # generate boxplot for reporter activity
    def reporter_boxplot(df_groups, channel, colors):
        for name, group in df_groups:
            fig = plt.figure(figsize=(25, 20))
            dataframe = group.iloc[:, 0:group.shape[1] - 4]
            # to show outliers:
            # flierprops = dict(markerfacecolor='0.5', markersize=4, linestyle='none')
            # ax = sns.boxplot(data=dataframe, color=color, whis=[5, 95], flierprops=flierprops)
            ax = sns.boxplot(data=dataframe, color=colors, whis=[5, 95], showfliers=False, linewidth=2.5)
            ax.spines['left'].set_linewidth(3)
            ax.spines['bottom'].set_linewidth(3)
            sns.despine()
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
            ax.set(ylim=(0, 1))
            sns.set_context("poster", rc={"font.size": 40, "axes.labelsize": 10})
            ax.xaxis.set_tick_params(width=3)
            ax.yaxis.set_tick_params(width=3)
            title = name.replace(f"{cell_reporter}__", "")
            total_cell_number = dataframe.shape[0]
            ax.set_title(f"{title}\nn={total_cell_number}", fontsize=80)
            fig.savefig(f"{output_path}/{channel}/boxplot/{title}.png")
            plt.close(fig)
        return


    # generate lineplot for reporter activity
    def reporter_lineplot(df_groups, channel):
        current_palette = sns.color_palette("Paired")
        cmap = ListedColormap(sns.color_palette(current_palette).as_hex())
        i = 2
        for name, group in df_groups:
            fig = plt.figure(figsize=(20, 20))
            sns.set_context("poster", rc={"font.size": 40, "axes.labelsize": 40})
            for name2, everything in df_groups:
                dataframe = everything.iloc[:, 0:everything.shape[1] - 3].mean(axis=0, skipna=True)
                sns.lineplot(data=dataframe, color="grey", lw=1)
                sns.despine()
            single_df = group.iloc[:, 0:group.shape[1] - 3].mean(axis=0, skipna=True)
            sns.set_context("talk")
            ax = sns.lineplot(data=single_df, color=cmap(i), lw=12)
            ax.spines['left'].set_linewidth(3)
            ax.spines['bottom'].set_linewidth(3)
            sns.despine()
            # ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
            ax.set(xlabel='time_h', ylabel='mean_reporter_activity')
            ax.xaxis.set_tick_params(width=3)
            ax.yaxis.set_tick_params(width=3)
            title = name.replace(f"{cell_reporter}__", "")
            total_cell_number = group.shape[0]
            ax.set_title(f"{title}\nn={total_cell_number}", fontsize=80)
            # ax.set(xlim=(0, 17))
            fig.savefig(f"{output_path}/{channel}/lineplot/{title}.png")
            plt.close(fig)
            i += 1
        return


    # execute all the functions above
    rc_c = file_combine(rc, f'Cell_Cytoplasm_Intensity_Mean_Ch={ch}_Img=1', f'cyto_{rc}')
    rc_n = file_combine(rc, f'Nucleus_Intensity_Mean_Ch={ch}_Img=1', f'nuc_{rc}')
    df_rc_c = pd.read_csv(rc_c)
    df_rc_n = pd.read_csv(rc_n)
    df_list = [df_rc_c, df_rc_n]
    df_list = [df.pipe(clean_header) for df in df_list]
    cellid_to_id(df_rc_n)
    df_rc_c = add_site_cnlookup(df_rc_c, "cell_cytoplasm_intensity_mean")
    df_rc_n = add_site_cnlookup(df_rc_n, "nucleus_intensity_mean")
    df_rc = cn_merge(df_rc_c, df_rc_n)
    plate_map = pd.read_csv(f"{source_path}/plate_map.csv").astype(str)
    df_rc_map_merged = df_pmap_merge(df_rc, plate_map)
    df_rc_true = true_signal(df_rc_map_merged)
    df_rc_time = add_time(df_rc_true)
    df_rc_na = display_na(df_rc_time)
    hmap_rc = df_to_hmap_nuc_only(df_rc_time, "cn_ratio", "reporter")
    hmap_rc_sorted = hmap_sort(hmap_rc)
    hmap_rc_groups = groups(hmap_rc_sorted)
    export_gps(hmap_rc_groups, rc, "reporter")
    df_rc_time.to_csv(f"{output_path}/{rc}/df_{rc}_time.csv", index=False, encoding="utf-8")
    df_rc_na.to_csv(f"{output_path}/{rc}/df_{rc}_na.csv", index=False)
    export_count(hmap_rc_groups.count(), f"hmap_{rc}_groups_count")
    heatmap(hmap_rc_groups, rc, df_open=0, df_end=4, colors="RdBu_r", data_identifier="reporter", minimum=0, maximum=1,
            destination="heatmap")
    reporter_boxplot(hmap_rc_groups, rc, color)
    reporter_lineplot(hmap_rc_groups, rc)
    export_csv(hmap_rc_sorted, rc)

    # quick check summary
    # logging.info("Statements here") -? top of script intialize logger file
    # logging.error"Statement here that only shows up for warning")
    logging.info(f"""
                PARAMETER SETTING:
                true signal (max, min) ------------ {rc}: {true_signal_highest}, {true_signal_lowest}
                signal max & min for color scale -- {rc}: {signal_max}, {signal_min}
                
                KEY DATAFRAME DIMENSIONS:
                cyto + nuc ------------------------ {rc}: {df_rc.shape}
                add col treatment ----------------- {rc}: {df_rc_map_merged.shape}
                true_signal ----------------------- {rc}: {df_rc_time.shape}
                hmap ------------------------------ {rc}: {hmap_rc.shape}
                sorted_hmap ----------------------- {rc}: {hmap_rc_sorted.shape}
                """)

    # ------------------------------------for cell division (cell division time & nucleus size)-------------------------
    # combine cell division files from different batches
    cell_div = file_combine(rc, 'Cell_Time_Since_First_Division', 'cell_division')
    # combine nucleus size files from different batches
    cell_nuc_size = file_combine(rc, "Nucleus_Number_of_Voxels", "cell_nuc_size")
    df_cell_div = pd.read_csv(cell_div)
    df_cell_n_size = pd.read_csv(cell_nuc_size)
    data_list = [df_cell_div, df_cell_n_size]
    data_list = [df.pipe(clean_header) for df in data_list]
    cellid_to_id(df_cell_n_size)
    df_cell_div = add_site_cnlookup(df_cell_div, "cell_time_since_first_division")
    df_cell_n_size = add_site_cnlookup(df_cell_n_size, "nucleus_number_of_voxels")

    # add "trackid" column
    def n_size_trackid(dataframe_n_size, dataframe_trackid):
        dataframe_n_size_2 = dataframe_n_size.copy()
        dataframe_n_size_2.drop(["trackid"], axis=1, inplace=True)
        dataframe2 = pd.merge(dataframe_n_size_2, dataframe_trackid[["cn_lookup", "trackid"]], on="cn_lookup",
                              how="outer")
        return dataframe2


    # case of "cell division time", keep only rows of first time since division
    def keep_first_time_div(dataframe):
        dataframe2 = dataframe.copy()
        indexnames = dataframe2[(dataframe2["cell_time_since_first_division"] != 0)].index
        dataframe2.drop(indexnames, axis=0, inplace=True)
        return dataframe2


    # create matrix dataframe for cell division time
    def df_to_hmap_div(hmap_sorted, df_div, values):
        div_matrix = df_div.pivot(index="track_index", columns="time_h", values=values)
        hmap_div_merge = pd.merge(left=hmap_sorted["cell__treatment"], right=div_matrix, left_index=True,
                                  right_index=True, how="left")
        hmap_div_merge.replace(to_replace="Cell", value=1, inplace=True)
        hmap_div_merge.fillna(value=0, inplace=True)
        hmap_div_merge[["cell", "treatment"]] = hmap_div_merge["cell__treatment"].str.split(pat="__", n=1, expand=True)
        hmap_div_merge['treatment'] = pd.Categorical(hmap_div_merge['treatment'], categories=treatment_order,
                                                     ordered=True)
        hmap_div_merge.sort_values(by=["treatment"], ascending=True, inplace=True)
        hmap_div_merge.drop(["cell", "treatment"], axis=1, inplace=True)
        return hmap_div_merge


    df_cell_div_map_merged = df_pmap_merge(df_cell_div, plate_map)
    df_cell_n_size_trackid = n_size_trackid(df_cell_n_size, df_rc)
    df_cell_n_size_map_merged = df_pmap_merge(df_cell_n_size_trackid, plate_map)
    df_cell_n_size_time = add_time(df_cell_n_size_map_merged)
    df_cell_div_first = keep_first_time_div(df_cell_div_map_merged)
    df_cell_div_time = add_time(df_cell_div_first)
    hmap_rc_div = df_to_hmap_div(hmap_rc_sorted, df_cell_div_time, "category")
    hmap_rc_n_size = df_to_hmap_nuc_only(df_cell_n_size_time, "nucleus_number_of_voxels", "n_size")
    hmap_rc_div_groups = groups(hmap_rc_div)
    hmap_rc_n_size_groups = groups(hmap_rc_n_size)
    heatmap(hmap_rc_div_groups, rc, colors=palette, df_open=1, df_end=0, minimum=0, maximum=1, data_identifier="div",
            destination="heatmap")
    heatmap(hmap_rc_n_size_groups, rc, colors=palette, df_open=1, df_end=0, minimum=0, maximum=1,
            data_identifier="n_size", destination="heatmap")
    export_gps(hmap_rc_div_groups, rc, "div")
    export_gps(hmap_rc_n_size_groups, rc, "n_size")

    # ---------------------------------------create heatmap based on cell division time---------------------------------
    # generate matrix dataframe for cell division time, only includes divided cells
    def hmap_div_sort_by_div(df_div, values):
        div_matrix = df_div.pivot(index="track_index", columns="time_h", values=values)
        div_matrix.replace(to_replace="Cell", value=1, inplace=True)
        div_matrix.fillna(value=0, inplace=True)
        # return column header index at the moment of cell division
        div_matrix["div"] = (div_matrix == 1).idxmax(axis=1)
        return div_matrix


    # create reporter matrix dataframe using cell division time as key
    # then sort by division time from late to early
    # cells without division will sit at the bottom
    def hmap_reporter_sort_by_div(hmap_div_sorted, hmap_sorted, sort_side):
        hmap_div_merge = pd.merge(left=hmap_div_sorted["div"], right=hmap_sorted, left_index=True,
                                  right_index=True, how=sort_side)
        # hmap_div_merge.sort_values(by=["div", "above_thres", "median"], ascending=False, inplace=True)
        hmap_div_merge[["cell", "treatment"]] = hmap_div_merge["cell__treatment"].str.split(pat="__", n=1, expand=True)
        hmap_div_merge['treatment'] = pd.Categorical(hmap_div_merge['treatment'], categories=treatment_order,
                                                     ordered=True)
        hmap_div_merge.sort_values(by=["treatment", "div", "median"], ascending=[True, False, False], inplace=True)
        hmap_div_merge.drop(["cell", "treatment"], axis=1, inplace=True)
        return hmap_div_merge


    # generate matrix for cell division time, sort by cell division time, includes non-divided cells
    def hmap_div_all(hmap_reporter_sort, hmap_div_sort):
        hmap_div_merge = pd.merge(left=hmap_reporter_sort["cell__treatment"], right=hmap_div_sort, left_index=True,
                                  right_index=True, how="left")
        hmap_div_merge.replace(to_replace="Cell", value=1, inplace=True)
        hmap_div_merge.fillna(value=0, inplace=True)
        return hmap_div_merge


    def pre_post_div(dataframe):
        dataframe.dropna(subset=[0.0], axis=0, inplace=True)
        dataframe = dataframe.iloc[:, :-1]
        dataframe["pre_div_median"] = ""
        dataframe["pos_div_median"] = ""
        # iterate through rows using the numerical order of rows as index
        for i in range(len(dataframe)):
            # set up split time for each row based on value in "div"
            # using get_loc to get the position of column header same as "div" on that row
            split_pt = dataframe.columns.get_loc(dataframe.iloc[i, 0])
            # split from 0.0 to col before split_pt & get median
            pre_div = dataframe.iloc[i, 1:split_pt].median(skipna=True)
            # split from col after split_pt to last time frame & get median
            pos_div = dataframe.iloc[i, split_pt + 1:-4].median(skipna=True)
            # assign pre-  & post-div median to each row when iterate
            dataframe.iloc[i, -2] = pre_div
            dataframe.iloc[i, -1] = pos_div
        dataframe[["pre_div_median", "pos_div_median"]] = dataframe[["pre_div_median", "pos_div_median"]].astype(float)
        dataframe["pos_vs_pre"] = dataframe["pos_div_median"] - dataframe["pre_div_median"]
        return dataframe


    def reporter_catplot(dataframe, channel):
        # for name, group in df_groups:
        dataframe[["cell", "treatment"]] = dataframe["cell__treatment"].str.split("__", expand=True)
        # fig = plt.figure(figsize=(20, 20))
        fig, ax = plt.subplots(figsize=(18, 18))  # width x height
        # current_palette = sns.color_palette("hus1")
        # cmap = ListedColormap(sns.color_palette(current_palette).as_hex())
        # sns.catplot(ax=ax, y="pos_vs_pre", x="treatment", kind="violin", inner="quartile", data=dataframe, scale="width")
        sns.catplot(ax=ax, kind="swarm", data=dataframe, x="treatment", hue="treatment", y="pos_vs_pre",
                    legend_out=True)
        ax.spines['left'].set_linewidth(3)
        ax.spines['bottom'].set_linewidth(3)
        ax.spines['top'].set_linewidth(0)
        ax.spines['right'].set_linewidth(0)
        # legend_labels = dataframe["treatment"].drop_duplicates().tolist()
        # ax.legend(labels=legend_labels)
        # ax.set(xlabel="Treatments", ylabel="pos-pre: cell division")
        # ax.set_xticklabels(ax.get_xticklabels(), rotation=60)
        # just hide x-axis tick labels but not x-axis title
        ax.get_xaxis().set_ticks([])
        ax.set(ylim=(-1, 1))
        sns.set_context("poster", font_scale=3)
        # rc={"font.size": 40, "axes.labelsize": 40})
        ax.xaxis.set_tick_params(width=3)
        ax.yaxis.set_tick_params(width=3)
        # total_cell_number = dataframe.shape[0]
        # ax.fig.suptitle(f"n={total_cell_number}", fontsize=60)
        fig.savefig(f"{output_path}/{channel}/catplot/{channel}.png")
        plt.close(fig)
        return


    # export the reporter median of pos-pre cell division for GraphPads
    def pre_pos_div_summary(dataframe, channel):
        dataframe[["cell", "treatment"]] = dataframe["cell__treatment"].str.split("__", expand=True)
        dataframe["index"] = dataframe.index
        dataframe2 = dataframe.pivot(index="index", columns="treatment", values="pos_vs_pre")
        dataframe2.to_csv(f"{output_path}/{channel}/pos_vs_pre_div_summary.csv", index=False)
        return dataframe2


    hmap_cell_div_sorted = hmap_div_sort_by_div(df_cell_div_time, "category")
    hmap_rc_sort_by_div = hmap_reporter_sort_by_div(hmap_cell_div_sorted, hmap_rc_sorted, sort_side="right")
    hmap_rc_div_sorted_groups = groups(hmap_rc_sort_by_div)
    heatmap(hmap_rc_div_sorted_groups, rc, df_open=1, df_end=4, colors="RdBu_r", data_identifier="reporter", minimum=0,
            maximum=1, destination="heatmap//sorted_by_div")
    hmap_div_all_sorted = hmap_div_all(hmap_rc_sort_by_div, hmap_cell_div_sorted)
    hmap_div_all_sorted_groups = groups(hmap_div_all_sorted)
    heatmap(hmap_div_all_sorted_groups, rc, colors=palette, df_open=1, df_end=1, minimum=0, maximum=1,
            data_identifier="div", destination="heatmap//sorted_by_div")
    export_gps(hmap_rc_div_sorted_groups, rc, "reporter_sbd")
    hmap_rc_div_only = hmap_reporter_sort_by_div(hmap_cell_div_sorted, hmap_rc_sorted, sort_side="left")
    hmap_rc_pre_post_div = pre_post_div(hmap_rc_div_only)
    hmap_rc_pre_post_div_groups = groups(hmap_rc_pre_post_div)
    export_gps(hmap_rc_pre_post_div_groups, rc, "pre_pos_div")
    reporter_catplot(hmap_rc_pre_post_div, rc)
    pre_pos_div_summary = pre_pos_div_summary(hmap_rc_pre_post_div, rc)
    gc.collect()
