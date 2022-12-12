from openbabel import pybel
import pandas as pd
import matplotlib.pyplot as plt


def read_file(format, filename):
    pybel.ob.obErrorLog.SetOutputLevel(0)
    mols = list(pybel.readfile(format, filename))
    return (mols)


def extract_metadata(mol_list, group_str):
    data = []
    mol_name = []
    smile = []
    GROUP_TYPE = []
    n_at = []
    for line in mol_list:
        data.append(line.data)
        mol_name.append(line.data["NAME"])
        GROUP_TYPE.append(line.data[group_str])
        smile.append(line.data["SMILES"])
        n_at.append(line.OBMol.NumAtoms())
    df = pd.DataFrame([mol_name, n_at, GROUP_TYPE, smile])
    pd.set_option('display.colheader_justify', 'center') 
    return (df)


def formating_dataframe(df, set_idx=False):
    # False for default index
    df = df.transpose()
    df.columns = ["Name", "N_at", "GROUP_TYPE", "Smile"]
    df = df.sort_values(by=['N_at'], ascending=False)
    if set_idx is True:
        df = df.set_index("GROUP_TYPE")
    return (df)


# def arrange_dataframe_pivot(df):
#     pivot = df.pivot(index="Name", columns="GROUP_TYPE", values="N_at")
#     pivot.index.name = "Name"
#     pivot.columns.name = "GROUP_TYPE"
#     return (pivot)


def clustering_and_statistics(df, GROUP_TYPE, stats=False):
    grouped = df.groupby(GROUP_TYPE)
    df_col = grouped[["N_at"]]
    size_df = grouped.size().to_frame(name="N_item/grp")
    max_at = df_col.max()
    df_merge = size_df.merge(max_at, how='inner', on='GROUP_TYPE')
    df_merge = df_merge.sort_values(by=['N_item/grp'], ascending=[False])
    df_merge.rename(columns = {'N_at':'Max_n_at'}, inplace = True)
    max_at.rename(columns = {'N_at':'Max_n_at'}, inplace = True)
    print("Number of groups found --> ", grouped.ngroups)
    if stats is True:
        print(grouped.size())
        print(max_at)
    grp_list = []
    for group in grouped:
        grp_list.append(group)
    return (grp_list, df_merge)


def full_formated_infputfile(format, file, group_str):
    mol_list = read_file(format, file)
    df = extract_metadata(mol_list, group_str)
    df = formating_dataframe(df) 
    return (df)


def plot_stats(df, plotname, col_name1, col_name2=None):
    if col_name2 is None:
        df.plot.bar(y=col_name1, use_index=True)
        plt.tick_params(left = True, right = False , labelleft = True ,
        labelbottom = False, bottom = True) 
        plt.xlabel("Groups")
        plt.xticks(rotation=30, horizontalalignment="center")
        plt.savefig(plotname +'.png')
    else:
        df.plot.bar(y=[col_name1, col_name2], use_index=True)
        plt.tick_params(left = True, right = False , labelleft = True ,
        labelbottom = False, bottom = True) 
        plt.xlabel("Groups")
        plt.xticks(rotation=30, horizontalalignment="center")
        plt.savefig(plotname +'.png')
        

if __name__ == "__main__":
    file = "RECETOX_GC-EI-MS_20201028.sdf"
    #file = "sample_test.sdf"
    format = "sdf"
    # group_str: Kingdom, Superclass, Class, Subclass, Parent Level 1, Parent Level 2, Parent Level 3
    group_str = "Subclass"
    my_df = full_formated_infputfile(format, file, group_str)
    grp_list, grp_stat = clustering_and_statistics(my_df, GROUP_TYPE="GROUP_TYPE", stats=False)
    print(grp_stat)
    #plot_stats(grp_stat, "plot1", "N_item/grp")
    #plot_stats(grp_stat, "plot2", "Max_n_at")
    plot_stats(grp_stat, "plot_Subclass", "N_item/grp", "Max_n_at")


 
    

