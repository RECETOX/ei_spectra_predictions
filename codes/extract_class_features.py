from openbabel import pybel
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def read_file(format, filename):
    pybel.ob.obErrorLog.SetOutputLevel(0)
    mols = list(pybel.readfile(format, filename))
    return (mols)


def build_dataframe_of_extracted_metadata(mol_list, group_str):
    data, mol_name, smile, GROUP_TYPE, n_at = [],[],[],[],[]
    for line in mol_list:
        data.append(line.data)
        mol_name.append(line.data["NAME"])
        GROUP_TYPE.append(line.data[group_str])
        smile.append(line.data["SMILES"])
        n_at.append(line.OBMol.NumAtoms())
    df = pd.DataFrame([mol_name, n_at, GROUP_TYPE, smile])
    pd.set_option('display.colheader_justify', 'center') 
    return (df)


def rearrange_dataframe(df):
    df = df.transpose()
    df.columns = ["NAME", "N_ATOMS", "GROUP_TYPE", "SMILE"] # df = df.sort_values(by=['N_at'], ascending=False)
    #df = df.set_index("GROUP_TYPE")
    return (df)


def inputdata_reformatted(format, file, group_str):
    mol_list = read_file(format, file)
    df = build_dataframe_of_extracted_metadata(mol_list, group_str)
    df = rearrange_dataframe(df) 
    return (df)


def grouped_by_item(df, GROUP_TYPE):
    grouped = df.groupby(GROUP_TYPE)
    return (grouped)

def group_stats(df, GROUP_TYPE):
    group_df = grouped_by_item(df, GROUP_TYPE)
    df_col = group_df[["N_ATOMS"]]
    size_df = group_df.size().to_frame(name="ITEMS/GRP")
    max_at = df_col.max()
    df_merge = size_df.merge(max_at, how='inner', on='GROUP_TYPE')
    df_merge = df_merge.sort_values(by=['ITEMS/GRP'], ascending=[False])
    df_merge.rename(columns = {'N_ATOMS':'MAX_N_ATOMS'}, inplace = True)
    max_at.rename(columns = {'N_ATOMS':'MAX_N_ATOMS'}, inplace = True)
    print("Number of groups found --> ", group_df.ngroups)
    # print(group_df.size())
    # print(max_at)  
    return (df_merge)


def plot_stats(df, plotname, col_name1, col_name2=None):
    df.plot.bar(y=[col_name1, col_name2], use_index=True)
    plt.tick_params(left = True, right = False , labelleft = True ,
    labelbottom = False, bottom = True) 
    plt.xlabel("Groups")
    plt.xticks(rotation=30, horizontalalignment="center")
    plt.savefig(plotname +'.png')
        


if __name__ == "__main__":
    file = "RECETOX_GC-EI-MS_20201028.sdf"
    #file = "sample.sdf"
    format = "sdf"
    # group_str: Kingdom, Superclass, Class, Subclass, Parent Level 1, Parent Level 2, Parent Level 3
    group_str = "Class"
    my_df = inputdata_reformatted(format, file, group_str)
    grp_stat = group_stats(my_df, GROUP_TYPE="GROUP_TYPE")
    print(grp_stat)

    # Saves SMILES of every group in separate files 
    my_df = my_df.set_index("SMILE")
    group_df = grouped_by_item(my_df, GROUP_TYPE="GROUP_TYPE")
    dic_g = group_df.groups
    keys = list(dic_g.keys())
    for k in keys:
        dic_val = dic_g[k].values
        outfile = open(group_str + '_' + k +'.smi', 'w+')
        for x in range(len(dic_val)):
            outfile.write(dic_val[x] + '\n')
        outfile.close()


    # Reads SMILES files and convert to XYZ
    for k in keys:
        input_smi = group_str + '_' + k +'.smi'
        outfilexyz = group_str + '_' + k +'.xyz'
        mol_ls = read_file("smi", input_smi)
        out = pybel.Outputfile("xyz", outfilexyz, overwrite=True)
        for line in mol_ls:
            line.make3D()
            out.write(line)




    # for name, group in group_df:
    #     plot = sns.scatterplot(data=group, x="NAME", y="N_ATOMS")
    #     plot.set(xticklabels=[]) 
    #     plot.set(xlabel=None)
    #     plot.tick_params(bottom=False)
    #     plt.savefig(name +'.png')




 
    

