from openbabel import pybel
import pandas as pd
import matplotlib.pyplot as plt


def read_file(format, filename):
    mols = list(pybel.readfile(format, filename))
    return (mols)


def extract_metadata(mol_list, group_str):
    data = []
    mol_name = []
    smile = []
    group_type = []
    n_at = []
    for line in mol_list:
        data.append(line.data)
        mol_name.append(line.data["NAME"])
        group_type.append(line.data[group_str])
        smile.append(line.data["SMILES"])
        n_at.append(line.OBMol.NumAtoms())

    df = pd.DataFrame([mol_name, n_at, group_type, smile])
    pd.set_option('display.colheader_justify', 'center') 
    return (df)


def formating_dataframe(df, set_idx=False):
    # False for default index
    df = df.transpose()
    df.columns = ["Name", "N_at", "Group_type", "Smile"]
    df = df.sort_values(by=['N_at'], ascending=False)
    if set_idx is True:
        df = df.set_index("Group_type")
    return (df)


def arrange_dataframe_pivot(df):
    pivot = df.pivot(index="Name", columns="Group_type", values="N_at")
    pivot.index.name = "Name"
    pivot.columns.name = "Group_type"
    return (pivot)


def clustering_and_statistics(df, group_type, stats=False):
    grouped = df.groupby(group_type)
    df_col = grouped[["N_at"]]
    if stats is True:
        #print("Number of groups found --> ", grouped.ngroups)
        #print("Number of items per group --> ", grouped.size(), df_col.max())
        print(type(grouped.ngroups))
        print(type(grouped.size().to_frame()))
        print(type(df_col.max()))
        print(grouped.ngroups)
        print(grouped.size().to_frame(name="N_comp"))
        print(df_col.max())

    grp_list = []
    for group in grouped:
        grp_list.append(group)
    
    return (grp_list)


def full_formated_infputfile(format, file, group_str):
    mol_list = read_file(format, file)
    df = extract_metadata(mol_list, group_str)
    df = formating_dataframe(df) 
    return (df)


if __name__ == "__main__":
    # file = "RECETOX_GC-EI-MS_20201028.sdf"
    file = "sample_test.sdf"
    format = "sdf"
    # group_str: Kingdom, Superclass, Class, Subclass, Parent Level 1, Parent Level 2, Parent Level 3
    group_str = "Superclass"
    my_df = full_formated_infputfile(format, file, group_str)
    # arrange_dataframe_pivot(my_df)  only works for default set of index  in formating_dataframe()
    # my_df = arrange_dataframe_pivot(my_df) 
    grouped = clustering_and_statistics(my_df, group_type="Group_type", stats=True)



# df.plot.bar(y="n_at", use_index=True)
# pivot.plot.bar()
# plt.savefig('test.png')


#df = df['n_at'].sort_values()
#df = df.groupby('class').groups

#print(df)
#print(df.get_group('Fatty Acyls'))
#print(df.get_group('Phenanthrenes and derivatives'))

    
    

