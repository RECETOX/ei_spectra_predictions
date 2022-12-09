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
    return (df)


def formating_dataframe(df, set_idx=False):
    # False for default index
    df = df.transpose()
    df.columns = ["name", "n_at", "group_type", "smile"]
    if set_idx is True:
        df = df.set_index("group_type")
    return (df)


def arrange_dataframe_pivot(df):
    pivot = df.pivot(index="name", columns="group_type", values="n_at")
    pivot.index.name = "name"
    pivot.columns.name = "group_type"
    return (pivot)


def clustering_by_grouptype(df, group_type):
    grouped = df.groupby(group_type)
    grp_list = []
    for group in grouped:
        grp_list.append(group)
    return (grp_list)


def formated_infputfile(format, file, group_str):
    mol_list = read_file(format, file)
    df = extract_metadata(mol_list, group_str)
    df = formating_dataframe(df) 
    return (df)


if __name__ == "__main__":
    # file = "RECETOX_GC-EI-MS_20201028.sdf"
    file = "sample_test.sdf"
    format = "sdf"
    group_str = "Class"
    my_df = formated_infputfile(format, file, group_str)
    # only when default index is set in formating_dataframe
    # my_df = arrange_dataframe_pivot(my_df) 

print(my_df)
# df.plot.bar(y="n_at", use_index=True)
# pivot.plot.bar()
# plt.savefig('test.png')


#df = df['n_at'].sort_values()
#df = df.sort_values(by=['n_at'], ascending=False)
#df = df.groupby('class').groups

#print(df)
#print(df.get_group('Fatty Acyls'))
#print(df.get_group('Phenanthrenes and derivatives'))

    
    

