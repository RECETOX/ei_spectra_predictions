import pandas as pd
import seaborn as sns
import numpy as np
from matplotlib import pyplot as plt
from rdkit import Chem
from itertools import combinations

def is_spectrum_for_compound(compund_name, spectrum_name):
    options = [compund_name + x for x in ["", "_isomer1", "_isomer2", " isomer 1", " isomer 2"]]
    return spectrum_name in options

def get_matching_rows(df, query_name, reference_name):
    return df[df.apply(lambda x: is_spectrum_for_compound(x[query_name], x[reference_name]), axis=1)]

def has_halogen_atoms(mol):
    # Check if the molecule contains any halogen atoms
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in ['F', 'Cl', 'Br', 'I']:
            return True
    return False

def has_atom(mol, atom):
    for mol_atom in mol.GetAtoms():
        if mol_atom.GetSymbol() == atom:
            return True
    return False

def has_organic_atoms(mol):
    # Check if the molecule contains any halogen atoms
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in ['C', 'O', 'N', 'H']:
            return True

    return False

def append_classes(df, left_on):
    molecules = Chem.SDMolSupplier("../../data/RECETOX_GC-EI-MS_20201028.sdf")
    class_names = pd.DataFrame({
        "class" : [m.GetProp("Class") for m in molecules],
        "superclass" : [m.GetProp("Superclass") for m in molecules],
        "molname" : [m.GetProp("NAME") for m in molecules],
        "n_atoms" : [m.GetNumAtoms() for m in molecules],
        "n_bonds" : [m.GetNumBonds() for m in molecules],
        "has_halogen": [has_halogen_atoms(m) for m in molecules],
        "smiles" : [m.GetProp("SMILES") for m in molecules],
        "Cl": [has_atom(m, 'Cl') for m in molecules],
        "Br": [has_atom(m, 'Br') for m in molecules],
        "F": [has_atom(m, 'F') for m in molecules],
        "S": [has_atom(m, 'S') for m in molecules],
        "P": [has_atom(m, 'P') for m in molecules],
        "Si": [has_atom(m, 'Si') for m in molecules],
        "C,O,N,H": [has_organic_atoms(m) for m in molecules],
        "N": [has_atom(m, 'N') for m in molecules],
    })
    merged_df = pd.merge(df, class_names, left_on=left_on, right_on='molname')
    return merged_df

# Define a function to map the true columns to a list of names
def get_true_names(row, df):
    return [col for col in df.columns[11:18] if row[col]]

# Function to split values with commas and create new rows
def split_and_add_rows(df, column_name, split_by):
    df_copy = df.copy()
    df_copy[column_name] = df_copy[column_name].str.split(split_by)
    df_copy = df_copy.explode(column_name).reset_index(drop=True)
    return df_copy

def generate_combinations(df, column_name):
    new_rows = []
    for index, row in df.iterrows():
        values = row[column_name].split(', ')
        all_combos = []
        
        for r in range(2, len(values) + 1):
            for combo in combinations(values, r):
                all_combos.append(', '.join(combo))
        
        for combo in all_combos:
            new_row = row.copy()
            new_row[column_name] = combo
            new_rows.append(new_row)
    
    return pd.DataFrame(new_rows).reset_index(drop=True)

def create_plot(df, path):
    sns.set_style(style='white')
    plt.figure(figsize=(17, 5))

    # Set the color palette
    colors = ['yellow', 'deepskyblue']
    sns.set_palette(sns.color_palette(colors))

    ax = sns.boxplot(x="true_names", y="value", hue="Number", 
                    data=df, hue_order=['CosineHungarian_0.01_0.0_1.0_matches',np.nan],
                    medianprops={'color': 'darkgreen', 'linewidth': 4.0},
                    flierprops={'marker': 'o', 'markersize': 10, 'markerfacecolor': 'none'})  # RUN PLOT   
    ax2 = ax.twinx()

    sns.boxplot(ax=ax2,x='true_names', y='value', hue='Number',
                data=df, hue_order=[np.nan, 'CosineHungarian_0.01_0.0_1.0_scores'], 
                medianprops={'color': 'b', 'linewidth': 4.0}, 
                flierprops={'marker': 'o', 'markersize': 10, 'markerfacecolor': 'none'})

    ax.legend_.remove()
    ax.set_ylim([0, 5])  # Set y-axis limits
    ax.yaxis.set_major_locator(plt.MultipleLocator(1))  # Set major tick marks
    ax.set_ylabel('Match values')  # Set y-axis label
    ax.yaxis.label.set_size(20)  # Set font size of y-axis label
    ax.set_xlabel('Chemical composition', fontsize=20)  # Set x-axis label and font size
    ax.tick_params(axis='x', labelsize=13)  # Set font size of x-axis tick labels
    ax.tick_params(axis='y', labelsize=13)  # Set font size of y-axis tick labels
    ax.yaxis.labelpad = 10
    ax.xaxis.labelpad = 10

    # Create a count for each x-axis label
    count_data = df['true_names'].value_counts().reset_index()
    count_data.columns = ['true_names', 'count']
    count_data = count_data.sort_values(by=['true_names'])
    count_data['count'] = count_data['count'] // 2

    # Remove the original x-axis labels
    ax.set_xticklabels([])

    # Add the count labels to the x-axis
    ax.set_xticks(np.arange(len(count_data)))
    ax.set_xticklabels(count_data['true_names'] + ' (' + count_data['count'].astype(str) + ')', rotation=45, ha='right')

    ax2.legend_.remove()
    ax2.set_ylim([0, 1])  # Set y-axis limits
    ax2.yaxis.set_major_locator(plt.MultipleLocator(0.2))  # Set major tick marks
    ax2.set_ylabel('Score values')  # Set y-axis label
    ax2.yaxis.label.set_size(20)  # Set font size of y-axis label
    ax2.tick_params(axis='y', labelsize=12)  # Set font size of y-axis tick labels
    ax2.yaxis.labelpad = 10

    # Change the legend labels
    handles, labels = ax.get_legend_handles_labels()
    labels[0] = 'Matches'
    labels[1] = 'Scores'
    ax.legend(handles, labels, loc='upper right', fontsize=14)

    plt.savefig(path, dpi=300, bbox_inches='tight')  # SAVE PLOT AS PNG FILE
    plt.show()
    plt.clf()
    plt.close()

def preprocess_data(merged_top5_same):
    # Concatenate the DataFrames in df1_list and add a 'value' column with the value 'matches'.
    df1 = merged_top5_same[['query', 'reference', 'true_names', 'CosineHungarian_0.01_0.0_1.0_matches']].copy()

    # Concatenate the DataFrames in df2_list and add a 'value' column with the value 'scores'.
    df2 = merged_top5_same[['query', 'reference', 'true_names', 'CosineHungarian_0.01_0.0_1.0_scores']].copy()

    # Concatenate df1 and df2 into a single DataFrame.
    df_cat = pd.concat([df1, df2])

    mdf = pd.melt(df_cat, id_vars=['query', 'reference', 'true_names'], var_name=['Number'])      # MELT

    #cleaning data
    mdf = split_and_add_rows(mdf, 'true_names', split_by=', C,O,N,H')
    mdf['true_names'] = mdf['true_names'].replace('', np.nan)
    mdf = mdf.dropna(subset=['value', 'true_names'])

    return mdf