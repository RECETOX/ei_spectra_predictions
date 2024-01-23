import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem.AllChem import AddHs
from itertools import combinations
from matchms.importing import load_from_msp

from matchms import set_matchms_logger_level
set_matchms_logger_level("ERROR")

def is_spectrum_for_compound(compund_name, spectrum_name):
    options = [compund_name + x for x in ["", "_isomer1", "_isomer2", " isomer 1", " isomer 2"]]
    return spectrum_name in options

def get_matching_rows(df, query_name, reference_name):
    return df[df.apply(lambda x: is_spectrum_for_compound(x[query_name], x[reference_name]), axis=1)]

def has_halogen_atoms(mol):
    """
    Check if a molecule contains any halogen atoms.

    Parameters:
    - mol (Chem.Mol): RDKit molecule object.

    Returns:
    - bool: True if the molecule has halogen atoms, False otherwise.
    """
    # Check if the molecule contains any halogen atoms
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in ['F', 'Cl', 'Br', 'I']:
            return True
    return False

def has_atom(mol, atom):
    """
    Check if a molecule contains a specific type of atom.

    Parameters:
    - mol (Chem.Mol): RDKit molecule object.
    - atom (str): Symbol of the atom to check.

    Returns:
    - bool: True if the molecule contains the specified atom, False otherwise.
    """
    for mol_atom in mol.GetAtoms():
        if mol_atom.GetSymbol() == atom:
            return True
    return False

def has_organic_atoms(mol):
    """
    Check if a molecule contains any organic atoms (C, O, N, H).

    Parameters:
    - mol (Chem.Mol): RDKit molecule object.

    Returns:
    - bool: True if the molecule has organic atoms, False otherwise.
    """
    # Check if the molecule contains any halogen atoms
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in ['C', 'O', 'H']:
            return True

    return False

def get_num_atoms_of_type(molecule, atom_type):
    return int(len([atom for atom in molecule.GetAtoms() if atom.GetSymbol() == atom_type]))

def append_classes(df, left_on):
    """
    Append molecular classes information to a DataFrame based on a specified column.

    Parameters:
    - df (pd.DataFrame): The input DataFrame.
    - left_on (str): The column to merge on.

    Returns:
    - pd.DataFrame: The input DataFrame with additional molecular classes information.
    """
    molecules = Chem.SDMolSupplier("../../data/recetox_gc-ei-ms_20201028_custom.sdf")
    class_names = pd.DataFrame({
        "class" : [m.GetProp("Class") for m in molecules],
        "superclass" : [m.GetProp("Superclass") for m in molecules],
        "subclass" : [m.GetProp("Subclass") for m in molecules],
        "molname" : [m.GetProp("NAME") for m in molecules],
        "n_atoms" : [AddHs(m).GetNumAtoms() for m in molecules],
        "n_bonds" : [AddHs(m).GetNumBonds() for m in molecules],
        "has_halogen": [has_halogen_atoms(m) for m in molecules],
        "smiles" : [m.GetProp("SMILES") for m in molecules],
        "Cl": [has_atom(m, 'Cl') for m in molecules],
        "Br": [has_atom(m, 'Br') for m in molecules],
        "F": [has_atom(m, 'F') for m in molecules],
        "S": [has_atom(m, 'S') for m in molecules],
        "P": [has_atom(m, 'P') for m in molecules],
        "Si": [has_atom(m, 'Si') for m in molecules],
        "C,O,N,H": [has_organic_atoms(m) or has_atom(m, 'N') for m in molecules],
        "N": [has_atom(m, 'N') for m in molecules],
        "aromatic_nitrogens": [int(m.GetProp("Aromatic Nitrogens")) for m in molecules],
        "molecular_complexity": [float(m.GetProp("Molecular Complexity")) for m in molecules],
        "molecular_flexibility": [float(m.GetProp("Molecular Flexibility")) for m in molecules],
        "rotatable_bonds": [int(m.GetProp("Rotatable Bonds")) for m in molecules],
        "amines": [int(m.GetProp("Amines")) for m in molecules],
        "amides": [int(m.GetProp("Amides")) for m in molecules],
        "non_ch_atoms": [int(m.GetProp("Non-C/H Atoms")) for m in molecules],
        "stereo_centers": [int(m.GetProp("Stereo Centers")) for m in molecules],
        "electronegative_atoms": [int(m.GetProp("Electronegative Atoms")) for m in molecules],
        "aromatic_atoms": [int(m.GetProp("Aromatic Atoms")) for m in molecules],
        "aromatic_rings": [int(m.GetProp("Aromatic Rings")) for m in molecules],
        "n_chlorine": [get_num_atoms_of_type(m, "Cl") for m in molecules],
        "n_bromine": [get_num_atoms_of_type(m, "Br") for m in molecules],
        "n_fluorine": [get_num_atoms_of_type(m, "F") for m in molecules],
        "n_sulfur": [get_num_atoms_of_type(m, "S") for m in molecules],
        "n_phosphorus": [get_num_atoms_of_type(m, "P") for m in molecules],
        "composition": [",".join(sorted(set([a.GetSymbol() for a in AddHs(m).GetAtoms()]))) for m in molecules]
    })
    merged_df = pd.merge(df, class_names, left_on=left_on, right_on='molname')
    return merged_df

# Define a function to map the true columns to a list of names
def get_true_names(row, df):
    """
    Map true columns to a list of names for a given row.

    Parameters:
    - row: The row in the DataFrame.
    - df (pd.DataFrame): The DataFrame.

    Returns:
    - list: List of true column names for the given row.
    """
    return [col for col in df.columns[12:19] if row[col]]

# Function to split values with commas and create new rows
def split_and_add_rows(df, column_name, split_by):
    """
    Split values in a DataFrame column by a specified delimiter and create new rows.

    Parameters:
    - df (pd.DataFrame): The input DataFrame.
    - column_name (str): The column to split and explode.
    - split_by (str): The delimiter to split values.

    Returns:
    - pd.DataFrame: DataFrame with additional rows after splitting and exploding the specified column.
    """
    df_copy = df.copy()
    df_copy[column_name] = df_copy[column_name].str.split(split_by)
    df_copy = df_copy.explode(column_name).reset_index(drop=True)
    return df_copy

def generate_combinations(df, column_name):
    """
    Generate combinations of values in a DataFrame column and create new rows.

    Parameters:
    - df (pd.DataFrame): The input DataFrame.
    - column_name (str): The column to generate combinations for.

    Returns:
    - pd.DataFrame: DataFrame with additional rows after generating combinations for the specified column.
    """
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

def preprocess_data(merged_top5_same, cols_to_keep):
    key_cols = ['query', 'reference'] + cols_to_keep
    
    matches_col = 'CosineHungarian_0.01_0.0_1.0_matches'
    scores_col = 'CosineHungarian_0.01_0.0_1.0_scores'
    
    # Concatenate the DataFrames in df1_list and add a 'value' column with the value 'matches'.
    df1 = merged_top5_same[key_cols + [matches_col]].copy()

    # Concatenate the DataFrames in df2_list and add a 'value' column with the value 'scores'.
    df2 = merged_top5_same[key_cols + [scores_col]].copy()

    # Concatenate df1 and df2 into a single DataFrame.
    df_cat = pd.concat([df1, df2])

    mdf = pd.melt(df_cat, id_vars=key_cols, var_name='Number')      # MELT
    mdf = mdf.dropna()
    return mdf

def clean_chemical_composition_data(mdf):
    mdf = split_and_add_rows(mdf, 'true_names', split_by=', C,O,N,H')
    mdf['true_names'] = mdf['true_names'].replace('', np.nan)
    mdf = mdf.dropna()
    return mdf

def load_spectra_metadata(file_path, metadata_column_name):
    spectra = list(load_from_msp(file_path))
    spectra_metadata = pd.DataFrame.from_dict([x.metadata for x in spectra])
    spectra_metadata.rename(columns={'compound_name': metadata_column_name}, inplace=True)
    spectra_names = spectra_metadata[metadata_column_name].to_list()
    return spectra, spectra_metadata, spectra_names


def append_spectrum_metadata(scores: pd.DataFrame):
    _ , spectra_metadata, _ = load_spectra_metadata("../data/filtered/simulated_matchms_filter_1%I_all_peaks.msp", 'query')
    _ , reference_spectra_metadata, _ = load_spectra_metadata("../data/experimental/RECETOX_GC-EI_MS_20201028.msp", 'reference')
    merged = scores.merge(spectra_metadata, on="query", how="inner")
    merged.rename(columns={'num_peaks': 'n_peaks_query'}, inplace=True)

    merged = merged.merge(reference_spectra_metadata, on="reference", how="inner")
    merged.rename(columns={'num_peaks': 'n_peaks_reference'}, inplace=True)

    numeric_columns = ['CosineHungarian_0.01_0.0_1.0_matches', 'n_peaks_query', 'n_peaks_reference']
    merged[numeric_columns] = merged[numeric_columns].apply(pd.to_numeric, errors='coerce')
    return merged

def normalize_df(df: pd.DataFrame, use_nist: bool = True, matches_norm_col: str = 'n_peaks_reference'):
    matches_col = 'CosineHungarian_0.01_0.0_1.0_matches'
    scores_col = 'CosineHungarian_0.01_0.0_1.0_scores'
    
    if use_nist:
        df[scores_col] = df[scores_col] * 1000
    
    if matches_norm_col:
        df[matches_col] = (df[matches_col] / df[matches_norm_col]) * 100
    
    return df