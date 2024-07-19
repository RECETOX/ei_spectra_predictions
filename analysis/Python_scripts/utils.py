import pandas as pd
import numpy as np
from rdkit import Chem
from typing import List, Tuple
from rdkit.Chem.AllChem import AddHs
from itertools import combinations
from matchms.importing import load_from_msp

from matchms import set_matchms_logger_level
set_matchms_logger_level("ERROR")


def is_spectrum_for_compound(compound_name: str, spectrum_name: str) -> bool:
    """Check if a given spectrum is for a specific compound.

    Args:
        compound_name (str): Name of the compound.
        spectrum_name (str): Name of the spectrum.

    Returns:
        bool: True if the spectrum is for the compound, False otherwise.
    """
    options = generate_compound_name_options(compound_name)
    return spectrum_name in options

def generate_compound_name_options(compound_name: str) -> List[str]:
    """Generate all possible options for a compound name combination.

    Args:
        compound_name (str): Compound name to check for.

    Returns:
        List[str]: possible combinations.
    """
    options = [
        compound_name + x
        for x in ["", "_isomer1", "_isomer2", " isomer 1", " isomer 2"]
    ]
    
    return options


def get_matching_rows(df: pd.DataFrame, query_name: str, reference_name: str) -> pd.DataFrame:
    """Get rows from a DataFrame where the spectrum matches the compound.

    Args:
        df (pd.DataFrame): DataFrame to filter.
        query_name (str): Column name in df for the compound name.
        reference_name (str): Column name in df for the spectrum name.

    Returns:
        pd.DataFrame: Filtered DataFrame where the spectrum matches
        the compound.
    """
    return df[df.apply(lambda x: is_spectrum_for_compound(x[query_name], x[reference_name]), axis=1)]


def has_halogen_atoms(mol: object) -> bool:
    """Check if a molecule contains any halogen atoms.

    Args:
        mol (object): RDKit molecule object.

    Returns:
        bool: True if the molecule has halogen atoms, False otherwise.
    """
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in ['F', 'Cl', 'Br', 'I']:
            return True
    return False


def has_atom(mol: object, atom: str) -> bool:
    """ Check if a molecule contains a specific type of atom.

    Args:
        mol (object): RDKit molecule object.
        atom (str): Symbol of the atom to check.

    Returns:
        bool: True if the molecule contains the specified atom, False otherwise.
    """
    for mol_atom in mol.GetAtoms():
        if mol_atom.GetSymbol() == atom:
            return True
    return False


def has_organic_atoms(mol: object) -> bool:
    """Check if a molecule contains any organic atoms (C, O, N, H).

    Args:
        mol (object): RDKit molecule object.

    Returns:
        bool: True if the molecule has organic atoms, False otherwise.
    """
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in ['C', 'O', 'H']:
            return True
    return False


def get_num_atoms_of_type(molecule: object, atom_type: str) -> int:
    """Get the number of atoms of a specific type in a molecule.

    Args:
        molecule (object): RDKit molecule object.
        atom_type (str): Symbol of the atom type to count.

    Returns:
        int: Number of atoms of the specified type in the molecule.
    """
    return int(len([atom for atom in molecule.GetAtoms() if atom.GetSymbol() == atom_type]))


def append_classes(df: pd.DataFrame, left_on: str) -> pd.DataFrame:
    """Append molecular classes information to a DataFrame based on a specified column.

    Args:
        df (pd.DataFrame): The input DataFrame.
        left_on (str): The column to merge on.

    Returns:
        merged_df (DataFrame): The input DataFrame with additional molecular classes information.
    """
    molecules = Chem.SDMolSupplier("../../analysis/data/experimental/recetox_gc-ei-ms_20201028_properties.sdf")
    class_names = pd.DataFrame({
        "class": [m.GetProp("Class") for m in molecules],
        "superclass": [m.GetProp("Superclass") for m in molecules],
        "subclass": [m.GetProp("Subclass") for m in molecules],
        "molname": [m.GetProp("NAME") for m in molecules],
        "n_atoms": [AddHs(m).GetNumAtoms() for m in molecules],
        "n_bonds": [AddHs(m).GetNumBonds() for m in molecules],
        "has_halogen": [has_halogen_atoms(m) for m in molecules],
        "smiles": [m.GetProp("SMILES") for m in molecules],
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


def sdf_to_dataframe(molecules: Chem.SDMolSupplier) -> pd.DataFrame:
    """Convert an SDF file to a pandas DataFrame.

    Args:
        sdf (Chem.SDMolSupplier): The SDF file to convert.

    Returns:
        pd.DataFrame: The converted DataFrame.
    """
    df =  pd.DataFrame({
        "n_atoms": [int(AddHs(m).GetNumAtoms()) for m in molecules],
        "aromatic_nitrogens": [int(m.GetProp("Aromatic Nitrogens")) for m in molecules],
        "molecular_complexity": [float(m.GetProp("Molecular Complexity")) for m in molecules],
        "molecular_flexibility": [float(m.GetProp("Molecular Flexibility")) for m in molecules],
        "rotatable_bonds": [int(m.GetProp("Rotatable Bonds")) for m in molecules],
        "stereo_centers": [int(m.GetProp("Stereo Centers")) for m in molecules],
        "electronegative_atoms": [int(m.GetProp("Electronegative Atoms")) for m in molecules],
    })
    return df

def load_sdf_into_df(filepath:str):
    db = Chem.SDMolSupplier(filepath, sanitize=True)
    return pd.DataFrame({
    "n_atoms": [int(AddHs(m).GetNumAtoms()) for m in db],
    "class": [m.GetProp("Class") for m in db],
    "superclass": [m.GetProp("Superclass") for m in db],
    "subclass": [m.GetProp("Subclass") for m in db],
    "inchikey": [str(m.GetProp("InChIKey")).split("=")[1] for m in db],
})


def get_true_names(row, df: pd.DataFrame) -> List[str]:
    """Map true columns to a list of names for a given row.

    Args:
       row (pd.Series): The row to map.
       df (pd.DataFrame): The input DataFrame.

    Returns:
       List[str]: List of true column names for the given row.
    """
    return [col for col in df.columns[12:19] if row[col]]


def split_and_add_rows(df: pd.DataFrame, column_name: str, split_by: str) -> pd.DataFrame:
    """Split values in a DataFrame column by a specified delimiter and create new rows.

    Args:
        df (pd.DataFrame): The input DataFrame.
        column_name (str): The column to split and explode.
        split_by (str): The delimiter to split values.

    Returns:
        df_copy (pd.DataFrame): DataFrame with additional rows after splitting and exploding the specified column.
    """
    df_copy = df.copy()
    df_copy[column_name] = df_copy[column_name].str.split(split_by)
    df_copy = df_copy.explode(column_name).reset_index(drop=True)
    return df_copy


def generate_combinations(df: pd.DataFrame, column_name: str) -> pd.DataFrame:
    """
    Generate combinations of values in a DataFrame column and create new rows.

     Args:
        df (pd.DataFrame): The input DataFrame.
        column_name (str): The column to generate combinations for.

    Returns:
        pd.DataFrame: DataFrame with additional rows after generating combinations for the specified column.
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


def preprocess_data(merged_top5_same: pd.DataFrame, cols_to_keep: List[str]) -> pd.DataFrame:
    """Preprocess a DataFrame by concatenating and reshaping it.

    Args:
        merged_top5_same (pd.DataFrame): The input DataFrame.
        cols_to_keep (List[str]): List of additional columns to keep.

    Returns:
        mdf (pd.DataFrame): The preprocessed DataFrame.
    """
    key_cols = ['query', 'reference'] + cols_to_keep

    matches_col = 'matches'
    scores_col = 'scores'

    # Concatenate the DataFrames in df1_list and add a 'value' column with the value 'matches'.
    df1 = merged_top5_same[key_cols + [matches_col]].copy()

    # Concatenate the DataFrames in df2_list and add a 'value' column with the value 'scores'.
    df2 = merged_top5_same[key_cols + [scores_col]].copy()

    # Concatenate df1 and df2 into a single DataFrame.
    df_cat = pd.concat([df1, df2])

    mdf = pd.melt(df_cat, id_vars=key_cols, var_name='Number')      # MELT
    mdf = mdf.dropna()
    return mdf


def clean_chemical_composition_data(mdf: pd.DataFrame) -> pd.DataFrame:
    """Clean chemical composition data in a DataFrame.

    Args:
        mdf (pd.DataFrame): The input DataFrame.

    Returns:
        mdf (pd.DataFrame): The cleaned DataFrame.
    """
    mdf = split_and_add_rows(mdf, 'true_names', split_by=', C,O,N,H')
    mdf['true_names'] = mdf['true_names'].replace('', np.nan)
    mdf = mdf.dropna()
    return mdf


def load_spectra_metadata(file_path: str, metadata_column_name: str) -> Tuple[List, pd.DataFrame, List[str]]:
    """Load spectral data from a file, extract metadata, and return the spectra, the metadata, and a list of names from the metadata.

    Args:
        file_path (str): Path to the file to load spectral data from.
        metadata_column_name (str): Name of the column in the metadata that contains the names.

    Returns:
        Tuple[List, pd.DataFrame, List[str]]: The spectra, the metadata, and the list of names.
    """
    spectra = list(load_from_msp(file_path))
    spectra_metadata = pd.DataFrame.from_dict([x.metadata for x in spectra])
    spectra_metadata.rename(columns={'compound_name': metadata_column_name}, inplace=True)
    spectra_names = spectra_metadata[metadata_column_name].to_list()
    return spectra, spectra_metadata, spectra_names


def append_spectrum_metadata(scores: pd.DataFrame) -> pd.DataFrame:
    """Append spectral metadata to a DataFrame of scores.

    Args:
        scores (pd.DataFrame): DataFrame of scores.

    Returns:
        merged (pd.DataFrame): The scores DataFrame with appended spectral metadata.
    """
    _, spectra_metadata, _ = load_spectra_metadata("../data/filtered/simulated_matchms_filter_1%I_all_peaks.msp", 'query')
    _, reference_spectra_metadata, _ = load_spectra_metadata("../data/experimental/RECETOX_GC-EI_MS_20201028.msp", 'reference')
    merged = scores.merge(spectra_metadata, on="query", how="inner")
    merged.rename(columns={'num_peaks': 'n_peaks_query'}, inplace=True)

    merged = merged.merge(reference_spectra_metadata, on="reference", how="inner")
    merged.rename(columns={'num_peaks': 'n_peaks_reference'}, inplace=True)

    numeric_columns = ['matches', 'n_peaks_query', 'n_peaks_reference']
    merged[numeric_columns] = merged[numeric_columns].apply(pd.to_numeric, errors='coerce')
    return merged


def normalize_df(df: pd.DataFrame, use_nist: bool = True, matches_norm_col: str = 'n_peaks_reference') -> pd.DataFrame:
    """Normalize the 'matches' and 'scores' columns in a DataFrame.

    Args:
        df (pd.DataFrame): DataFrame to normalize.
        use_nist (bool): Whether to multiply the 'scores' column by 1000.
        matches_norm_col (str): Column to normalize the 'matches' column by.

    Returns:
        df (pd.DataFrame): The normalized DataFrame.
    """
    matches_col = 'matches'
    scores_col = 'scores'

    if use_nist:
        df[scores_col] = df[scores_col] * 1000

    if matches_norm_col:
        df[matches_col] = (df[matches_col] / df[matches_norm_col]) * 100
    return df


def load_matchms_scores():
    matchms_scores = pd.read_csv("../data/output_matching/matchms/matchms_tol_0.0035_1%I_all_peaks_with_0s_only_matching.tsv", sep="\t")
    matchms_scores.rename(columns={'CosineHungarian_0.0035_0.0_1.0_scores': 'scores'}, inplace=True)
    matchms_scores.rename(columns={'CosineHungarian_0.0035_0.0_1.0_matches': 'matches'}, inplace=True)
    matchms_scores = append_classes(matchms_scores, 'query')
    matchms_scores = append_spectrum_metadata(matchms_scores)
    return matchms_scores