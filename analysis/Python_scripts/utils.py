import pandas as pd
from rdkit import Chem

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
#        "C,O,N,H": [has_organic_atoms(m) for m in molecules],
    })
    merged_df = pd.merge(df, class_names, left_on=left_on, right_on='molname')
    return merged_df

# Define a function to map the true columns to a list of names
def get_true_names(row, df):
    return [col for col in df.columns[11:17] if row[col]]