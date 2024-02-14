from rdkit import Chem
from rdkit.Chem import AllChem
from pathlib import Path
from jinja2 import Environment, FileSystemLoader
from typing import List, Dict, Tuple
import os
import argparse


def calculate_spin_multiplicity(molecule: object) -> int:
    """Calculate spin multiplicity of a molecule.

    The spin multiplicity is calculated from the number of free radical electrons using Hund's
    rule of maximum multiplicity defined as 2S + 1 where S is the total electron spin. The
    total spin is 1/2 the number of free radical electrons in a molecule.

    Args:
        molecule (object): RDKit molecule object.

    Returns:
        spin_multiplicity (int): Spin multiplicity.
    """
    num_radical_electrons = sum(atom.GetNumRadicalElectrons() for atom in molecule.GetAtoms())
    total_electronic_spin = num_radical_electrons / 2
    spin_multiplicity = int(2 * total_electronic_spin + 1)
    return spin_multiplicity


def get_method(multiplicity: int) -> str:
    """Get the quantum chemistry method based on spin multiplicity.

    Args:
        multiplicity (int): Spin multiplicity.

    Returns:
        str: Quantum chemistry method.
    """
    return 'ROHF' if multiplicity % 2 == 0 else 'RHF'


def read_file(file_path: str) -> List[str]:
    """Read lines from a file.

    Args:
        file_path (str): File path.

    Returns:
        lines (List[str]): List of lines.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()
    return lines


def read_file_rdkit(file_path: str) -> object:
    """Read an SDF file using RDKit.

    Args:
        file_path (str): File path.

    Returns:
        molfile (object): RDKit molecule object.
    """
    molfile = Chem.SDMolSupplier(file_path)
    return molfile


def get_props(molecule: object) -> Tuple[str, str, int, str]:
    """Get molecular properties which include chemical class, inchikey, n_atoms, and molecule name.

    Args:
        molecule (object): RDKit molecule object.

    Returns:
        tuple: Tuple of molecular properties (chemical class, InChIKey, n_atoms, molecule name).
    """
    molecule = Chem.AddHs(molecule)
    chem_class = molecule.GetProp("Class").replace(" ", "_")
    inchikey = molecule.GetProp("InChIKey")
    n_atoms = molecule.GetNumAtoms()
    molname = molecule.GetProp("NAME")
    return chem_class, inchikey, n_atoms, molname


def get_rms(molecule: object, c1: int, c2: int) -> float:
    """Get RMS value between two conformers.

    Args:
        molecule (object): RDKit molecule object.
        c1 (int): Conformer index 1.
        c2 (int): Conformer index 2.

    Returns:
        rms (float): RMS value.
    """
    rms = AllChem.GetBestRMS(molecule, molecule, c1, c2)
    return rms


def generate_3D_mol(mol: object) -> object:
    """ Generate 3D conformers for a molecule.

    Generates 3D conformers for a given RDKit molecule object, performs
    energy minimization, and filters the conformers to meet certain
    energy window and RMS constraints.

    This script was originally written by David Koes, University of Pittsburgh:
    https://github.com/dkoes/rdkit-scripts/blob/master/rdconf.py
    It is licensed under the MIT licence.

    Args:
        mol (object): RDKit molecule object.

    Returns:
        new_conf (object): RDKit molecule object with 3D conformers.
    """
    max_confs = 20
    sample = 1
    rms_par = 0.7
    energy = 10

    if mol is not None:
        AllChem.SanitizeMol(mol)
        mol = AllChem.AddHs(mol)
        cids = AllChem.EmbedMultipleConfs(mol, int(sample * max_confs), AllChem.ETKDG())

        cenergy = []
        for conf in cids:
            converged = not AllChem.UFFOptimizeMolecule(mol, confId=conf)
            cenergy.append(AllChem.UFFGetMoleculeForceField(mol, confId=conf).CalcEnergy())

        sorted_cids = sorted(cids, key=lambda cid: cenergy[cid])

        if len(sorted_cids) > 0:
            mine = cenergy[sorted_cids[0]]
        else:
            mine = 0

        written = {}
        final = []
        for conf in sorted_cids:
            if len(written) >= max_confs:
                break
            passed = True
            for seenconf in written.keys():
                rms = get_rms(mol, seenconf, conf)
                if (rms < rms_par) or (energy > 0 and cenergy[conf] - mine > energy):
                    passed = False
                    break
            if passed:
                written[conf] = True
                final.append(mol.GetConformer(conf))
    new_conf = final[0]
    return new_conf


def read_parameters(param_file: str, keys_to_search: List[str]) -> Dict[str, str]:
    """
    Read and extract parameter values from a file based on a list of keys.

    For each key, the function looks for a line in the file that contains the key and
    extracts the third word in the line as the value for that key. If the key is "ntraj"
    and the line contains only two words, the value is set to an empty string.

    Args:
        param_file (str): The path to the parameter file.
        keys_to_search (list): List of keys to search.

    Returns:
        parameters (Dict[str, str]): Dictionary of extracted parameters.
    """
    parameters = {}
    lines = read_file(param_file)
    for line in lines:
        for key in keys_to_search:
            if key in line:
                if key == "ntraj":
                    if len(line.split()) == 2:
                        parameters[key] = ""
                    else:
                        parameters[key] = line.split()[2]
                else:
                    parameters[key] = line.split()[2]
    return parameters


# dir function argument is not used in the function
def load_template(filename: str) -> object:
    """Load a template from a file.

    Args:
        filename (str): Template file name.

    Returns:
        load_template (object): Template object.
    """
    script_directory = os.path.dirname(os.path.abspath(__file__))
    project_directory = os.path.dirname(script_directory)
    # Construct the template directory path based on the current script directory
    template_dir = os.path.join(project_directory, 'templates')

    file_loader = FileSystemLoader(template_dir)
    environment = Environment(loader=file_loader)
    load_template = environment.get_template(filename)
    return load_template


def write_from_template(parameters: dict, template: object, file: str) -> None:
    """Write file from a template.

    Args:
        parameters (dict): dict containing parameters.
        template (object): Template object.
        file (str): Output file path.
    """
    content = template.render(**parameters)

    with open(file, 'w') as message:
        message.write(content)


def write_gamess_input(multiplicity: int, mol: object, molname: str, mol_input_path: str) -> None:
    """Write GAMESS input file for a given molecule.

    This function generates a 3D conformation of the molecule, adds Hydrogen atoms,
    and writes the necessary GAMESS input options. It then iterates over the atoms in
    the molecule, writing each atom's symbol, atomic number, and x, y, and z coordinates to the file.

    Args:
        multiplicity (int): Spin multiplicity.
        mol (object): RDKit molecule object.
        molname (str): Molecular name.
        mol_input_path (str): Output file path.
    """
    conf = generate_3D_mol(mol)
    mol = AllChem.AddHs(mol)
    opt = f""" $CONTRL SCFTYP={get_method(multiplicity)} MULT={multiplicity} NPRINT=-5 RUNTYP=OPTIMIZE $END\n $STATPT OPTTOL=0.0005 NSTEP=100 $END\n $BASIS GBASIS=N31 NGAUSS=6 $END\n $SYSTEM MWORDS=128 $END """
    with open(mol_input_path, 'w') as outfile:
        outfile.write(f"{opt}{os.linesep}")
        outfile.write(f"{os.linesep}")
        outfile.write(f" $DATA{os.linesep}")
        outfile.write(f"{molname}{os.linesep}")
        outfile.write(f"C1{os.linesep}")
        for (i, a) in enumerate(mol.GetAtoms()):
            pos = conf.GetAtomPosition(i)
            outfile.write(f"{a.GetSymbol()}\t{a.GetAtomicNum()}{pos.x:>16.10f}{pos.y:>16.10f}{pos.z:>16.10f}{os.linesep}")
        outfile.write(f"$END{os.linesep}")


listarg = argparse.ArgumentParser()
listarg.add_argument('--sdf_filename', type=str)
listarg.add_argument('--params_filename', type=str)
listarg.add_argument('--project_dirname', type=str)
args = listarg.parse_args()

if __name__ == "__main__":

    mol_file = read_file_rdkit(args.sdf_filename)

    for mol in mol_file:
        chem_class, inchikey, n_atoms, molname = get_props(mol)

        # TODO: Extract function to construct the paths as they are often re-used
        proj_dir = Path(args.project_dirname)
        Path.mkdir(proj_dir, parents=True, exist_ok=True)

        mol_dir = proj_dir / "classes" / chem_class / inchikey / "Optimization"
        Path.mkdir(mol_dir, parents=True, exist_ok=True)
        mol_input_path = mol_dir / (inchikey + ".inp")

        spectra_dir = proj_dir / "classes" / chem_class / inchikey / "Spectra"
        Path.mkdir(spectra_dir, parents=True, exist_ok=True)

        multiplicity = calculate_spin_multiplicity(mol)
        if not Path(mol_input_path).exists():
            Path.touch(mol_input_path, exist_ok=True)
            write_gamess_input(multiplicity, mol, molname, mol_input_path)

        spectrum_input_path = spectra_dir / ("qcxms" + ".in")
        qcxms_template = load_template("templates", "qcxms_input_template.in")
        if not Path(spectrum_input_path).exists():
            qcxms_params = read_parameters(args.params_filename, ["QC_Program", "QC_Level", "ntraj", "tmax", "tinit", "ieeatm"])
            write_from_template(parameters=qcxms_params, template=qcxms_template, file=spectrum_input_path)

        mol_pbs_path = mol_dir / (inchikey + ".pbs")
        pbs_template = load_template("templates", "optimization_gamess_template.pbs")
        if not Path(mol_pbs_path).exists():
            pbs_params = read_parameters(args.params_filename, ["WALLTIME", "NCPUS", "MEM", "SCRATCH_LOCAL", "USER_EMAIL"])
            write_from_template(parameters={"MOLNAME": inchikey, **pbs_params}, template=pbs_template, file=mol_pbs_path)
