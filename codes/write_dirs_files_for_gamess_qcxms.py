from rdkit import Chem
from rdkit.Chem import AllChem
from pathlib import Path
from jinja2 import Environment, FileSystemLoader
import os
import argparse


def calculate_spin_multiplicity(molecule):
    """Calculate spin multiplicity of a molecule. The spin multiplicity is calculated
    from the number of free radical electrons using Hund's rule of maximum
    multiplicity defined as 2S + 1 where S is the total electron spin. The
    total spin is 1/2 the number of free radical electrons in a molecule.

    Args:
        molecule (object): RDKit molecule object.

    Returns:
        int: Spin multiplicity.
    """
    num_radical_electrons = sum(atom.GetNumRadicalElectrons() for atom in molecule.GetAtoms())
    total_electronic_spin = num_radical_electrons / 2
    spin_multiplicity =  int(2 * total_electronic_spin + 1)

    return spin_multiplicity


def get_method(multiplicity):
    """Get the quantum chemistry method based on spin multiplicity.

    Args:
        multiplicity (int): Spin multiplicity.

    Returns:
        str: Quantum chemistry method.
    """
    return 'ROHF' if multiplicity % 2 == 0 else 'RHF'


def read_file(file_path):
    """Read lines from a file.

    Args:
        file_path (str): File path.

    Returns:
        list: List of lines.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    return lines


def read_file_rdkit(file_path):
    molfile = Chem.SDMolSupplier(file_path)
    return molfile


def get_props(molecule):
    """Get molecular properties which include chemical class, inchikey, n_atoms, and molecule name.

    Args:
        molecule (object): RDKit molecule object.

    Returns:
        tuple: Tuple of molecular properties.
    """
    molecule = Chem.AddHs(molecule)
    chem_class = molecule.GetProp("Class").replace(" ", "_")
    inchikey = molecule.GetProp("InChIKey")
    n_atoms = molecule.GetNumAtoms()
    molname = molecule.GetProp("NAME")
    return chem_class, inchikey, n_atoms, molname


def get_rms(molecule, c1, c2):
    """Get RMS value between two conformers.

    Args:
        molecule (object): RDKit molecule object.
        c1 (int): Conformer index 1.
        c2 (int): Conformer index 2.

    Returns:
        float: RMS value.
    """
    rms = AllChem.GetBestRMS(molecule, molecule, c1, c2)
    return rms


def generate_3D_mol(mol):
    """
    This script was originally written by David Koes, University of Pittsburgh:
    https://github.com/dkoes/rdkit-scripts/blob/master/rdconf.py
    It is licensed under the MIT licence.

    Given a smiles file, generate 3D conformers.
    Energy minimizes and filters conformers to meet energy window 
    and rms constraints.
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


def read_parameters(param_file, keys_to_search):
    """
    Read and extract parameter values from a file based on a list of keys.

    Parameters:
    - param_file (str): The path to the parameter file.
    - keys_to_search (list): A list of keys to search for in each line of the file.

    Returns:
    dict: A dictionary containing the extracted parameter values corresponding to the provided keys.
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


def load_template(dir, filename):
    """Load a template from a file.

    Args:
        dir (str): Directory path.
        filename (str): Template file name.

    Returns:
        object: Template object.
    """
    script_directory = os.path.dirname(os.path.abspath(__file__))
    project_directory = os.path.dirname(script_directory)
    # Construct the template directory path based on the current script directory
    template_dir = os.path.join(project_directory, 'templates')

    file_loader = FileSystemLoader(template_dir)
    environment = Environment(loader=file_loader)
    load_template = environment.get_template(filename)
    return load_template


def write_from_template(parameters, template, file):
    """Write file from a template.

    Args:
        parameters (dict): dict containing parameters.
        template (object): Template object.
        file (str): Output file path.
    """
    content = template.render(**parameters)
    
    with open(file, 'w') as message:
        message.write(content)


def write_gamess_input(multiplicity, mol, molname, mol_input_path):
    """Write GAMESS input file.

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
