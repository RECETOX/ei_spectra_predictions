from rdkit import Chem
from rdkit.Chem import AllChem
from pathlib import Path
import os
import argparse

from jinja2 import Environment, FileSystemLoader


def CalculateSpinMultiplicity(Mol):
    """Calculate spin multiplicity of a molecule. The spin multiplicity is calculated
    from the number of free radical electrons using Hund's rule of maximum
    multiplicity defined as 2S + 1 where S is the total electron spin. The
    total spin is 1/2 the number of free radical electrons in a molecule.

    Arguments:
    Mol (object): RDKit molecule object.

    Returns:
    int : Spin multiplicity.
    """
    NumRadicalElectrons = 0
    for Atom in Mol.GetAtoms():
        NumRadicalElectrons += Atom.GetNumRadicalElectrons()

    TotalElectronicSpin = NumRadicalElectrons/2
    SpinMultiplicity = 2 * TotalElectronicSpin + 1

    return int(SpinMultiplicity)


def get_method(multiplicity):
    if multiplicity % 2 == 0 :
        method = 'ROHF'
    else:
        method = 'RHF'
    return method


def read_file(file):
    with open(file, 'r') as f:
        lines = f.readlines()
    return lines


def read_file_rdkit(file):
    molfile = Chem.SDMolSupplier(file)
    return molfile


def get_props(mol):
    mol = Chem.AddHs(mol)
    chem_class = mol.GetProp("Class").replace(" ", "_")
    inchikey = mol.GetProp("InChIKey")
    n_atoms = mol.GetNumAtoms()
    molname = mol.GetProp("NAME")
    return chem_class, inchikey, n_atoms, molname


def getRMS(mol, c1, c2):
    rms = AllChem.GetBestRMS(mol, mol, c1, c2)
    return rms


def generate_3D_mol(mol):
    """
    This script was originally written by David Koes, University of Pittsburgh:
    https://github.com/dkoes/rdkit-scripts/blob/master/rdconf.py
    It is licensed under the MIT licence.

    Given a smiles file, generate 3D conformers.
    Energy minimizes and filters conformers to meet energy window and rms constraints.
    """
    maxconfs = 20
    sample = 1
    rmspar = 0.7
    energy = 10
    
    if mol is not None:
        AllChem.SanitizeMol(mol)
        mol = AllChem.AddHs(mol)
        cids = AllChem.EmbedMultipleConfs(mol, int(sample * maxconfs), AllChem.ETKDG())

        cenergy = []
        for conf in cids:
            converged = not AllChem.UFFOptimizeMolecule(mol, confId=conf)
            cenergy.append(AllChem.UFFGetMoleculeForceField(mol, confId=conf).CalcEnergy())

        sortedcids = sorted(cids, key=lambda cid: cenergy[cid])
        if len(sortedcids) > 0:
            mine = cenergy[sortedcids[0]]
        else:
            mine = 0

        written = {}
        final = []
        for conf in sortedcids:
            if len(written) >= maxconfs:
                break
            passed = True
            for seenconf in written.keys():
                rms = getRMS(mol, seenconf, conf)
                if (rms < rmspar) or (energy > 0 and cenergy[conf] - mine > energy):
                    passed = False
                    break
            if passed:
                written[conf] = True
                final.append(mol.GetConformer(conf))
    new_conf = final[0]
    return new_conf


def read_parameters_pbs(param_file):   
    list_pbs_keys = [] 
    lines = read_file(param_file)
    for line in lines:
        if "WALLTIME" in line:
            list_pbs_keys.append(line.split()[2])
        if "NCPUS" in line:
            list_pbs_keys.append(line.split()[2])
        if "MEM" in line:
            list_pbs_keys.append(line.split()[2])
        if "SCRATCH_LOCAL" in line:
            list_pbs_keys.append(line.split()[2])
        if "USER_EMAIL" in line:
            list_pbs_keys.append(line.split()[2])
    return list_pbs_keys


def read_parameters_qcxms(param_file):   
    list_qcxms_keys = [] 
    lines = read_file(param_file)
    for line in lines:
        if "QC_Program" in line:
            list_qcxms_keys.append(line.split()[2])
        if "QC_Level" in line:
            list_qcxms_keys.append(line.split()[2])
        if "ntraj" in line:
            if len(line.split()) == 2:
                list_qcxms_keys.append("")
            else:
                list_qcxms_keys.append(line.split()[2])
        if "tmax" in line:
            list_qcxms_keys.append(line.split()[2])
        if "tinit" in line:
            list_qcxms_keys.append(line.split()[2])
        if "ieeatm" in line:
            list_qcxms_keys.append(line.split()[2])
    return list_qcxms_keys


def load_template(dir, filename):
    template_path = Path(dir).resolve()
    file_loader = FileSystemLoader(template_path)
    environment = Environment(loader=file_loader)
    load_template = environment.get_template(filename)
    return load_template


def write_pbs_from_template(mylist, molname, template_name, file):
    content = template_name.render(MOLNAME=molname, 
                                   WALLTIME=mylist[0],
                                   NCPUS=mylist[1],
                                   MEM=mylist[2],
                                   SCRATCH_LOCAL=mylist[3],
                                   USER_EMAIL=mylist[4])
    with open(file, 'w') as message:
        message.write(content)


def write_qcxms_input_from_template(mylist, template_name, file):
    content = template_name.render(QC_Program=mylist[0],
                                   QC_Level=mylist[1],
                                   ntraj=mylist[2],
                                   tmax=mylist[3],
                                   tinit=mylist[4],
                                   ieeatm=mylist[5])
    with open(file, 'w') as message:
        message.write(content)


def write_gamess_input(multiplicity, mol, molname, mol_input_path):
    conf = generate_3D_mol(mol)
    mol = AllChem.AddHs(mol)
    opt = f""" $CONTRL SCFTYP={get_method(multiplicity)} MULT={multiplicity} NPRINT=-5 RUNTYP=OPTIMIZE $END\n $STATPT OPTTOL=0.0005 NSTEP=100 $END\n $BASIS GBASIS=N31 NGAUSS=6 $END """
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

    molfile = read_file_rdkit(args.sdf_filename)

    for mol in molfile:
        chem_class, inchikey, n_atoms, molname = get_props(mol)
        
        proj_dir = Path(args.project_dirname)
        Path.mkdir(proj_dir, parents=True, exist_ok=True)

        moldir = proj_dir / "classes" / chem_class / inchikey / "Optimization"
        Path.mkdir(moldir, parents=True, exist_ok=True)
        mol_input_path = moldir / (inchikey + ".inp")
        Path.touch(mol_input_path)

        spectradir = proj_dir / "classes" / chem_class / inchikey / "Spectra"
        Path.mkdir(spectradir, parents=True, exist_ok=True)
        spectrum_input_path = spectradir / ("qcxms" + ".in")
                
        multiplicity = CalculateSpinMultiplicity(mol)
        write_gamess_input(multiplicity, mol, molname, mol_input_path)

        spectrum_input_path = spectradir / ("qcxms" + ".in")
        qcxms_template = load_template("templates", "qcxms_input_template.in")
        write_qcxms_input_from_template(read_parameters_qcxms(args.params_filename), qcxms_template, spectrum_input_path)


        mol_pbs_path = moldir / (inchikey + ".pbs")
        pbs_template = load_template("templates", "optimization_gamess_template.pbs")
        write_pbs_from_template(read_parameters_pbs(args.params_filename), inchikey, pbs_template, mol_pbs_path)