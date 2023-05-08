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


def read_parameters(param_file):   
    mylist = [] 
    lines = read_file(param_file)
    for line in lines:
        if "WALLTIME" in line:
            mylist.append(line.split()[2])
        if "NCPUS" in line:
            mylist.append(line.split()[2])
        if "MEM" in line:
            mylist.append(line.split()[2])
        if "SCRATCH_LOCAL" in line:
            mylist.append(line.split()[2])
        if "USER_EMAIL" in line:
            mylist.append(line.split()[2])
    return mylist


def write_pbs_from_template(mylist, molname, templatename, file):
    content = templatename.render(MOLNAME=molname, WALLTIME=mylist[0],
                                  NCPUS=mylist[1], MEM=mylist[2], 
                                  SCRATCH_LOCAL=mylist[3], USER_EMAIL=mylist[4])
    with open(file, 'w') as message:
        message.write(content)


listarg = argparse.ArgumentParser()
listarg.add_argument('--sdf_input', type=str) 
listarg.add_argument('--params_input', type=str) 
args = listarg.parse_args()

if __name__ == "__main__":

    template_path = Path("templates").resolve()
    file_loader = FileSystemLoader(template_path)
    environment = Environment(loader=file_loader)

    molfile = Chem.SDMolSupplier(args.sdf_input)

    for mol in molfile:
        mol = Chem.AddHs(mol)
        chem_class = mol.GetProp("Class").replace(" ", "_")
        inchikey = mol.GetProp("InChIKey")
        name = mol.GetProp("NAME")

        moldir = Path.cwd() / "classes" / chem_class / inchikey / "Optimization"
        Path.mkdir(moldir, parents=True, exist_ok=True)
        mol_input_path = moldir / (inchikey + ".inp")
        Path.touch(mol_input_path)

        spectradir = Path.cwd() / "classes" / chem_class / inchikey / "Spectra"
        Path.mkdir(spectradir, parents=True, exist_ok=True)
        spectrum_input_path = spectradir / ("qcxms" + ".in")
        Path.touch(spectrum_input_path)

        AllChem.EmbedMolecule(mol, maxAttempts=10000, useRandomCoords=False)
        conf = mol.GetConformer()
        multiplicity = CalculateSpinMultiplicity(mol)
        opt = f""" $CONTRL SCFTYP={get_method(multiplicity)} \
            MULT={multiplicity} NPRINT=-5 RUNTYP=OPTIMIZE $END\n \
            $STATPT OPTTOL=0.0005 NSTEP=100 NPRT=-2 $END\n $BASIS \
            GBASIS=N31 NGAUSS=6 $END """

        with open(mol_input_path, 'w') as outfile:
            outfile.write(f"{opt}{os.linesep}")
            outfile.write(f"{os.linesep}")
            outfile.write(f" $DATA{os.linesep}")
            outfile.write(f"{name}{os.linesep}")
            outfile.write(f"C1{os.linesep}")
            
            for (i, a) in enumerate(mol.GetAtoms()):
                pos = conf.GetAtomPosition(i)
                outfile.write(f"{a.GetSymbol()}\t{a.GetAtomicNum()}{pos.x:>16.10f}{pos.y:>16.10f}{pos.z:>16.10f}{os.linesep}")
            
            outfile.write(f"$END{os.linesep}")

        pbs_template = environment.get_template("optimization_gamess_template.pbs")
        mol_pbs_path = moldir / (inchikey + ".pbs")

        write_pbs_from_template(read_parameters(args.params_input), inchikey, pbs_template, mol_pbs_path)