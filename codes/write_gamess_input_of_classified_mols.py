from openbabel import pybel
import os
import argparse
from jinja2 import Environment, FileSystemLoader


def read_file(file):
    with open(file, 'r') as f:
        lines = f.readlines()
    return lines

def sort_mols_into_classes(mols) -> dict:
    molecules_by_class = {}
    for mol in mols:
        if molecules_by_class.get(mol.data["Class"]) is None:
            molecules_by_class[mol.data["Class"]] = set()
    
        molecules_by_class[mol.data["Class"]].add(mol)
    return molecules_by_class


def get_method(x):
    if (x.OBMol.GetTotalSpinMultiplicity() % 2) == 0 :
        method = 'ROHF'
    else:
        method = 'RHF'
    return method


listarg = argparse.ArgumentParser()
listarg.add_argument('--filename', type=str) 
args = listarg.parse_args()

def write_pbs_script(template, WALLTIME, NCPUS, MEM, SCRATCH_LOCAL, USER_EMAIL, file):
    content = template.render(MOLNAME=inchikey, WALLTIME=WALLTIME,NCPUS=NCPUS,MEM=MEM,SCRATCH_LOCAL=SCRATCH_LOCAL,USER_EMAIL=USER_EMAIL)
    with open(file, 'w') as message:
        message.write(content)

if __name__ == "__main__":
    file_loader=FileSystemLoader(os.path.join("bin", "templates", ""))
    environment = Environment(loader=file_loader)
    template = environment.get_template("optimization_gamess_template.pbs")

    with open("gamess_submission_parameters.in", 'r') as f:
        lines = f.readlines()
        for line in lines:
            if "WALLTIME" in line:
                WALLTIME = line.split()[2]
            if "NCPUS" in line:
                NCPUS = line.split()[2]
            if "MEM" in line:
                MEM = line.split()[2]
            if "SCRATCH_LOCAL" in line:
                SCRATCH_LOCAL = line.split()[2]
            if "USER_EMAIL" in line:
                USER_EMAIL = line.split()[2]


    mols_from_sdf: list[pybel.Molecule] = list(pybel.readfile("sdf", args.filename))
    classified_mols = sort_mols_into_classes(mols_from_sdf)
    if not os.path.exists("classes"):
        os.mkdir("classes")
        for chem_class in classified_mols.keys():
            outdir = os.path.join("classes", chem_class.replace(" ", "_"))
            os.mkdir(outdir)

            for mol in classified_mols[chem_class]:
                inchikey = mol.data["InChIKey"]
                molecule_dir = os.path.join(outdir, inchikey)
                optim = os.path.join(molecule_dir, "Optimization")
                mass_spec = os.path.join(molecule_dir, "Spectra_simulation")
                os.mkdir(molecule_dir)
                os.mkdir(optim)
                os.mkdir(mass_spec)
                
                with open(os.path.join(optim, inchikey + ".inp"), 'w') as outfile:
                    mol.make3D()
                    opt = f''' $CONTRL SCFTYP={get_method(mol)} MULT={mol.OBMol.GetTotalSpinMultiplicity()} NPRINT=-5 RUNTYP=OPTIMIZE $END\n $STATPT OPTTOL=0.0005 NSTEP=100 NPRT=-2 $END\n $BASIS  GBASIS=N31 NGAUSS=6 $END'''
                    outfile.write(mol.write("inp", opt={"k": opt}))
                    
                    write_pbs_script(template, WALLTIME, NCPUS, MEM, SCRATCH_LOCAL, USER_EMAIL, inchikey, optim)
                        
  