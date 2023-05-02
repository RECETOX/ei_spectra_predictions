from jinja2 import Environment, FileSystemLoader
import pandas as pd 
from openbabel import pybel
import os


def read_file(file):
    with open(file, 'r') as f:
        lines = f.readlines()
    return lines

def get_method(x):
    if (x.OBMol.GetTotalSpinMultiplicity() % 2) == 0 :
        method = 'ROHF'
    else:
        method = 'RHF'
    return method

def sort_mols_into_classes(mols) -> dict:
    molecules_by_class = {}
    for mol in mols:
        if molecules_by_class.get(mol.data["Class"]) is None:
            molecules_by_class[mol.data["Class"]] = set()
    
        molecules_by_class[mol.data["Class"]].add(mol)
    return molecules_by_class


def get_n_atoms(file):
    lines = read_file(file)
    for num, line in enumerate(lines):
        if "TOTAL NUMBER OF ATOMS" in line:
            n_atoms = int(line.split()[5])
    return n_atoms


mols_from_sdf: list[pybel.Molecule] = list(pybel.readfile("sdf", "sample.sdf"))
classified_mols = sort_mols_into_classes(mols_from_sdf)

for chem_class in classified_mols.keys():
    outdir = os.path.join("classes", chem_class.replace(" ", "_"))

    for mol in classified_mols[chem_class]: 
        inchikey = mol.data["InChIKey"]
        molecule_dir = os.path.join(outdir, inchikey)
        optim = os.path.join(molecule_dir, "Optimization")
        mass_spec = os.path.join(molecule_dir, "Spectra_simulation")
        log_file = os.path.abspath(os.path.join(optim, inchikey + ".log"))

        file_loader = FileSystemLoader(os.path.join("bin", "templates", ""))
        environment = Environment(loader=file_loader)
        
        if os.path.exists(log_file):

            message_1 = "EXECUTION OF GAMESS TERMINATED -ABNORMALLY-"
            message_2 = "EXECUTION OF GAMESS TERMINATED NORMALLY"
            message_3 = "EQUILIBRIUM GEOMETRY LOCATED"
            message_4 = "ALWAYS THE LAST POINT COMPUTED!"
            MOLNAME = inchikey
            lines_log = read_file(log_file)
            n_atoms = get_n_atoms(log_file)  

            for num, line in enumerate(lines_log):
                if message_1 in line:
                    print(f"Abnormal execution, check your input: {optim}")
                
                if message_2 and message_4 in line:
                    data = []
                    start_coord_index = num + 4
                    end_coord_index = start_coord_index + n_atoms
                    for i in range(start_coord_index, end_coord_index):
                        data.append(lines_log[i])
                        
                    lines_param = read_file("all_parameters.in")
                    for line in lines_param:
                        if "OPTTOL" in line:
                            OPTTOL = line.split()[2]
                        if "NSTEP" in line:
                            NSTEP = line.split()[2]
                        if "GBASIS" in line:
                            GBASIS = line.split()[2]
                        if "NGAUSS" in line:
                            NGAUSS = line.split()[2]

                    template = environment.get_template("gamess_input_template.inp")
                    content = template.render(SCFTYP=get_method(mol), MULT=mol.OBMol.GetTotalSpinMultiplicity(), OPTTOL= OPTTOL, NSTEP=NSTEP, GBASIS=GBASIS, NGAUSS=NGAUSS, MOLNAME=inchikey, COORDINATES=data)
                    print(f"Write INP file with least coordinates: {optim}")
                    if not os.path.exists(os.path.join(optim, inchikey + "_start.inp")):
                        os.rename(os.path.join(optim, inchikey + ".inp"), os.path.join(optim, inchikey + "_start.inp"))                  
                    with open(os.path.join(optim, inchikey + ".inp"),'w') as file:
                        file.write(content)

                if message_2 and message_3 in line:
                    data = []
                    start_coord_index = num + 4
                    end_coord_index = start_coord_index + n_atoms
                    for i in range(start_coord_index, end_coord_index):
                        data.append(lines_log[i].split())
                    
                    atomic_number = []
                    element_symbol = []
                    x_coord = []
                    y_coord = []
                    z_coord = []
                    for line in data:
                        element_symbol.append(line[0])
                        atomic_number.append(line[1])
                        x_coord.append(line[2])
                        y_coord.append(line[3])
                        z_coord.append(line[4])
                        
                    df = pd.DataFrame([element_symbol, x_coord, y_coord, z_coord])
                    df = df.transpose()
                    df[0] = df[0].map("{0:<4s}".format)
                    df[1] = pd.to_numeric(df[1], downcast="float").map("{0:>16.10f}".format)
                    df[2] = pd.to_numeric(df[2], downcast="float").map("{0:>16.10f}".format)
                    df[3] = pd.to_numeric(df[3], downcast="float").map("{0:>16.10f}".format)
                    
                    result = df[0] + df[1] + df[2] + df[3]

                    template = environment.get_template("structure_input_template.xyz")
                    content = template.render(N_ATOMS=n_atoms, MOLNAME=MOLNAME, COORDINATES=result)
                    print(f"Write XYZ file with optimized coordinates: {mass_spec}")
                    with open(os.path.join(mass_spec, inchikey + ".xyz"),'w') as file:
                        file.write(content)
        else:
            print(f"LOG file does not exist:{optim}")




