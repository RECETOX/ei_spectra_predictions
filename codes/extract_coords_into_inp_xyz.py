from pathlib import Path
from write_dirs_files_for_gamess_qcxms import CalculateSpinMultiplicity
from write_dirs_files_for_gamess_qcxms import get_method, get_props
from write_dirs_files_for_gamess_qcxms import read_file_rdkit, load_template
from write_dirs_files_for_gamess_qcxms import read_file
import pandas as pd
import argparse


def read_parameters_gamess(param_file):   
    mylist = [] 
    lines = read_file(param_file)
    for line in lines:
        if "OPTTOL" in line:
            mylist.append(line.split()[2])
        if "NSTEP" in line:
            mylist.append(line.split()[2])
        if "GBASIS" in line:
            mylist.append(line.split()[2])
        if "NGAUSS" in line:
            mylist.append(line.split()[2])
    return mylist


def write_inp_from_template(mylist, multiplicity, molname, template_name,file, coord):
    content = template_name.render(SCFTYP=get_method(multiplicity), 
                                   MULT=multiplicity, 
                                   OPTTOL= mylist[0], NSTEP=mylist[1],
                                   GBASIS=mylist[2], NGAUSS=mylist[3], 
                                   MOLNAME=molname, COORDINATES=coord)
    with open(file, 'w') as message:
        message.write(content)


def write_xyz_from_template(molname, n_atoms, template_name, file, coord):
    content = template_name.render(N_ATOMS=n_atoms, MOLNAME=molname, COORDINATES=coord)
    with open(file, 'w') as message:
        message.write(content)


listarg = argparse.ArgumentParser()
listarg.add_argument('--sdf_input', type=str) 
args = listarg.parse_args()

if __name__ == "__main__":

    molfile = read_file_rdkit(args.sdf_input)

    for mol in molfile:
        chem_class, inchikey, n_atoms, molname = get_props(mol)
        multiplicity = CalculateSpinMultiplicity(mol)

        moldir = Path("../classes").resolve() / chem_class / inchikey / "Optimization"
        print(moldir)
        log_file = moldir / (inchikey + ".log")
        print(log_file)
        print(Path(log_file))
        mol_input_path = moldir / (inchikey + ".inp")

        spectradir = Path("../classes").resolve() / chem_class / inchikey / "Spectra"
        mol_xyz_path = spectradir / (inchikey + ".xyz")

        message_1 = "EXECUTION OF GAMESS TERMINATED -ABNORMALLY-"
        message_2 = "EXECUTION OF GAMESS TERMINATED NORMALLY"
        message_3 = "EQUILIBRIUM GEOMETRY LOCATED"
        message_4 = "ALWAYS THE LAST POINT COMPUTED!"
        if Path(log_file).exists():
            lines_log = read_file(log_file)

            for num, line in enumerate(lines_log):
                if message_1 in line:
                    print(f"Abnormal execution, check your input: {log_file}")
                
                if message_2 and message_4 in line:
                    data = []
                    start_coord_index = num + 4
                    end_coord_index = start_coord_index + n_atoms
                    for i in range(start_coord_index, end_coord_index):
                        data.append(lines_log[i])

                    mylist = read_parameters_gamess("all_parameters.in")
                    inp_template = load_template("../templates", "gamess_input_template.inp")
                    
                    start_inp_file = moldir / (inchikey + "_start.in")

                    if not Path(start_inp_file).exists():
                        Path(mol_input_path).rename(Path(start_inp_file))
                    
                    print(f"Write INP file with least coordinates: {mol_input_path}")
                    write_inp_from_template(mylist, multiplicity, molname, inp_template, mol_input_path, data)

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
                    
                    xyz_template = load_template("../templates", "structure_input_template.xyz")
                    
                    print(f"Write XYZ file with optimized coordinates: {mol_xyz_path}")
                    write_xyz_from_template(molname, n_atoms, xyz_template, mol_xyz_path, result)
        else:
            print("hi")
            #print(f"LOG file does not exist:{log_file}")
            