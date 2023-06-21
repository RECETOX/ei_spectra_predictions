from rdkit import Chem
from pathlib import Path
from write_dirs_files_for_gamess_qcxms import CalculateSpinMultiplicity
from write_dirs_files_for_gamess_qcxms import get_method, get_props
from write_dirs_files_for_gamess_qcxms import read_file_rdkit, load_template
from write_dirs_files_for_gamess_qcxms import read_file, write_pbs_from_template
from write_dirs_files_for_gamess_qcxms import read_parameters_pbs
import pandas as pd
import argparse
import os


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


def write_gamess_input_from_template(mylist, multiplicity, molname, template_name, file, coord, text):
    content = template_name.render(SCFTYP=get_method(multiplicity), 
                                   MULT=multiplicity, 
                                   OPTTOL=mylist[0], NSTEP=mylist[1],
                                   GBASIS=mylist[2], NGAUSS=mylist[3], 
                                   MOLNAME=molname, COORDINATES=coord, text=text, NPRT=-2)
    with open(file, 'w') as message:
        message.write(content)


def write_xyz_from_template(molname, n_atoms, template_name, file, coord):
    content = template_name.render(N_ATOMS=n_atoms, MOLNAME=molname, COORDINATES=coord)
    with open(file, 'w') as message:
        message.write(content)


listarg = argparse.ArgumentParser()
listarg.add_argument('--sdf_filename', type=str)
listarg.add_argument('--params_filename', type=str) 
listarg.add_argument('--project_dirname', type=str)
args = listarg.parse_args()

if __name__ == "__main__":

    molfile = read_file_rdkit(args.sdf_filename)
    count_abn = 0
    count_inp = 0
    count_xyz = 0

    info_log = f"info_extract_coords_{args.project_dirname}.log"
    if Path(info_log).exists():
        os.remove(f"info_extract_coords_{args.project_dirname}.log")

    for mol in molfile:
        mol = Chem.AddHs(mol)
        chem_class, inchikey, n_atoms, molname = get_props(mol)
        multiplicity = CalculateSpinMultiplicity(mol)
        proj_dir = Path(args.project_dirname).resolve()

        moldir = proj_dir / "classes" / chem_class / inchikey / "Optimization"
        log_file = moldir / (inchikey + ".log")
        mol_input_path = moldir / (inchikey + ".inp")

        spectradir = proj_dir / "classes" / chem_class / inchikey / "Spectra"
        mol_xyz_path = spectradir / (inchikey + ".xyz")

        #TODO: Move those into global constants in this script
        message_1 = "EXECUTION OF GAMESS TERMINATED -ABNORMALLY-"
        message_2 = "EXECUTION OF GAMESS TERMINATED NORMALLY"
        message_3 = "EQUILIBRIUM GEOMETRY LOCATED"
        message_4 = "ALWAYS THE LAST POINT COMPUTED!"
        message_5 = "COORDINATES OF ALL ATOMS ARE (ANGS)"
        if Path(log_file).exists():
            lines_log = read_file(log_file)

            for num, line in enumerate(lines_log):
                if message_1 in line: #TODO: Each case should be in an individual function which is documented
                    for num, line in enumerate(lines_log):
                        if message_5 in line:
                            data = []
                            start_coord_index = num + 3
                            end_coord_index = start_coord_index + n_atoms
                            for i in range(start_coord_index, end_coord_index):
                                data.append(lines_log[i])

                            mol_pbs_path = moldir / (inchikey + ".pbs")
                            pbs_template = load_template("templates", "optimization_gamess_template.pbs")
                            write_pbs_from_template(read_parameters_pbs(args.params_filename), inchikey, pbs_template, mol_pbs_path)

                            mylist = read_parameters_gamess(args.params_filename)
                            inp_template = load_template("templates", "gamess_input_template.inp")

                            start_inp_file = moldir / (inchikey + "_start.inp")

                            if not Path(start_inp_file).exists():
                                Path(mol_input_path).rename(Path(start_inp_file))

                            with open(f"info_extract_coords_{args.project_dirname}.log", 'a') as f:
                                    f.write(f"Abnormal execution, write INP file  : {Path(log_file).parent}")
                                    f.write(f"{os.linesep}")
                            print(f"Abnormal execution, write INP file  : {log_file}")
                            write_gamess_input_from_template(mylist, multiplicity, molname, inp_template, mol_input_path, data, message_1)
                            count_abn += 1
                if message_2 in line:
                    for num, line in enumerate(lines_log):
                        if message_4 in line:
                            data = []
                            start_coord_index = num + 4
                            end_coord_index = start_coord_index + n_atoms
                            for i in range(start_coord_index, end_coord_index):
                                data.append(lines_log[i])

                            mol_pbs_path = moldir / (inchikey + ".pbs")
                            pbs_template = load_template("templates", "optimization_gamess_template.pbs")
                            write_pbs_from_template(read_parameters_pbs(args.params_filename), inchikey, pbs_template, mol_pbs_path)

                            mylist = read_parameters_gamess(args.params_filename)
                            inp_template = load_template("templates", "gamess_input_template.inp")
                            start_inp_file = moldir / (inchikey + "_start.inp")

                            if not Path(start_inp_file).exists():
                                Path(mol_input_path).rename(Path(start_inp_file))

                            with open(f"info_extract_coords_{args.project_dirname}.log", 'a') as f:
                                f.write(f"Write INP file with last coordinates: {Path(mol_input_path).parent}")
                                f.write(f"{os.linesep}")
                            print(f"Write INP file with last coordinates: {mol_input_path}")
                            write_gamess_input_from_template(mylist, multiplicity, molname, inp_template, mol_input_path, data, message_2)
                            count_inp += 1
                if message_3 in line:
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
                    xyz_template = load_template("templates", "structure_input_template.xyz")

                    if not Path(mol_xyz_path).exists():
                        with open(f"info_extract_coords_{args.project_dirname}.log", 'a') as f:
                            f.write(f"Write XYZ file with optimized coordinates: {mol_xyz_path}")
                            f.write(f"{os.linesep}")
                        print(f"Write XYZ file with optimized coordinates: {Path(mol_xyz_path).parent}")
                        write_xyz_from_template(molname, n_atoms, xyz_template, mol_xyz_path, result)
                    count_xyz += 1
        else:
            with open(f"info_extract_coords_{args.project_dirname}.log", 'a') as f:
                f.write(f"LOG file does not exist:{Path(log_file).parent}")
                f.write(f"{os.linesep}")
            print(f"LOG file does not exist:{Path(log_file).parent}")
    
    with open(f"info_extract_coords_{args.project_dirname}.log", 'a') as f:
        f.write(f"Number of abnormal molecules: {count_abn}")
        f.write(f"{os.linesep}")
        f.write(f"Number of INP molecules: {count_inp}")
        f.write(f"{os.linesep}")
        f.write(f"Number of XYZ molecules: {count_xyz}")
        f.write(f"{os.linesep}")
        f.write(f"Total of processed molecules= {count_xyz+count_inp+count_abn}")


