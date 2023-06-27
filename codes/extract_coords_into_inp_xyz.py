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


PBS_TEMPLATE = load_template("templates", "optimization_gamess_template.pbs")
INP_TEMPLATE = load_template("templates", "gamess_input_template.inp")


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

def coords_as_dataframe(element_symbol, x_coord, y_coord, z_coord):
    df = pd.DataFrame([element_symbol, x_coord, y_coord, z_coord])
    df = df.transpose()
    df[0] = df[0].map("{0:<4s}".format)
    df[1] = pd.to_numeric(df[1], downcast="float").map("{0:>16.10f}".format)
    df[2] = pd.to_numeric(df[2], downcast="float").map("{0:>16.10f}".format)
    df[3] = pd.to_numeric(df[3], downcast="float").map("{0:>16.10f}".format)
    result = df[0] + df[1] + df[2] + df[3]
    return result

def read_coordinates(n_atoms, gamess_output, line_number):
    start, end = compute_start_end_indices(n_atoms, line_number, 4)
    
    data = []
    for i in range(start, end):
        data.append(gamess_output[i].split())

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
    return element_symbol,x_coord,y_coord,z_coord


def extract_geometry_to_xzy(project_dirname, n_atoms, molname, mol_xyz_path, gamess_output, line_number):
    """Extract xyz coordinates from GAMESS output.

    Args:
        project_dirname (str): Name of the project directory.
        n_atoms (int): Number of atoms in the molecule.
        molname (str): Identifier for the molecule.
        mol_xyz_path (str): Path where to write the xzy file
        gamess_output (list[str]): Content of GAMESS output log file.
        line_number (int): Line number of the line currently being processed

    Returns:
        _type_: _description_
    """
    element_symbol, x_coord, y_coord, z_coord = read_coordinates(n_atoms, gamess_output, line_number)

    result = coords_as_dataframe(element_symbol, x_coord, y_coord, z_coord)
    xyz_template = load_template("templates", "structure_input_template.xyz")

    if not Path(mol_xyz_path).exists():
        with open(f"info_extract_coords_{project_dirname}.log", 'a') as f:
            f.write(f"Write XYZ file with optimized coordinates: {mol_xyz_path}")
            f.write(f"{os.linesep}")
        print(f"Write XYZ file with optimized coordinates: {Path(mol_xyz_path).parent}")
        write_xyz_from_template(molname, n_atoms, xyz_template, mol_xyz_path, result)
        return 1
    return 0


def compute_start_end_indices(n_atoms, num, offset):
    start_coord_index = num + offset
    end_coord_index = start_coord_index + n_atoms
    return start_coord_index, end_coord_index


def read_coords_raw(gamess_output, start_coord_index, end_coord_index):
    data = []
    for i in range(start_coord_index, end_coord_index):
        data.append(gamess_output[i])
    return data


def prepare_resubmission(params_filename, inchikey, n_atoms, molname, multiplicity, moldir, mol_input_path, message, gamess_output, num, offset):
    start_coord_index, end_coord_index = compute_start_end_indices(n_atoms, num, offset)
    data = read_coords_raw(gamess_output, start_coord_index, end_coord_index)

    mol_pbs_path = moldir / (inchikey + ".pbs")
    write_pbs_from_template(read_parameters_pbs(params_filename), inchikey, PBS_TEMPLATE, mol_pbs_path)

    mylist = read_parameters_gamess(params_filename)
    start_inp_file = moldir / (inchikey + "_start.inp")

    if not Path(start_inp_file).exists():
        Path(mol_input_path).rename(Path(start_inp_file))
    write_gamess_input_from_template(mylist, multiplicity, molname, INP_TEMPLATE, mol_input_path, data, message)
    return 1


def write_summary_log(project_dirname, count_abn, count_inp, count_xyz):
    """
    """
    with open(f"info_extract_coords_{project_dirname}.log", 'a') as f:
        f.write(f"Number of abnormal molecules: {count_abn}")
        f.write(f"{os.linesep}")
        f.write(f"Number of INP molecules: {count_inp}")
        f.write(f"{os.linesep}")
        f.write(f"Number of XYZ molecules: {count_xyz}")
        f.write(f"{os.linesep}")
        f.write(f"Total of processed molecules= {count_xyz+count_inp+count_abn}")



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
        gamess_log = moldir / (inchikey + ".log")
        mol_input_path = moldir / (inchikey + ".inp")

        spectradir = proj_dir / "classes" / chem_class / inchikey / "Spectra"
        mol_xyz_path = spectradir / (inchikey + ".xyz")

        writer = Chem.SDWriter(str(spectradir / f'{inchikey}.sdf'))
        writer.write(mol)
        writer.close()

        #TODO: Move those into global constants in this script
        message_1 = "EXECUTION OF GAMESS TERMINATED -ABNORMALLY-"
        message_2 = "EXECUTION OF GAMESS TERMINATED NORMALLY"
        message_3 = "EQUILIBRIUM GEOMETRY LOCATED"
        message_4 = "ALWAYS THE LAST POINT COMPUTED!"
        message_5 = "COORDINATES OF ALL ATOMS ARE (ANGS)"

        if Path(gamess_log).exists():
            gamess_output = read_file(gamess_log)

            for num, line in enumerate(gamess_output):
                if message_1 in line: #TODO: Each case should be in an individual function which is documented
                    for num, line in enumerate(gamess_output):
                        if message_5 in line:
                            count_abn += prepare_resubmission(args.params_filename, inchikey, n_atoms, molname, multiplicity, moldir, mol_input_path, message_1, gamess_output, num, 3)

                            with open(info_log, 'a') as f:
                                    f.write(f"Abnormal execution, write INP file  : {Path(gamess_log).parent}")
                                    f.write(f"{os.linesep}")
                            print(f"Abnormal execution, write INP file  : {gamess_log}")

                if message_2 in line:
                    for num, line in enumerate(gamess_output):
                        if message_4 in line:
                            count_inp += prepare_resubmission(args.params_filename, inchikey, n_atoms, molname, multiplicity, moldir, mol_input_path, message_2, gamess_output, num, 4)

                            with open(info_log, 'a') as f:
                                f.write(f"Write INP file with last coordinates: {Path(mol_input_path).parent}")
                                f.write(f"{os.linesep}")
                            print(f"Write INP file with last coordinates: {mol_input_path}")

                if message_3 in line:
                    count_xyz += extract_geometry_to_xzy(args.project_dirname, n_atoms, molname, mol_xyz_path, gamess_output, num)
        else:
            with open(info_log, 'a') as f:
                f.write(f"LOG file does not exist:{Path(gamess_log).parent}")
                f.write(f"{os.linesep}")
            print(f"LOG file does not exist:{Path(gamess_log).parent}")
    
    write_summary_log(args.project_dirname, count_abn, count_inp, count_xyz)


