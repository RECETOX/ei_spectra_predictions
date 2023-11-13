from rdkit import Chem
from pathlib import Path
from write_dirs_files_for_gamess_qcxms import calculate_spin_multiplicity, get_method, get_props
from write_dirs_files_for_gamess_qcxms import read_file_rdkit, load_template, read_file, write_from_template, read_parameters
import pandas as pd
import argparse
import os


PBS_TEMPLATE = load_template("templates", "optimization_gamess_template.pbs")
INP_TEMPLATE = load_template("templates", "gamess_input_template.inp")
XYZ_TEMPLATE = load_template("templates", "structure_input_template.xyz")


def coords_as_dataframe(element_symbol, x_coord, y_coord, z_coord):
    """Convert coordinates DataFrame.

    Args:
        element_symbol (list): List of element symbols.
        x_coord (list): List of x-coordinates.
        y_coord (list): List of y-coordinates.
        z_coord (list): List of z-coordinates.

    Returns:
        Formatted DataFrame.
    """
    df = pd.DataFrame([element_symbol, x_coord, y_coord, z_coord]).transpose()
    df[0] = df[0].map("{0:<4s}".format)
    df[1] = pd.to_numeric(df[1], downcast="float").map("{0:>16.10f}".format)
    df[2] = pd.to_numeric(df[2], downcast="float").map("{0:>16.10f}".format)
    df[3] = pd.to_numeric(df[3], downcast="float").map("{0:>16.10f}".format)
    result = df[0] + df[1] + df[2] + df[3]

    return result


def read_coordinates(n_atoms, gamess_output, line_number):
    """Read coordinates from the GAMESS output.

    Args:
        n_atoms (int): Number of atoms in the molecule.
        gamess_output (list): Content of GAMESS output log file.
        line_number (int): Line number of the line currently being processed.

    Returns:
        tuple: Tuple containing lists of element symbols, x-coordinates, y-coordinates, and z-coordinates.
    """
    start, end = compute_start_end_indices(n_atoms, line_number, 4)

    data = [line.split() for line in gamess_output[start:end]]

    element_symbol, atomic_number, x_coord, y_coord, z_coord = zip(*data)

    return element_symbol, x_coord, y_coord, z_coord


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

    if not Path(mol_xyz_path).exists():
        with open(f"info_extract_coords_{project_dirname}.log", 'a') as f:
            f.write(f"Write XYZ file with optimized coordinates: {mol_xyz_path}")
            f.write(f"{os.linesep}")
        print(f"Write XYZ file with optimized coordinates: {Path(mol_xyz_path).parent}")
        write_from_template(parameters={"N_ATOMS": n_atoms, "MOLNAME": molname, "COORDINATES": result}, template=XYZ_TEMPLATE, file=mol_xyz_path)
        return 1
    return 0


def compute_start_end_indices(n_atoms, num, offset):
    """Compute start and end indices for coordinates in GAMESS output.

    Args:
        n_atoms (int): Number of atoms in the molecule.
        num (int): Line number.
        offset (int): Offset value.

    Returns:
        tuple: Tuple containing start and end indices.
    """
    start_coord_index = num + offset
    end_coord_index = start_coord_index + n_atoms
    
    return start_coord_index, end_coord_index


def read_coords_raw(gamess_output, start_coord_index, end_coord_index):
    """Read raw coordinates from GAMESS output.

    Args:
        gamess_output (list): Content of GAMESS output log file.
        start_coord_index (int): Start index for coordinates.
        end_coord_index (int): End index for coordinates.

    Returns:
        list: List of raw coordinates.
    """
    return [gamess_output[i] for i in range(start_coord_index, end_coord_index)]


def prepare_resubmission(params_filename, inchikey, n_atoms, molname, multiplicity, mol_dir, mol_input_path, message, gamess_output, num, offset):
    """Prepare resubmission by updating input files.

    Args:
        params_filename (str): Path to parameters file.
        inchikey (str): InChIKey identifier.
        n_atoms (int): Number of atoms in the molecule.
        molname (str): Identifier for the molecule.
        multiplicity (int): Spin multiplicity.
        mol_dir (Path): Path to the molecule directory.
        mol_input_path (Path): Path to the molecule input file.
        message (str): Resubmission message.
        gamess_output (list): Content of GAMESS output log file.
        num (int): Line number.
        offset (int): Offset value.

    Returns:
        int: 1 if resubmission is prepared, 0 otherwise.
    """
    start_coord_index, end_coord_index = compute_start_end_indices(n_atoms, num, offset)
    data = read_coords_raw(gamess_output, start_coord_index, end_coord_index)

    mol_pbs_path = mol_dir / (inchikey + ".pbs")
    pbs_params = read_parameters(params_filename, ["WALLTIME", "NCPUS", "MEM", "SCRATCH_LOCAL", "USER_EMAIL"])
    write_from_template({"MOLNAME": inchikey, **pbs_params}, PBS_TEMPLATE, mol_pbs_path)

    inp_params = read_parameters(params_filename, ["OPTTOL", "NSTEP", "GBASIS", "NGAUSS"])
    start_inp_file = mol_dir / (inchikey + "_start.inp")

    if not Path(start_inp_file).exists():
        Path(mol_input_path).rename(Path(start_inp_file))

    write_from_template({"SCFTYP": get_method(multiplicity), "MULT": multiplicity, **inp_params, "MOLNAME": molname, "COORDINATES": data, "text": message, "NPRT":-2}, INP_TEMPLATE, mol_input_path)

    return 1


def write_summary_log(project_dirname, count_abn, count_inp, count_xyz):
    """Write a summary log with counts of abnormal, INP, and XYZ molecules.

    Args:
        project_dirname (str): Name of the project directory.
        count_abn (int): Count of abnormal molecules.
        count_inp (int): Count of INP molecules.
        count_xyz (int): Count of XYZ molecules.
    """
    log_file_path = f"info_extract_coords_{project_dirname}.log"
    with open(log_file_path, 'a') as f:
        f.write(f"Number of abnormal molecules: {count_abn}\n")
        f.write(f"Number of INP molecules: {count_inp}\n")
        f.write(f"Number of XYZ molecules: {count_xyz}\n")
        f.write(f"Total of processed molecules= {count_xyz + count_inp + count_abn}\n")


listarg = argparse.ArgumentParser()
listarg.add_argument('--sdf_filename', type=str)
listarg.add_argument('--params_filename', type=str) 
listarg.add_argument('--project_dirname', type=str)
args = listarg.parse_args()


if __name__ == "__main__":
    mol_file = read_file_rdkit(args.sdf_filename)
    count_abn = 0
    count_inp = 0
    count_xyz = 0

    info_log = f"info_extract_coords_{args.project_dirname}.log"
    if Path(info_log).exists():
        os.remove(info_log)

    for mol in mol_file:
        mol = Chem.AddHs(mol)
        chem_class, inchikey, n_atoms, molname = get_props(mol)
        multiplicity = calculate_spin_multiplicity(mol)
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
