from openbabel import openbabel, pybel

def read_file(format, inputfilename):
    mol = next(pybel.readfile(format, inputfilename))
    return (mol)

def convert(mol):
    mol.make3D()
    #mol.OBMol.AddHydrogens()
    return(mol)

def write_file(mol, format, outfilename):
    return mol.write(format, outfilename, overwrite=True)

if __name__ == "__main__":
    format_in = "smi"
    format_out = "xyz"
    inputfilename = "sample.smi"
    outfilename = "sample_conv.xyz"
    mol = read_file(format_in, inputfilename)
    mol = convert(mol)
    write_file(mol, format_out, outfilename)



