from openbabel import pybel

def read_file(format, filename):
    pybel.ob.obErrorLog.SetOutputLevel(0)
    mols = list(pybel.readfile(format, filename))
    return (mols)

def convert(format_in, inputfilename, format_out, outfilename):
    mols = read_file(format_in, inputfilename)
    out = pybel.Outputfile(format_out, outfilename+"."+format_out, overwrite=True)
    for line in mols:
        line.make3D()
        out.write(line)
    #print(mol.OBMol.NumAtoms())
    

if __name__ == "__main__":
    format_in = "sdf"
    format_out = "xyz"
    file = "sample."+ format_in
    convert(format_in, file, format_out, "sample_out")


