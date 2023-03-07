from openbabel import pybel
from pathlib import Path

def read_file(format, filename):
    pybel.ob.obErrorLog.SetOutputLevel(0)
    mols = list(pybel.readfile(format, filename))
    return (mols)

def convert(format_in, inputfilename, format_out, outfilename):
    mols = read_file(format_in, inputfilename)
    out = pybel.Outputfile(format_out, outfilename, overwrite=True)
    for line in mols:
        line.make3D()
        out.write(line)    

if __name__ == "__main__":
    file = Path("sample.sdf")
    outfile = file.with_suffix('.xyz')
    format_in = file.suffix[1:]
    format_out = outfile.suffix[1:]

    convert(format_in, file.name, format_out, outfile.name)


