from openbabel import pybel

mol = next(pybel.readfile("smi", "sample.smi"))
mol.make3D()
mol.OBMol.AddHydrogens()

mol.write("xyz", "sample_converted.xyz", overwrite=True)