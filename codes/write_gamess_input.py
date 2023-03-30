#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from openbabel import pybel
import os


# In[ ]:


def sort_mols_into_classes(mols) -> dict:
    molecules_by_class = {}
    for mol in mols:
        if molecules_by_class.get(mol.data["Class"]) is None:
            molecules_by_class[mol.data["Class"]] = set()
    
        molecules_by_class[mol.data["Class"]].add(mol)
    return molecules_by_class


# In[ ]:


mols_from_sdf: list[pybel.Molecule] = list(pybel.readfile("sdf", "../data/sample.sdf"))


# In[ ]:


def get_method(x):
    if (x.OBMol.GetTotalSpinMultiplicity() % 2) == 0 :
        method = 'ROHF'
    else:
        method = 'RHF'
    return method


# In[ ]:


classified_mols = sort_mols_into_classes(mols_from_sdf)
if not os.path.exists("classes"):
    os.mkdir("classes")
    for chem_class in classified_mols.keys():
        outdir = os.path.join("classes", chem_class)
        os.mkdir(outdir)
        for mol in classified_mols[chem_class]:
            inchikey = mol.data["InChIKey"]
            molecule_dir = os.path.join(outdir, inchikey)
            os.mkdir(molecule_dir)
            with open(os.path.join(molecule_dir, inchikey + ".inp"), 'w') as outfile:
                mol.make3D()
                opt = f''' $CONTRL SCFTYP={get_method(mol)} MULT={mol.OBMol.GetTotalSpinMultiplicity()} RUNTYP=OPTIMIZE $END\n $STATPT OPTTOL=0.0005 NSTEP=100 $END\n $BASIS  GBASIS=N31 NGAUSS=6 $END'''
                outfile.write(mol.write("inp", opt={"k": opt}))

