#!/usr/bin/env python
# coding: utf-8

# In[1]:


from jinja2 import Environment, FileSystemLoader
import numpy as np
import pandas as pd 


# In[2]:


opt_f='bromide.log'
not_opt_f='coumarin.log'
with open(not_opt_f, 'r') as f:
    lines = f.readlines()


# In[3]:


# Check if optimization was succesfull 
opt = False
for line in lines:
    if  "**** THE GEOMETRY SEARCH IS NOT CONVERGED! ****" in line:
        opt = True


# In[4]:


# If optimization is successfull 
if not opt:
    data = [] 
    for num, line in enumerate(lines):
        if "TOTAL NUMBER OF ATOMS" in line:
            n_atoms = int(line.split()[5])
        if "EQUILIBRIUM" in line:
            start_coord_index = num + 4
            end_coord_index = start_coord_index + n_atoms
            for i in range(start_coord_index, end_coord_index):
                data.append(lines[i].split())
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


# In[5]:


file_loader=FileSystemLoader("../templates")
environment = Environment(loader=file_loader)
if not opt:
    template = environment.get_template("structure_input_template.xyz")

    MOLNAME="bromide"
    content = template.render(N_ATOMS=n_atoms, MOLNAME=MOLNAME, COORDINATES=result)
    print(content)
    
    # write into a file
    # filename = MOLNAME+'.xyz'
    # with open(filename, mode="w", encoding="utf-8") as message:
    #    message.write(content)


# In[6]:


# If optimization is NOT successfull 

if  opt:
	template = environment.get_template("gamess_input_template.inp")
	data = []
	for num, line in enumerate(lines):
		if "TOTAL NUMBER OF ATOMS" in line:
			n_atoms = int(line.split()[5])
		if "ALWAYS THE LAST POINT COMPUTED!" in line:
			start_coord_index = num + 4
			end_coord_index = start_coord_index + n_atoms
			for i in range(start_coord_index, end_coord_index):
				data.append(lines[i])
				
	with open("gamess_parameters.in", 'r') as f:
		lines = f.readlines()
	for line in lines:
		if "SCFTYP" in line:
			SCFTYP = line.split()[2]
		if "MULT" in line:
			MULT = line.split()[2]
		if "OPTTOL" in line:
			OPTTOL = line.split()[2]
		if "NSTEP" in line:
			NSTEP = line.split()[2]
		if "GBASIS" in line:
			GBASIS = line.split()[2]
		if "NGAUSS" in line:
			NGAUSS = line.split()[2]
		if "MOLNAME" in line:
			MOLNAME = line.split()[2]


	content = template.render(SCFTYP=SCFTYP, MULT=MULT, OPTTOL=OPTTOL,
			NSTEP=NSTEP, GBASIS=GBASIS, NGAUSS=NGAUSS, MOLNAME=MOLNAME, COORDINATES=data)
	print(content)
	
	# write into a file
	# filename = MOLNAME+'.inp'
	# with open(filename, mode="w", encoding="utf-8") as message:
	# 	message.write(content)


# In[7]:


with open("cluster_parameters.in", 'r') as f:
    lines = f.readlines()
for line in lines:
    if "WALLTIME" in line:
        WALLTIME = line.split()[2]
    if "BIN" in line:
        BIN = line.split()[2]
    if "USER_EMAIL" in line:
        USER_EMAIL = line.split()[2]

template = environment.get_template("array_sub_template.sh")

content = template.render(WALLTIME=WALLTIME, BIN=BIN, USER_EMAIL=USER_EMAIL)
print(content)

