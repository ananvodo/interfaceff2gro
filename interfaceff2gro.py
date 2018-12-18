#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 23:43:41 2018

@author: Andres Vodopivec
"""


import MDAnalysis
import sys
import subprocess
import os



def atom_param_from_prm(file):
    '''Extracts atomtypes from .prm charmm forcefield file. l is argument being a list containing the lines of forcefield file.'''
    #Getting the the atoms info index range
    for i in range(len(file)):
        file[i] = file[i].strip()
        if file[i][:9] == "NONBONDED":
            atom_upper_index = i
        if file[i][:5] == "HBOND":
            atom_lower_index = i
    
    for i in range(atom_upper_index, atom_lower_index):
        if file[i][:10] == "!INTERFACE":
            atom_upper_index_modified = i
    
    atoms_prm = []
    
    for i in range(atom_upper_index_modified, atom_lower_index):
        if len(file[i]) > 20 and file[i][:1] != '!':
                line = file[i][:40]
                line = line.strip()
                part = line.split()
                atoms_prm.append([str(part[0]), float(part[1]), float(part[2]), float(part[3])])
    
    return atoms_prm



def atoms_from_lammps(file):
    '''Extracts atomtypes from .data lammps file. l is argument being a list containing the lines of .data file.'''
#Getting the upper index and lower index to get the part I need from the file
    for i in range(len(file)):
        file[i] = file[i].strip()
        if file[i][:11] == "Pair Coeffs":
            atomtype_upper_index = i + 1  #the +1 is to get rid of the title lines
        if file[i][:11] == "Bond Coeffs":
            atomtype_lower_index = i - 1  #the -1 is to get rid of the title lines
   
    atoms_datafile = []     
    for i in range(atomtype_upper_index, atomtype_lower_index):
        if len(file[i]) > 20:
            line = file[i]
            line = line.strip()
            part = line.split()
            atoms_datafile.append([str(part[4]).upper()])
    
    return atoms_datafile



def bonds_from_prm(file):
    '''Extracts bondtypes from .prm charmm forcefield file. l is argument being a list containing the lines of forcefield file.'''
    #Getting the the bondtype info index range
    for i in range(len(file)):
        file[i] = file[i].strip()
        if file[i][:5] == "BONDS":
            bonds_upper_index = i
        if file[i][:6] == "ANGLES":
            bonds_lower_index = i

    for i in range(bonds_upper_index, bonds_lower_index):
        if file[i][:10] == "!INTERFACE":
            bonds_upper_index_modified = i
    
    #Getting the bonds info
    bonds_prm = []
    for i in range(bonds_upper_index_modified, bonds_lower_index):
        if len(file[i]) > 20 and file[i][:1] != '!':
                line = file[i]
                line = line.strip()
                part = line.split()
                bonds_prm.append([str(part[0]), str(part[1]), float(part[2]), float(part[3])])
    
    return bonds_prm



def bond_from_lammps(file):
    '''Extracts bondtypes from .data lammps file. l is argument being a list containing the lines of .data file.'''
    #Getting the the bondtype info index range
    for i in range(len(file)):
        file[i] = file[i].strip()
        if file[i][:11] == "Bond Coeffs":
            bond_upper_index = i + 2   #the +2 is to get rid of the empty lines
        if file[i][:12] == "Angle Coeffs":
            bond_lower_index = i - 1   #the -1 is to get rid of the empty lines
    #Getting the bondtype info
    bonds_datafile = []
    for i in range(bond_upper_index, bond_lower_index):
        line = file[i][48:]
        line = line.strip()
        part = line.split("-")
        bonds_datafile.append([str(part[0]).upper(), str(part[1]).upper()])
    
    return bonds_datafile



def angles_from_prm(file):
    '''Extracts angletypes from .prm charmm forcefield file. l is argument being a list containing the lines of forcefield file.'''
#Getting the the angle info index range
    for i in range(len(file)):
        file[i] = file[i].strip()
        if file[i][:6] == "ANGLES":
            angles_upper_index = i   
        if file[i][:9] == "DIHEDRALS":
            angles_lower_index = i 
    
    for i in range(angles_upper_index, angles_lower_index):
        if file[i][:10] == "!INTERFACE":
            angles_upper_index_modified = i
    
    #Getting the angles info
    angles_prm = []
    for i in range(angles_upper_index_modified, angles_lower_index):
        if len(file[i]) > 30 and file[i][:1] != '!':
                line = file[i][:35]
                line = line.strip()
                part = line.split()
                angles_prm.append([str(part[0]), str(part[1]), str(part[2]), float(part[3]), float(part[4])])
                
    return angles_prm



def angle_from_lammps(file):
    '''Extracts angletypes from .data lammps file. l is argument being a list containing the lines of .data file.'''
    #Getting the the atomtype info index range
    for i in range(len(file)):
        file[i] = file[i].strip()
        if file[i][:12] == "Angle Coeffs":
            angletype_upper_index = i + 2   #the +2 is to get rid of the empty lines
        if file[i][:15] == "Dihedral Coeffs":
            angletype_lower_index = i - 1   #the -1 is to get rid of the empty lines
    
    #Getting the angletype info  
    angletype_datafile = []
    for i in range(angletype_upper_index, angletype_lower_index):
        line = file[i][48:]
        line = line.strip()
        part = line.split("-")
        angletype_datafile.append([str(part[0]).upper(), str(part[1]).upper(), str(part[2]).upper()])
   
    return angletype_datafile



def create_lammps_data_file():
    '''Creates lampps *.data file from Materials studio .car .mdf and .frc file. To do this a utlity from Lammps package is used. J is the name of *.car file (e.g., something.car). 
    The car and mdf file needs to have the same name. The frc file need to be located in ./msi2lmp/frc_files. '''
    # The argument inserted in the terminal is ./interfaceff2gro.py mylammpsfile.data. So sys.argv[0] is /full_path/interfaceff2gro.py and sys.argv[1] is mylammpsfile.data
    bash_command = os.path.dirname(sys.argv[0]) + '/msi2lmp/src/msi2lmp.exe ' + sys.argv[1][:-4] + ' -class 2 -frc ' + os.path.dirname(sys.argv[0]) + '/msi2lmp/frc_files/pcff_interface.frc'
#    g = subprocess.Popen(bash_command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    g = subprocess.Popen(bash_command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(g.stdout.read())


#######################################################################################################################################################################
#######################################################################################################################################################################

#--------------------------------------------------------------------------------
# Creating the data file from msi2lmp.exe and loading the data file and prm file
#--------------------------------------------------------------------------------

    
# Create a data file from car and mdf files.
try:
    create_lammps_data_file()
except:
    print('First command-line argument is not a proper .car file, or .mdf file is missing/having different name than the .car file, or msi2lmp.exe is not located in ' + os.path.dirname(sys.argv[0]) + '/msi2lmp/src.')
    print('Usage: ./interfaceff2gro.py nameoffile.car.')
    sys.exit()

# Getting all the atom from the created data file.
try:
    u = MDAnalysis.Universe(sys.argv[1][:-4] + '.data')
except:
    print('The .data file from .car and .mdf files has not been created.')
    sys.exit()

# Open the forcefield file
try:
    forcefield = open(os.path.dirname(sys.argv[0]) + '/forcefield/charmm27_interface_v1_5.prm','r').readlines();
except:
    print('Frocefield .prm file missing. It should be located in ' + os.path.dirname(sys.argv[0]) + '/forcefield.')
    sys.exit()
    
# Read the created lammps data file
lammps = open(sys.argv[1][:-4] + '.data','r').readlines()


#######################################################################################################################################################################
#######################################################################################################################################################################

#----------------------------------------------------------------
# CREATING DICTIONARIES OF MASSES AND ATOMIC NUMBER OF ALL ATOMS
#----------------------------------------------------------------



# Dictionary for atomic numbers
#------------------------------

element2atomicNumber= {}
element2atomicNumber['H']='1'
element2atomicNumber['He']='2'
element2atomicNumber['C']='6'
element2atomicNumber['N']='7'
element2atomicNumber['O']='8'
element2atomicNumber['F']='9'
element2atomicNumber['Ne']='10'
element2atomicNumber['Na']='11'
element2atomicNumber['Mg']='12'
element2atomicNumber['Al']='13'
element2atomicNumber['P']='15'
element2atomicNumber['S']='16'
element2atomicNumber['Cl']='17'
element2atomicNumber['K']='19'
element2atomicNumber['Ca']='20'
element2atomicNumber['Fe']='26'
element2atomicNumber['Zn']='30'
element2atomicNumber['Br']='35'
element2atomicNumber['I']='53'
element2atomicNumber['Cs']='55'
element2atomicNumber['Si']='14'
element2atomicNumber['Ag']='47'
element2atomicNumber['Au']='79'
element2atomicNumber['Cu']='29'
element2atomicNumber['Ni']='28'
element2atomicNumber['Pb']='82'
element2atomicNumber['Pd']='46'
element2atomicNumber['Pt']='78'

# Dictionary for masses
#----------------------

mass2element= {}
mass2element[1.007970]='H'
mass2element[1.008000]='H'
mass2element[15.999400]='O'
mass2element[40.080000]='Ca'
mass2element[22.989770]='Na'
mass2element[26.981530]='Al'
mass2element[26.982000]='Al'
mass2element[28.086000]='Si'
mass2element[30.973800]='P'
mass2element[39.098300]='K'
mass2element[107.868000]='Ag'
mass2element[196.967000]='Au'
mass2element[63.546000]='Cu'
mass2element[58.710000]='Ni'
mass2element[207.200000]='Pb'
mass2element[106.400000]='Pd'
mass2element[195.090000]='Pt'

    
    

#######################################################################################################################################################################
#######################################################################################################################################################################
    
    
#------------------
# WRITING GRO FILE
#------------------
# For more information about gromacs format please refer to: http://manual.gromacs.org/current/online/gro.html
# For more information about lammps format please refer to: https://lammps.sandia.gov/doc/2001/data_format.html


atoms_for_gro = [] # Atom information for gro file (grabbed from LAMMPS)

for each_atom in zip(u.atoms.indices,u.atoms.masses,u.atoms.positions):
    atoms_for_gro.append(each_atom)

gro_file = open(sys.argv[1][:-4]+'.gro','w')

gro_file.write('Periodic slab: SURF, t= 0.0\n')
gro_file.write('%6d\n'%(len(atoms_for_gro)))

gro_format = '%5d%5s%5s%5d%8.3f%8.3f%8.3f\n' # Format of .itp file.
convert_dist = 10 # [Angstom / nm]

for each_atom in atoms_for_gro:
    resid_num = 1 # Only one residue, thus it is a constant
    resid_name = 'SURF'
    atom_num = each_atom[0] + 1 # Gromacs enumerates atoms according to index+1
    element = mass2element[each_atom[1]] # The atomic mass is the key to find the corresponding atom (or element) in the dictionary mass2element
    x_pos, y_pos, z_pos = each_atom[2][0]/convert_dist, each_atom[2][1]/convert_dist, each_atom[2][2]/convert_dist # all coordinates are in nm
    
    gro_file.write(gro_format%(resid_num, resid_name, element, atom_num, x_pos, y_pos, z_pos))

# Writing the box size
box_size_x, box_size_y, box_size_z = u.dimensions[0]/convert_dist, u.dimensions[1]/convert_dist, u.dimensions[2]/convert_dist # all dimensions are in nm
gro_file.write('%10.5f%10.5f%10.5f\n'%(box_size_x, box_size_y, box_size_z))
gro_file.close()




#######################################################################################################################################################################
#######################################################################################################################################################################

#------------------
# WRITING ITP FILE
#------------------

# For more information about gromacs format please refer to: http://www.pi.iccom.cnr.it/sites/default/files/joyce/joyce/node25.html
# For more information about lammps format please refer to: https://lammps.sandia.gov/doc/2001/data_format.html
# For more information about charmm prm format please refer to:http://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node25.html

itp_file = open(sys.argv[1][:-4]+'.itp','w')
#itp_file = open('test.itp','w')

kcal2kJ = 4.184 # conversion constant between kcal and kJ
    




    
# writing [atomtype] section
#---------------------------

atom_lammps = atoms_from_lammps(lammps) # Extract atoms present in lammps .data file from "Pair Coeffs" part. 
atom_lammps_dict = {} # Create a dict correlating index from .data file with typename
for i in range(len(atom_lammps)):
    atom_lammps_dict.update({i+1:atom_lammps[i][0]})


itp_file.write('[ atomtypes ]\n')
itp_file.write(';name	at.num	mass	charge	ptype	sigma	epsilon\n')


atomtypes_lmps = [] # Extracting the atomtypes info from "Atom" part in LAMMPS data file.
for m in set(zip(u.atoms.types, u.atoms.masses)): # zip creates a tuple, set gives an unordered collection with no duplicate elements
    atomtypes_lmps.append(m)



atom_prm = atom_param_from_prm(forcefield) # Extract all atoms, atom epsilion and sigma from .prm forcefield file.
atom_prm_Rmin_div_2, atom_prm_epsilion = {},{} # Create dicts having episilion and sigma values related to atomtype
for atom in atom_prm:
    atom_prm_Rmin_div_2[atom[0]] = atom[3] # [Angstrom]
    atom_prm_epsilion[atom[0]] = atom[2] # [kcal/mol]



for each_atom in atomtypes_lmps:
    atom_name = atom_lammps_dict[int(each_atom[0])]
    atomic_num = element2atomicNumber[mass2element[each_atom[1]]]
    atom_mass = each_atom[1]
    atom_charge = 0
    particle_type = "A"
    atom_sigma = atom_prm_Rmin_div_2[atom_name] * (2 / convert_dist)/ (2**(1/6)) # [nm]  Dividing by 2^(1/6) because of difference in CHARMM and GMX LJ formula
    atom_epsilion = -(kcal2kJ * atom_prm_epsilion[atom_name]) # [kJ/mol]   Epsilion is positive in gmx format
    
    if atom_mass > 10:
        cformat_atomtypes = '%-s\t%-s\t%8.5f\t%-.6f\t%s\t%.12f\t%.6f \n'
        itp_file.write(cformat_atomtypes%(atom_name, atomic_num, atom_mass, atom_charge, particle_type, atom_sigma, atom_epsilion))
    
    else:
        cformat_atomtypes_2 = '%-s\t%-s\t%8.6f\t%-.6f\t%s\t%.12f\t%.6f \n'
        itp_file.write(cformat_atomtypes_2%(atom_name, atomic_num, atom_mass, atom_charge, particle_type, atom_sigma, atom_epsilion))







# writing [bondtype] section
#---------------------------

bonds_using_atom_name_lmps = bond_from_lammps(lammps) # Doing a list of all bonds using atom name present in LAMMPS file.

bonds_prm = bonds_from_prm(forcefield) # Extract all bonds and bonding parameters from .prm forcefield file. 

if len(bonds_using_atom_name_lmps) > 0:
    itp_file.write('\n')
    itp_file.write('[ bondtypes ]\n')
    itp_file.write('; i	j	func	b0	kb\n')
    cformat_bondtypes = '%-s\t%-s\t%d\t%-.4f\t%-.2f\n'
    
    for bond_lammps in bonds_using_atom_name_lmps:
        
        for bond_prm in bonds_prm:
            
            if (bond_lammps[0] == bond_prm[0] and bond_lammps[1] == bond_prm[1]) or (bond_lammps[1] == bond_prm[0] and bond_lammps[0] == bond_prm[1]):
                atom_1 = bond_prm[0]
                atom_2 = bond_prm[1]
                function_type = 1 # all bonds are harmonic form
                equ_dist = bond_prm[3] / convert_dist # [nm]
                force_const = bond_prm[2] * 2 * (convert_dist**2) * kcal2kJ # converting [kcal/mol/A**2] to [kJ/mol/nm**2]. Include factor 2 (see definitions in Gromacs manual)
                
                itp_file.write(cformat_bondtypes%(atom_1, atom_2, function_type, equ_dist, force_const))






# writing [angletype] section
#----------------------------

angles_using_atom_name_lmps = angle_from_lammps(lammps) # Extract angles from .data file

angles_prm = angles_from_prm(forcefield) # Extract all angles and angle parameters from .prm forcefield file.


if len(angles_using_atom_name_lmps) > 0:
    itp_file.write('\n')
    itp_file.write('[ angletypes ]\n')
    itp_file.write('; i	j	k	func	th0	cth\n')
    cformat_angletypes = '%-s\t%-s\t%-s\t%d\t%-.2f\t%-.4f\t%.1f\t%.1f\n'
    
    for angle_lmps in angles_using_atom_name_lmps:
        
        for angle_prm in angles_prm:
            atom_1 = angle_prm[0]
            atom_2 = angle_prm[1]
            atom_3 = angle_prm[2]
            function_type = 5
            equ_angle = angle_prm[4]
            force_const = angle_prm[3] * 2 * kcal2kJ # Include factor 2 (see definitions in Gromacs manual)
            
            if (angle_lmps[0] == angle_prm[0] and angle_lmps[1] == angle_prm[1]and angle_lmps[2] == angle_prm[2]) or (angle_lmps[2] == angle_prm[0] and angle_lmps[0] == angle_prm[2] and angle_lmps[1] == angle_prm[1]):
                itp_file.write(cformat_angletypes%(atom_1, atom_2, atom_3, function_type, equ_angle, force_const, 0, 0))





# There are no Dihedral nor improper dihedrals 
#---------------------------------------------




# writing [molecule] section
#---------------------------

itp_file.write('\n')
itp_file.write('[ moleculetype ]\n')
itp_file.write('; molname    nrexcl\n')
nrexcl = 3
cformat_molecule = '%-5s%4d\n'
itp_file.write(cformat_molecule%('SURF', nrexcl))






# writing [atom] section
#-----------------------

itp_file.write('\n')
itp_file.write('[ atoms ]\n')
itp_file.write(';   nr       type  resnr residue  atom   cgnr     charge       mass\n')
cformat_atoms = '%6d%11s%7d%7s%7s%7d%11.3f%11.3f   ;\n'


atoms_lammps = [] # Doing a list of all atoms in LAMMPS file.
for each_atom in zip(u.atoms.indices, u.atoms.types, u.atoms.masses, u.atoms.charges):
    atoms_lammps.append(each_atom)


for atom_lmps in atoms_lammps:
    atom_num = atom_lmps[0] + 1
    atom_type = atom_lammps_dict[int(atom_lmps[1])]
    resid_num = 1 # There is only one residue
    resid_name = 'SURF'
    element = mass2element[atom_lmps[2]] # reference name according to periodic table
    cgnr = atom_lmps[0] + 1 # Charge group number
    atom_carge = atom_lmps[3]
    atom_mass = atom_lmps[2]
    
    itp_file.write(cformat_atoms%(atom_num, atom_type, resid_num, resid_name, element, cgnr, atom_carge, atom_mass))






# writing [bond] section
#-----------------------

bonds_using_atom_index_lmps = [] # Doing a list of all bonds using atom index number present in LAMMPS file.
try:
    for each_bond in u.bonds.indices:
        bonds_using_atom_index_lmps.append(each_bond)
except:
    print('Bonds between atoms have not been found.')



cformat_bonds = '%5d%6d%6d \n'

if len(bonds_using_atom_index_lmps) > 0:
    itp_file.write('\n')
    itp_file.write('[ bonds ]\n')
    itp_file.write(';  ai    aj funct            c0            c1            c2            c3\n')
    
    for bond_lmps in bonds_using_atom_index_lmps:
        atom_1 = bond_lmps[0] + 1
        atom_2 = bond_lmps[1] + 1
        function_type = 1
        
        itp_file.write(cformat_bonds%(atom_1, atom_2, function_type))
  







# writing [angle] section
#------------------------

cformat_angles = '%5d%6d%6d%6d \n'


angles_using_atom_index_lmps = [] # Doing a list of all angles in LAMMPS file.
try:
    for each_angle in u.angles.indices:
        angles_using_atom_index_lmps.append(each_angle)
except:
    print('Angles between atoms have not been found.')


if len(angles_using_atom_index_lmps) > 0:
    itp_file.write('\n')
    itp_file.write('[ angles ]\n')
    itp_file.write(';  ai    aj    ak funct            c0            c1            c2            c3\n')
    
    for angle_lmps in angles_using_atom_index_lmps:
        atom_1 = angle_lmps[0] + 1
        atom_2 = angle_lmps[1] + 1
        atom_3 = angle_lmps[2] + 1
        function_type = 5
        itp_file.write(cformat_angles%(atom_1, atom_2, atom_3, function_type))







# Wriritng position restraints
#-----------------------------
        
itp_file.write('\n')
itp_file.write('; Include Position restraint file\n')
itp_file.write('#ifdef POSRES\n')
itp_file.write('#include "posre_SURF.itp"\n')
itp_file.write('#endif\n')
itp_file.write('\n')







# Closing file
#-------------

itp_file.close()




#######################################################################################################################################################################
#######################################################################################################################################################################


#------------------
# WRITING TOP FILE
#------------------


top_file = open('system.top','w')
top_file.write('#include "charmm27.ff/forcefield.itp"\n')
top_file.write('#include "%s.itp"\n'%sys.argv[1][:-4])
top_file.write('#include "charmm27.ff/ions.itp"\n')
top_file.write('#include "charmm27.ff/tip4p.itp"\n')
top_file.write('\n')
top_file.write('[ system ]\n')
top_file.write('Surface\n')
top_file.write('\n')
top_file.write('[ molecules ]\n')
top_file.write('SURF 1 ;1 periodic slab\n')
top_file.write('\n')
top_file.close()




#######################################################################################################################################################################
#######################################################################################################################################################################


# Printing output to terminal
print('Three files have been created: %s.itp, %s.gro and system.top.'%(sys.argv[1][:-4], sys.argv[1][:-4]))
