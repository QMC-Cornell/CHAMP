#!/usr/bin/env python
import sys
import re
from string import *
import os
##############################################################################
program_name = "xyz_2_ab_initio_input.py"
## 
##	Requires XYZ file in format: ATOM CHARGE X Y Z
##		 
##		 
##	Written by DOMINIK DOMIN (dominik.domin@gmail.com)
##	Date of last modification February 15, 2011
version = 1.02
################################################################################
##  CONVERTS XYZ FILES INTO INPUT FOR GAUSSIAN OR GAMESS OR NWCHEM OR MOLPRO  ##
################################################################################

from abinitio_input_write import *

if version_number < 1 and subversion_number>=0:
	print "******************************************************************************************"
	print "DEPENDENCY ERROR: abinitio_input_write.py is out of date and incompatible with this script"
	print "******************************************************************************************"
	sys.exit(-5)

#default D3_G_BASIS + D3_G_ECP
st = sys.argv[0]
ri = str(os.path.basename(st))
ng = st.index(ri)
st = str(st[0:ng])
##
current_path = os.getenv('PWD')
basis_path = os.getenv('D3_G_BASIS')
ecp_path = os.getenv('D3_G_ECP')
if basis_path == None:
	basis_path = st+"BFD_GAUSSIAN_BASIS_AND_ECP/"
	if (str(os.path.exists(basis_path)) == "False"):
		basis_path = current_path
if ecp_path == None:
	ecp_path = st+"BFD_GAUSSIAN_BASIS_AND_ECP/"
	if (str(os.path.exists(ecp_path)) == "False"):
		ecp_path = current_path

## HELP MENU
if (len(sys.argv)>1):
	if (sys.argv[1]=="h" or sys.argv[1]=="help" or sys.argv[1]=="-h" or sys.argv[1]=="-H" or sys.argv[1]=="H"):
		print " ------------------HELP------------------------------------------------------- "
		print "CORRECT USAGE: %s xyz_file units Charge Multiplicity Program" % program_name
		print "  Where xyz_file is the name of file containing xyz data "
		print "     in ATOM CHARGE X Y Z format"
		print "    AND "
		print "  Units can be 'A' for Angstrom or 'B' for Bohr"
		print "  Charge must be an integer  "
		print "  Multiplicity (2S+1) must be an integer "
		print "  Program must be GAUSSIAN, GAMESS, NWCHEM, or MOLPRO"
		print " ------------------------------------------------------------------------------ "
		print " NOTE: This script is good for single determinant functions."
		print "       You must also have in the same directory as your xyz_file the GAUSSIAN atomic basis files named just like your atoms"
		print "           If you do not provide the GAUSSIAN atomic basis files you will get the following error message:"
		print "                                                                         \"IOError: [Errno 2] No such file or directory\"\n"
		print " If you supply GAUSSIAN ECP files then the calculation will be a pseudopotential calculation"
		print " PLEASE BE SURE TO SET D3_G_BASIS AND D3_G_ECP IN YOUR ENVIRONMENTAL VARIABLES"
		print "												*** SCRIPT VERSION %1.2f ***" % version
		sys.exit(0)
if (len(sys.argv)<6):
	print "SYNTAX ERROR: too few arguments"
	print "-------------------------------------------------------------------------------"
	print "CORRECT USAGE: %s  xyz_file units Charge Multiplicity Program" % program_name
	print "  *** FOR FURTHER HELP type: %s help ***" % program_name
	print "_______________________________________________________________________________"
	sys.exit(0)


## UNITS of the xyz file are read in at one sets unitfactor variable

### NOTE: Angstrom units are used to store geometries. Please not this is completely opposite of my previous scripts! 
### Warning if mixing and matching from other scripts please keep in mind that the geometries here are in ANGSTROM AND NOT IN BOHR units!

unitfactor = 0.0
if sys.argv[2] == "B" or sys.argv[2] =="b":
	unitfactor = 1.0 * 0.529177
if sys.argv[2] == "A" or sys.argv[2] =="a":
	unitfactor = 1.0 
if sys.argv[2] != "A":
	if sys.argv[2] != "a":
		if sys.argv[2] != "B":
			if sys.argv[2] != "b":	
				print "ERROR: Units must be specified 'A' (Angstrom) or 'B' (Bohr)"
				sys.exit(0)
if unitfactor == 0.0:
	print "ERROR unitfactor = 0"
	
print " Welcome to the XYZ file -> ab initio input file GENERATOR version %1.2f\n X Y Z coordinates will be multiplied by factor = %f " % (version,unitfactor)

## PERIODIC TABLE
element_list = ['X', 'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al', 'Si','P','S','Cl','Ar','K','Ca','Sc', 'Ti','V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Uut','Uuq','Uup','Uuh','Uus','Uuo'] #Note that 'X' is just a dummy atom useful for Basis Set Superposition Error and floating basis sets.

## VARIABLES
charge = int(atof(str(sys.argv[3])))
ename=[] #ELEMENT (atom) NAME
zcharges =[] #Charge of atom
nmul = int(atof(sys.argv[4])) #Spin Multiplicity
nelectrons = -charge #Number of Electrons
temptemptemp = str(sys.argv[5])
program_format = lower(temptemptemp)

if program_format[0:3] =="gau":
	print "Gaussian format has been selected"
	program_format = "gaussian"
	basis_path = str(basis_path)+"/G09_BASIS/"
	ecp_path = str(ecp_path)+"/G09_ECP/"
if program_format[0:3]=="gam":
 	print "GAMESS format has been selected"
	program_format = "gamess"
	basis_path = str(basis_path)+"/GAMESS_BASIS/"
	ecp_path = str(ecp_path)+"/GAMESS_ECP/"
if program_format[0:2]=="nw":
	program_format = "nwchem"
	basis_path = str(basis_path)+"/NWCHEM_BASIS/"
	ecp_path = str(ecp_path)+"/NWCHEM_ECP/"
if program_format[0:2]=="mo":
	print "MOLPRO format has been selected"
	basis_path = str(basis_path)+"/MOLPRO_BASIS/"
	ecp_path = str(ecp_path)+"/MOLPRO_ECP/"
	program_format = "molpro"
else:
	if program_format[0:3]!="mol" and program_format[0:3]!="nwc" and program_format[0:3]!="gam" and program_format[0:3]!="gau":
		print "ERROR: unrecognized software format  %s" % str(program_format)
		print "     Formats recognized: GAUSSIAN, NWCHEM, GAMESS, and MOLPRO" 
		sys.exit(-1)

xcor=[] # X coordinates in Bohr
ycor=[] # Y coordinates in Bohr
zcor=[] # Z coordinates in Bohr
gbas = [] #temp gamess basis list
gpseudo = [] #temp gamess basis list
name1 = str(sys.argv[1])
name2 = name1.replace("cc-pv","")
tempname2 = name2.replace(".xyz","")
name2 = tempname2.replace(".","")
nalpha = 0 #number of alpha electrons
nbeta = 0 # number of beta electrons
print " Molecule has Charge of %d and Spin Multiplicity of %d \n" % (-nelectrons,nmul)

## Reading in the xyz file
lines = open(sys.argv[1], 'r').readlines()
for line in lines:
	splitoid = split(line)
	if len(splitoid)!=0:
		if len(splitoid)==4:
			for a in range(len(element_list)):
				if splitoid[0]==element_list[a]:
					ename.append(splitoid[0])
					zcharges.append(a)
					nelectrons = nelectrons + a
					xcor.append((atof(str(splitoid[1]))*unitfactor))
					ycor.append((atof(str(splitoid[2]))*unitfactor))
					zcor.append((atof(str(splitoid[3]))*unitfactor))
					break
				else:
					if a == len(element_list)-1:
						print "****************************************"
						print "*** ERROR: Could not find element %s ***" % splitoid[0]
						print "****************************************"
						print " Perhaps the name of the element is not appropriately capitalized?"
						sys.exit(-1)
		else:
			if len(splitoid)<5:
				print("**ERROR IN XYZ FILE**\n TOO FEW ARGUMENTS\n")
				print str(splitoid)
				print("---------------------------------------------")
				print(" XYZ FILE must be in format:   ELEMENT_NAME  CHARGE  X  Y  Z ")
				sys.exit(-1)
			ename.append(splitoid[0])
			zcharges.append(int(atof(str(splitoid[1]))))
			nelectrons = nelectrons + int(atof(str(splitoid[1])))
			xcor.append((atof(str(splitoid[2]))*unitfactor))
			ycor.append((atof(str(splitoid[3]))*unitfactor))
			zcor.append((atof(str(splitoid[4]))*unitfactor))

nalpha = (nmul - 1 + nelectrons)/2
nbeta = nelectrons - nalpha
unique_pseudo_list_ename = []
unique_list_ename = []
for i in range(len(ename)):
	if unique_list_ename.count(str(ename[i]))==0:
		unique_list_ename.append(str(ename[i]))
		try:
			open(str(ename[i])+".ecp",'r')
			unique_pseudo_list_ename.append(str(ename[i])+".ecp")
			print "    Found ecp for %s in %s" % (str(ename[i]),current_path)
		except IOError:
			try:
				temp_name = str(ecp_path)+"/"+str(ename[i])+".ecp"
				temp2 = temp_name.replace("//","/")
				temp_name = temp2.replace("//","/")
				open(temp_name,'r')
				unique_pseudo_list_ename.append(str(temp_name))
				print "    Found ecp for %s in %s" % (str(ename[i]),ecp_path)
			except IOError:
				print "Could not find ECP for %s" % str(ename[i])
				print "Warning: Assuming no ECP for %s" % str(ename[i])
unique_atom_list_ename = []
for j in unique_list_ename:
	try:
		open(str(j),'r')
		unique_atom_list_ename.append(str(current_path+"/"+j))
		print "    Found atom %s basis in %s" % (str(j),str(current_path))
	except IOError:
		try:
			temp_name = str(basis_path)+"/"+str(j)
			temp2 = temp_name.replace("//","/")
			temp_name = temp2.replace("//","/")
			open(str(temp_name),'r')
			print "    Found atom %s basis in %s" % (str(j),str(basis_path))
			unique_atom_list_ename.append(str(temp_name))
		except IOError:
			print "ERROR: Could not find basis set for atom %s in neither %s nor %s" % (str(j),str(current_path),str(basis_path))
			print "*** ABORTING %s" % program_name
			sys.exit(-1)

pseudo_yes = "N"
pseudo = []
if len(unique_pseudo_list_ename)>0:
	pseudo_yes="T"


if program_format == "gaussian":
	write_gaussian(name1,name2,nmul,pseudo_yes,charge,ename,unique_atom_list_ename,unique_pseudo_list_ename,pseudo,xcor,ycor,zcor)
if program_format == "gamess":
	write_gamess(name1,name2,nmul,pseudo_yes,charge,ename,unique_atom_list_ename,unique_pseudo_list_ename,pseudo,xcor,ycor,zcor,zcharges,ecp_path,basis_path)
if program_format == "nwchem":
	mem = 1800
	write_nwchem(mem,name1,name2,nmul,pseudo_yes,charge,ename,unique_atom_list_ename,unique_pseudo_list_ename,pseudo,xcor,ycor,zcor,zcharges,ecp_path,basis_path)
if program_format == "molpro":
	mem=200
	write_molpro(mem,name1,name2,nmul,pseudo_yes,charge,ename,unique_atom_list_ename,unique_pseudo_list_ename,pseudo,xcor,ycor,zcor,zcharges,ecp_path,basis_path)




