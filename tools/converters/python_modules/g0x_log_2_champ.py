#!/usr/bin/python

#########################################################################################################################################
program_name = "g0x_log_2_champ.py"																										#
version_num = 0.04																														#
#																																		#
#  This script extracts orbitals from Gaussian 03 log files																				#
#  and it will create a champ output file																								#
#																																		#
arguments_list = "gaussian_03_output.log (N/M)"																								#
usage_description = "This script generates CHAMP input files based Gaussian 03/09 log files. \nIt just generates Hartree-Fock.\n An optional second argument for chosing whether purely numerical (N) or mixed numerical-analytical basis set (M) will be printed.\n"
#
#  Author: Dominik Domin	
#  Last Modified on  October 21, 2010
#########################################################################################################################################

import sys
import re
from string import *
from orbital_tools_v1 import *
from general_tools_v1 import *
from gaussian_read_v1 import *
from pseudovb_input_read_v1 import *
from write_champ_v1 import *
import os.path


############################################################################################################################################
##### MAIN DRIVER FOR SCRIPT ######                                                                                                 
############################################################################################################################################
#start of main program

print "** %s version %1.2lf **" % (program_name,version_num)	
if len(sys.argv)>1:
	if sys.argv[1]=="help" or sys.argv[1]=="-h" or sys.argv[1]=="--h" or sys.argv[1]=="HELP" or sys.argv[1]=="-H" or sys.argv[1]=="-help":
		script_help(program_name, arguments_list, usage_description)

if len(sys.argv)<2:
	print "ERROR: insufficiant arguments"
	print "Usage: %s %s" % (program_name,arguments_list)
	print "If you need further help type -h"
	sys.exit(0)

#test that all modules are the required versions ####################################################################
#     axx[1] contains version number and axx[3] contains the revision number
awc = split(write_champ_version_string)
agr = split(gaussian_read_version_string)
apir = split(pseudovb_input_read_version_string)
agt = split(general_tools_version_string)
aot = split(orbital_tools_version_string)
a=0 
if str(awc[1])!=str(1) or int(awc[3])< int(5):
	print "FATAL ERROR: write_champ module is out of date or incompatible with this script"
	a = 1
if str(agr[1])!=str(1) or int(agr[3])< int(10):
	print str(agr)
	print "FATAL ERROR: gaussian_read module is out of date or incompatible with this script"
	a = 1
if str(agt[1])!=str(1) or int(agt[3])< int(5):
	print "FATAL ERROR: general_tools module is out of date or incompatible with this script"
	a = 1
if str(aot[1])!=str(1) or int(aot[3])< int(2):
	print "FATAL ERROR: orbital_tools module is out of date or incompatible with this script"
	a = 1
if a == 1:
	sys.exit(0)
#end of testing of module versions ###############################################################################	

base_file_name = sys.argv[1]
numericalORmixed = "M"
if len(sys.argv)==3:
	numericalORmixed = str(sys.argv[2]).upper()
	if numericalORmixed == "N":
		print "Purely numerical/grid basis selected"
	if numericalORmixed =="M":
		print "Mixed numerical and analytical basis selected"
	if numericalORmixed !="N" and numericalORmixed !="M":
		print "ERROR: %s is not a valid basis type selection"
		print "     Please choose between N for numerical/grid"
		print "     and M for mixed numerical and analytical"
		print "     basis sets."
		sys.exit(0) 
else:
	print "Warning: no basis set type choosen. Assuming Mixed basis type"
loglines=[]
index_start_stop = []
everything = []

if base_file_name[-3:]=="log":
	loglines = open(base_file_name,"r").readlines()
else:
	loglines = open(base_file_name+".log","r").readlines()

gv=check_log_gaussian_version(loglines)
if gv==0:
	print "ERROR: %s is most likely not a GAUSSIAN(R) log file." % base_file_name
	sys.exit(0)
if gv < 10000:	
	print "Parsing Gaussian version 0%d" % int(gv/1000)
else:
	print "Parsing Gaussian version	%d" % int(gv/1000)
#reading in geometry

restrained_return_start_stop_index("-- Stationary point found.","Standard orientation:","Rotational constants",index_start_stop,loglines)
a = int(index_start_stop[0])+4
b = int(index_start_stop[1])-2
if a!=4:
	index_start_stop=[]
	index_start_stop.append(a)
	index_start_stop.append(b)
else:
	index_start_stop=[]
	return_start_stop_index("Input orientation:","Distance matrix", index_start_stop,loglines)
	index_start_stop[0]=index_start_stop[0]+4
	index_start_stop[1]=index_start_stop[1]-2
	if index_start_stop[0]==4:
		index_start_stop=[]	
		return_start_stop_index("Z-Matrix orientation","Stoichiometry", index_start_stop,loglines)
		index_start_stop[0]=index_start_stop[0]+4
		index_start_stop[1]=index_start_stop[1]-3

molecular_geometry = read_geometry(index_start_stop,loglines)
index_start_stop=[]

#reading in parameters 
ae=0
be=0
te=0.0
temp = []
temp = read_energy_number_electrons(loglines)
ae = temp[0]
be = temp[1]
te = temp[2]

#reading in SCF orbitals
rsss = []
return_start_stop_index("Molecular Orbital Coefficients","DENSITY MATRIX", index_start_stop,loglines)
if index_start_stop[0]!= 0:
	if index_start_stop[1]==0:
		index_start_stop = []
		return_start_stop_index("Molecular Orbital Coefficients","Density Matrix:", index_start_stop,loglines)
		rsss = read_orbitals_scf("Eigenvalues --","\n",index_start_stop,loglines)
	else:
		rsss = read_orbitals_scf("EIGENVALUES --","\n",index_start_stop,loglines)

scf_orbitals = []
occupation_numbers = []
bf_labels = []
energies = []
sym_labels = []
if rsss!=None and len(rsss)!=0:
        scf_orbitals = rsss[0]
        occupation_numbers = rsss[1]
	bf_labels = rsss[2]
	energies = rsss[3]
	sym_labels = rsss[4]
index_start_stop=[]

#reading in Basis Set
return_start_stop_index("AO basis set in the form of general basis input","^\n",  index_start_stop,loglines)
index_start_stop[0]=index_start_stop[0]-1
atoms = read_basis_set(index_start_stop,loglines)
unique_atoms = return_unique_list(atoms)
#reading in Pseudopotential
index_start_stop = []
restrained_return_start_stop_index("Pseudopotential Parameters","Number     Electrons     Momentum     of R      Exponent        Coefficient","Integral buffers will be",index_start_stop,loglines)
if index_start_stop[0]!=0:
	ecps = read_pseudopotentials(index_start_stop,loglines)
else:
	ecps=[]
ase= []
ap = []
for ecp in ecps:
	s = split(str(ecp[-1]))
	se = int(s[-2])
	ase.append(se)
	
unique_ecps = return_unique_list(ecps)
for ecp in unique_ecps:
	s = split(str(ecp[-1]))
	p = int(s[-1])+1
	ap.append(p)

i = 0 
for ecp in unique_ecps:
	i = i + 1
	write_champ_pseudopotentials("pseudopot%d" % (i),ecp)

#create csfs
#creates a simple HF determinant
determinants = ""
if len(occupation_numbers)==0:
	index_start_stop=[]
	return_start_stop_index("Occupied","Virtual",  index_start_stop,loglines)
	occupation_numbers=parse_occupation_guess(index_start_stop,loglines)

for a in range(ae):
	if occupation_numbers[i]=="O":
		determinants = determinants + "%d " % (a+1)
determinants = determinants + "   "
for b in range(be):
	if occupation_numbers[i]=="O":
		determinants = determinants + "%d " % (b+1)
determinants = determinants + " \n "
vb_weights = "true"
dets_in_csfs = "1\n1.0\n"
csf_coef=" 1.000 "
csfs=[]
csfs.append(determinants)
csfs.append(csf_coef)
csfs.append(dets_in_csfs)
csfs.append(vb_weights)



##########################################################################
# CREATE CHAMP INPUT FILE
##########################################################################
if base_file_name[-4:]==".log":
	champfile=str(base_file_name[:-4])+"-champ.inp"
else:	
	champfile=str(base_file_name[:base_file_name.index(".")-1])+"-champ.inp"

if (os.path.exists(champfile)):
	print "WARNING overwriting the file %s which already exists" % (champfile)
print "champ file %s" % (champfile)
champ_file=open(champfile,"w")
champ_file.write("! %s version %s - %s%s%s%s%s%s%s%s%s%s\n\n" % (program_name,version_num,awc[1],awc[3],agr[1],agr[3],apir[1],apir[3],agt[1],agt[3],aot[1],aot[3]))
champ_file.close()
control_group =None
if control_group ==None:
	dlft=champ_default()
	mc = dlft.mc
	seed = dlft.seed
	etrial = dlft.etrial
	nstep =dlft.nstep
	nblk=dlft.nblk
	nblkeq=dlft.nblkeq
	nconf=dlft.nconf
	nconf_new=dlft.nconf_new
	tau=dlft.tau
	isite=dlft.isite
etrial = te
ecps = []
if len(ap) > 0:
	ecps.append(6)
	ecps.append(6)
	ecps.append(ap)
	ecps.append(ase)

write_champ_control(champfile,mc,seed,etrial,nstep,nblk,nblkeq,nconf,nconf_new,tau,isite)
write_champ_geometry(champfile,molecular_geometry,ecps)
write_champ_wavefunction(champfile,ae,be)

#basis_type = "mixed_analytical_numerical"
if numericalORmixed == "N":
	basis_type = "numerical"
if numericalORmixed == "M":
	basis_type = "mixed_analytical_numerical"
which_analytical_basis = "Gaussian"
basis = []
basis.append(unique_atoms)
basis.append(bf_labels)
write_champ_basis(champfile,basis_type,which_analytical_basis,basis)
nn = 0
parameters = [] #
options = [] #

#num_radial_points = parameters[0]
#gridtype = parameters[1]
#gridr0 = parameters[2]
#gridarg = parameters[3]
parameters.append(1600)
parameters.append(3)
parameters.append(7.0)
parameters.append(1.00300)
#nrbas=options[0]
#igrid=options[1]
#nr=options[2]
#exp_h_bas=options[3]
#r0_bas=options[4]
#icusp=options[5]

#5 3 1600 1.003000 7.000000 0 

for basis in unique_atoms:
	options = []
	nn = nn +1
	abasis = []
	c1 = []
	c2 = []
	for a in basis:
		c1.append("GTO")
		c2.append(a)
	abasis.append(c1)
	abasis.append(c2)	
	radial_grid=compute_radial_grids(abasis, parameters)
	nrbas=len(radial_grid)-1
	options.append(nrbas)
	igrid=parameters[1]
	options.append(igrid)
	nr=len(radial_grid[0])
	options.append(nr)
	exp_h_bas=parameters[3]
	options.append(exp_h_bas)
	r0_bas=parameters[2]
	options.append(r0_bas)
	icusp=0
	options.append(icusp)
	write_champ_basis_grid("basis.%d" % nn,radial_grid,options)

if len(scf_orbitals)==0:
	print "**************************************************************"
	print "*** ERROR NO SCF ORBITALS FOUND IN GAUSSIAN LOG FILE *********"
	print "** TERIMINATING SCRIPT WITHOUT WRITING A CHAMP INPUT FILE ****"
	print "**************************************************************"
	sys.exit(-9)
write_champ_orbitals(champfile,scf_orbitals,energies,sym_labels)

write_champ_csfs(champfile,csfs)

params = champ_jastrow_parameters(unique_atoms)
new_jastrow = JASTROW(4, 4, 0.8, params)
write_champ_jastrow(champfile,new_jastrow)

write_champ_stop(champfile)

#now based on options in pseudo-vb file will print additional
optimize=True
cusps=False
correlated_sampling=False
forces=False

params = "jastrow csfs orbitals"
pvar = 0.01
ethresh = 0.001

if optimize==False and cusps==True:
	write_champ_cusps(champfile,options)
	
if correlated_sampling == True:
	write_champ_correlated_sampling(champfile,correlated_sampling_vector,geometry,istrech,alfstr)

if forces == True:
	write_champ_forces(champfile,force_types, directions)

if optimize==True:
	write_champ_optimization(champfile,pvar,ethresh,params)



#END OF PROGRAM	
