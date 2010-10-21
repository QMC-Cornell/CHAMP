#!/usr/bin/python

##################################################################################
program_name = "g0x_log_2_xmvb.py"
version_num = 0.07
#
#  This script extracts orbitals from Gaussian 03/09 log files
#  and it will create BFI, FIL.GUS, and XMI files.
# 
#  This is the new generation script using Dominik's Version 1 python modules.
# 
#  Author: Dominik Domin
#  Last Modified on October 21, 2010
##################################################################################

import sys
import re
from string import *
from general_tools_v1 import *
from orbital_tools_v1 import *
from gaussian_read_v1 import *
from xmvb_write_v1 import *
from pseudovb_input_read_v1 import *

#====================================================================================================================================================
### HELP MENU strings
usage_description = "The script uses Gaussian 03/09 LOG files and creates xmvb input files (notably: fil.gus, .bfi and .xmi) "
arguments_list = " GAUSSIAN03/09_log_file optional_pseudo_vb_input"
#====================================================================================================================================================

#test that all modules are the required versions
#     axx[1] contains version number and axx[3] contains the revision number
agr = split(gaussian_read_version_string)
agt = split(general_tools_version_string)
aot = split(orbital_tools_version_string)
axw = split(xmvb_write_version_string)
apir = split(pseudovb_input_read_version_string)

a=0 
if str(agr[1])!=str(1) or int(agr[3])< int(4):
	print "FATAL ERROR: gaussian_read module is out of date or incompatible with this script"
	a = 1
if str(agt[1])!=str(1) or int(agt[3])< int(2):
	print "FATAL ERROR: general_tools module is out of date or incompatible with this script"
	a = 1
if str(aot[1])!=str(1) or int(aot[3])< int(2):
	print "FATAL ERROR: orbital_tools module is out of date or incompatible with this script"
	a = 1
if str(axw[1])!=str(1) or int(axw[3])< int(2):
	print "FATAL ERROR: orbital_tools module is out of date or incompatible with this script"
	a = 1
#if str(apir[1])!=str(1) or str(apir[3])< str(2):
#        print "FATAL ERROR: pseudovb_input_read module is out of date or incompatible with this script"
#        a = 1
if a == 1:
	sys.exit(0)
#end of testing of module versions	


#variables for getting stuff from Gaussian Log file.
ck = "" #%chk point file
rf = "" #%rwf file
it = "" #%int file
stoic = "" #stoichiometry label
fram = ""  #framework group
cha = -999999 #charge
mul = -999999 #multiplicity
bf = -9999999 #number ofbasis functions
ae = -9999999 #alpha electrons
be = -9999999 #beta electrons

#start of main program

print "** %s version %1.2lf **" % (program_name,version_num)	
if len(sys.argv)<2:
	print "ERROR: insufficiant arguments"
	print "Usage: %s %s" % (program_name,arguments_list)
	print "If you need further help type -h"
	sys.exit(0)
	
if sys.argv[1]=="help" or sys.argv[1]=="-h" or sys.argv[1]=="--h" or sys.argv[1]=="HELP" or sys.argv[1]=="-H" or sys.argv[1]=="-help":
	script_help(program_name, arguments_list, usage_description)

base_file_name = sys.argv[1]
vb_input_file_name = "NONE"
if len(sys.argv)==3:
	vb_input_file_name = str(sys.argv[2])
loglines=[]
index_start_stop = []
everything = []

if base_file_name[-3:]=="log":
	loglines = open(base_file_name,"r").readlines()
	tmp_bfn = base_file_name[:-3]
	base_file_name = tmp_bfn
else:
	loglines = open(base_file_name+".log","r").readlines()


# Extracting some labels we need
if (0==0): #some excuse not to have to change indentation
	for line in loglines:
		if re.search("%chk=",line)!=None:
			ck = line
		if re.search("%rwf=",line)!=None:
			rf = line
		if re.search("%int=",line)!=None:
			it = line
		if re.search("Stoichiometry",line)!=None:
			s = split(line)
			stoic = str(s[1])
		if re.search(" Framework group",line)!=None:
			s = split(line)
			fram = str(s[2])
		if re.search("Charge =",line)!=None and re.search("Multiplicity =",line)!=None:
			s = split(line)
			cha = int(s[2])
			mul = int(s[5])
		if re.search("basis functions",line)!=None:
			if re.search("primitive gaussians",line)!=None:
				s = split(line)
				bf = int(s[0])
		if re.search("alpha electrons",line)!=None and re.search("beta electrons",line)!=None:
			s = split(line)
			ae = int(s[0])
			be = int(s[3])

	ck=str(ck[6:-1])+".chk"
	rf=str(rf[6:-1])+".rwf"
	it=str(it[6:-1])+".int"
	if ck==".chk" or rf==".rwf" or it==" .int":
		print "************************************************************"
		print "  ***  RERUN Gaussian with %chk, %int, %rwf commands! ***"
		print "    or you must add them by hand to the generated files!!!! "
		print "************************************************************" 
	stoic=stoic[:-1]
	fram=fram[:-1]
	if (ae-be)+1 != mul:
		print "ERROR: multiplicity %d incompatible with %d alpha electrons %d beta electrons" % (mul,ae,be)
		sys.exit(9)
	if bf<=0:
		print "ERROR: number of basis function is less than or equal to zero!"
		sys.exit(8)


#nbo/nlmo orbitals read in
rsss = []
return_start_stop_index("Molecular Orbital Coefficients","DENSITY MATRIX", index_start_stop,loglines)
rsss = read_orbitals_scf("EIGENVALUES --","\n",index_start_stop,loglines)
scf_orbitals = []
occupation_numbers = []
if rsss!=None and len(rsss)!=0:
	scf_orbitals = rsss[0]
	occupation_numbers = rsss[1]
index_start_stop=[]

#nbo/nlmo analysis
return_start_stop_index("\(Occupancy\)   Bond orbital/ Coefficients/ Hybrids", "NHO Directionality and", index_start_stop,loglines)
index_start_stop[0] = index_start_stop[0]+1
index_start_stop[1] = index_start_stop[1]-4 #potentially tricky off-by-one check if an error occurs
nbo_analysis = []
nbo_analysis = analyze_nbos(index_start_stop,loglines)
index_start_stop=[]

nlmo_analysis = []
return_start_stop_index("NLMO/Occupancy/Percent from Parent NBO/ Atomic Hybrid Contributions", "PNLMOs in the AO basis:", index_start_stop,loglines)
index_start_stop[0] = index_start_stop[0]+1
index_start_stop[1] = index_start_stop[1]-3 # #potentially tricky off-by-one check if an error occurs
nlmo_analysis = []
nlmo_analysis = analyze_nbos(index_start_stop,loglines)
index_start_stop=[]



#nbo/nlmo orbitals read in
return_start_stop_index(" PNBOs in the AO basis:"," PNBO overlap matrix:", index_start_stop,loglines)
pnbo_orbitals = read_orbitals_nbo("      ---------- -------","\n",index_start_stop,loglines)
index_start_stop=[]

return_start_stop_index(" NLMOs in the AO basis:", " NLMO overlap matrix:", index_start_stop,loglines)
if index_start_stop[1]==0:
	print "Warning: NLMO overlap matrix not printed assuming that NLMO are the last vectors"
	index_start_stop=[]
	return_start_stop_index(" NLMOs in the AO basis:", " Leave Link  607 at", index_start_stop,loglines)
	if index_start_stop[1]==0:
		i=int(index_start_stop[0])
		for line in loglines[i:]:
			i = i + 1
			if len(split(line))==1:
				index_start_stop[1]=i
				break

nlmo_orbitals = read_orbitals_nbo("      ---------- -------","\n",index_start_stop,loglines)
index_start_stop=[]


return_start_stop_index(" NBOs in the AO basis:", " Second Order Perturbation Theory Analysis of Fock Matrix in NBO Basis", index_start_stop,loglines)
nbo_orbitals = read_orbitals_nbo("      ---------- -------","\n",index_start_stop,loglines)
index_start_stop=[]

return_start_stop_index(" PNLMOs in the AO basis:", " PNLMO overlap matrix:", index_start_stop,loglines)
pnlmo_orbitals = read_orbitals_nbo("      ---------- -------","\n",index_start_stop,loglines)

cfscf=orbitals_standard_format_to_compact_format(scf_orbitals)
sfpnbo=pnbo_orbitals[2:] # check if not off by one!!!
sfpnlmo=pnlmo_orbitals[2:] # check if not off by one!!!
sfnlmo=nlmo_orbitals[2:] # check if not off by one!!!
sfnbo=nbo_orbitals[2:] # check if not off by one!!!
cfpnbo=orbitals_standard_format_to_compact_format(sfpnbo)
cfpnlmo=orbitals_standard_format_to_compact_format(sfpnlmo)
cfnlmo=orbitals_standard_format_to_compact_format(sfnlmo)
cfnbo=orbitals_standard_format_to_compact_format(sfnbo)
nelectron=ae+be
nbasis=bf
norb=max(ae,be)
nmul=mul
if nmul<1:
	print "ERROR: incorrect multiplicity!"
	sys.exit(7)
if nmul>=1:
	print "script assumes RHF/ROHF trial wave function"
	if nmul!= (ae-be)+1:
		print "Warning Error in Multiplicity: 2*(Nalpha - Nbeta)+1 != Nmultiplicity"
struct =[] #useless at this point but leaving here for possibility of VB structures
orb =[]

if len(occupation_numbers) == 0:
	print "warning autogenerating orbital occupation numbers rather than reading the NBO analysis"
	for a in range(max(ae,be)):
		occupation_numbers.append('O')
	for a in range(bf):
		occupation_numbers.append('V')
	print "Orbital Occupations are assumed to be %s" % str(occupation_numbers)

if len(sfpnbo)!=0:
	write_fil_gus_format(cfpnbo,"fil.gus.pnbo")
	write_bfi_file(base_file_name+"pnbo.bfi",rf,ck,it,len(sfpnbo))
	print "warning using RHF/ROHF based orbital occupation numbers"
	orb = orbitals_standard_format_to_numbers(sfpnbo,occupation_numbers)	
	write_xmi_file(base_file_name+"pnbo.hf",nelectron,nbasis,norb,nmul,struct,orb,ae)
if len(sfnbo)!=0:
	write_fil_gus_format(cfnbo,"fil.gus.nbo")
	write_bfi_file(base_file_name+"nbo.bfi",rf,ck,it,len(sfnbo))
	orb = orbitals_standard_format_to_numbers(sfnbo,occupation_numbers)
	write_xmi_file(base_file_name+"nbo.hf",nelectron,nbasis,norb,nmul,struct,orb,ae)
if len(sfpnlmo)!=0:
	write_fil_gus_format(cfpnlmo,"fil.gus.pnlmo")
	write_bfi_file(base_file_name+"pnlmo.bfi",rf,ck,it,len(sfpnlmo))
	orb = orbitals_standard_format_to_numbers(sfpnlmo,occupation_numbers)
	write_xmi_file(base_file_name+"pnlmo.hf",nelectron,nbasis,norb,nmul,struct,orb,ae)
if len(sfnlmo)!=0:
	write_fil_gus_format(cfnlmo,"fil.gus.nlmo")
	write_bfi_file(base_file_name+"nlmo.bfi",rf,ck,it,len(sfnlmo))
	orb = orbitals_standard_format_to_numbers(sfnlmo,occupation_numbers)
	write_xmi_file(base_file_name+"nlmo.hf",nelectron,nbasis,norb,nmul,struct,orb,ae)
if len(cfscf)!=0:
	write_fil_gus_format(cfscf,"fil.gus.scf")
	write_bfi_file(base_file_name+"scf.bfi",rf,ck,it,len(cfscf))
	orb = orbitals_standard_format_to_numbers(scf_orbitals,occupation_numbers)
	write_xmi_file(base_file_name+"scf.hf",nelectron,nbasis,norb,nmul,struct,orb,ae)

if vb_input_file_name != "NONE":
	inactive_orbitals =[]
	print "****************\nParsing pseudo-VB input file\n****************\n"
	vb_input=pvb_input_file_read(vb_input_file_name,nbo_analysis)
        cutoffpercentage = vb_input[0]
	spec_vb =  vb_input[1]
	spec_mo =  vb_input[2]
	spec_inactive = vb_input[3]
	spec_hl = vb_input[8]
	#using 	nlmo_analysis
	print "cleaning up pnlmo based on cutoff of %s" % cutoffpercentage
	spec_inactive.sort()
	inactive_orbitals_listing = spec_inactive
	inactive_orbital_occupancy = []
	orbital_identification = pnlmo_orbitals[0]
	for i in inactive_orbitals_listing: #finding our PNLMO inactive orbitals
        	clean_up_atoms = []
        	if nlmo_analysis[0][i-1][0]!=i:
                	print "ERROR: orbital %d not found NLMO ANALYSIS" % (i)
                	sys.exit(0)
        	inactive_orbital_occupancy.append(float(nlmo_analysis[0][i-1][1]))
        	for j in range(4,len(nlmo_analysis[0][i-1])):
                	if j%2==1:
                        	if float(nlmo_analysis[0][i-1][j]) <= float(cutoffpercentage):
                                	clean_up_atoms.append(nlmo_analysis[0][i-1][j-1]) #cleaning them up
        	ianctive_orbital = []
        	inactive_orbital = pnlmo_orbitals[i+1]
        	inactive_orbitals.append(vector_atoms_clean_up(inactive_orbital,orbital_identification,clean_up_atoms))
	cpnmlo=inactive_orbitals
	inactive_orbitals=[]
	cfcpnmlo=orbitals_standard_format_to_compact_format(cpnmlo)
	if len(cfcpnmlo)!=0:
        	write_fil_gus_format(cfcpnmlo,"fil.gus.Cpnlmo")
        	write_bfi_file(base_file_name+"Cpnlmo.bfi",rf,ck,it,len(cfcpnmlo))
        	orb = orbitals_standard_format_to_numbers(cfcpnmlo,occupation_numbers)
        	write_xmi_file(base_file_name+"Cpnlmo.hf",nelectron,nbasis,norb,nmul,struct,orb,ae)
	#cleaining up nlmo
	print "cleaning up nlmo based on cutoff of %s " % cutoffpercentage
        spec_inactive.sort()
        inactive_orbitals_listing = spec_inactive
        inactive_orbital_occupancy = []
        orbital_identification = nlmo_orbitals[0]
        for i in inactive_orbitals_listing: #finding our PNLMO inactive orbitals
                clean_up_atoms = []
                if nlmo_analysis[0][i-1][0]!=i:
                        print "ERROR: orbital %d not found NLMO ANALYSIS" % (i)
                        sys.exit(0)
                inactive_orbital_occupancy.append(float(nlmo_analysis[0][i-1][1]))
                for j in range(4,len(nlmo_analysis[0][i-1])):
                        if j%2==1:
                                if float(nlmo_analysis[0][i-1][j]) <= float(cutoffpercentage):
                                        clean_up_atoms.append(nlmo_analysis[0][i-1][j-1]) #cleaning them up
                ianctive_orbital = []
                inactive_orbital = nlmo_orbitals[i+1]
                inactive_orbitals.append(vector_atoms_clean_up(inactive_orbital,orbital_identification,clean_up_atoms))
        cnmlo=inactive_orbitals
        inactive_orbitals=[]
        sfcnlmo=orbitals_standard_format_to_compact_format(cnmlo)
        if len(sfcnlmo)!=0:
                write_fil_gus_format(sfcnlmo,"fil.gus.Cnlmo")
                write_bfi_file(base_file_name+"Cnlmo.bfi",rf,ck,it,len(sfcnlmo))
                orb = orbitals_standard_format_to_numbers(sfcnlmo,occupation_numbers)
                write_xmi_file(base_file_name+"Cnlmo.hf",nelectron,nbasis,norb,nmul,struct,orb,ae)
	#using 	nbo_analysis


sys.exit()


