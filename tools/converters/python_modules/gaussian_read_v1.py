#!/usr/bin/env python

gaussian_read_version_string ="version 1 revision 10"
###################################################################################################
### A collection of g98/g03/g09 related functions that are useful for writing scripts. 		###
### Authors: Dominik Domin              	                                        	###
### Last updated: October 21, 2010  	                                                	###
###################################################################################################

import sys
import re
from string import *
from general_tools_v1 import *

#====== version check =================================================================================
agt = split(general_tools_version_string)
if str(agt[1])==str(1) and str(agt[3])<str(5):
	print "FATAL ERROR: general_tools module is out of date or incompatible with gaussian_read_v1"
	sys.exit()
#end version check ====================================================================================

#=======================================================

#====================================================================================================================================================
def check_log_gaussian_version(loglines):
	'''This function checks the log file to make sure that it is made from a Gaussian(R) program and it returns the version number.
	'''
	is_gaussian = 0
	version_number = 0
	version_string = " "
	test_phrase = "the Gaussian"
	test_phrase_version = "Revision"
	for line in loglines:
		if re.search(test_phrase,line)!=None:
			is_gaussian = 1
		if is_gaussian == 1:
			if re.search(test_phrase_version,line)!=None:
				splitoid = split(line.replace(","," "))
				if len(splitoid)==4:
					version_string = str(splitoid[1])+str(splitoid[3])
					version_string = version_string.replace(".","")
					version_string = version_string.replace("A","1")
					version_string = version_string.replace("B","2")
					version_string = version_string.replace("C","3")
					version_string = version_string.replace("D","4")
					version_string = version_string.replace("E","5")
					version_string = version_string.replace("F","6")
					version_string = version_string.replace("G","7")
					version_string = version_string.replace("H","8")
					version_string = version_string.replace("I","9")
					version_string = version_string.replace("R","0")
					version_number = int(version_string)
					break
	return version_number
#end check_log_gaussian_version()

# =============================================================================================
def read_energy_number_electrons(loglines):
	'''Reads in number of alpha and beta electrons along with total SCF energy.
	
	'''
	ae = 0
	be = 0
	te = 0.0
	for line in loglines:
		if re.search("SCF Done:",line)!=None:
			t = line.replace("="," ")
			s = t.replace("SCF Done:" ," ")
			t = s.replace("A.U. after", " ")
			m = split(t)
			te = float(m[1])
		if re.search("alpha electrons",line)!=None and re.search("beta electrons",line)!=None:
			s = split(line)
			ae = int(s[0])
			be = int(s[3])
		if te != 0.0 and ae != 0:
			break
	everything = []
	everything.append(ae)
	everything.append(be)
	everything.append(te)
	return everything
#end read_energy_number_electrons()

#====================================================================================================================================================
def read_orbitals_scf(start_phrase,new_orbitals_phrase,index_start_stop,loglines):
	"""Reads the SCF orbitals from a Gaussian 09 log file and returns the orbitals.
	
	Currently returns a list of list which contains the orbitals. Does not return
	atomic orbital basis labels although easily can do so by adding bf_labels to the
	returned superstructure.	Should modify to return an orbitals class. 
	Need to see how implemented in Brian's Zori scripts.
	"""
	i = 0
	eigenvalues_array = []
	occ_array = []
	vector1=[]
	vector2=[]
	vector3=[]
	vector4=[]
	vector5=[]
	vector6=[]
	vector7=[]
	vector8=[]
	orbitals_list =[]
	numcolumns = 0
	atom_num_list = []
	label_list = []
	occ_array = []
	basis_number2=0
	bf_labels=[]
	sym_labels = []
	for i in range(index_start_stop[0],index_start_stop[1]):
		splitoid=split(loglines[i])
		s1 = split(loglines[i][:20])
		s2 = split(loglines[i][20:])
		if len(s1)==0 and len(s2)!=0:
			if str(s2[0])=='O' or str(s2[0])=='V' or str(s2[0])[-1:]=='O' or str(s2[0])[-1:]=='V' :
				numcolumns = len(s2)
				for a in range(numcolumns):
					if len(s2[a])==1:
						occ_array.append(s2[a])
					else:
						occ_array.append(s2[a][-1:])
						sym_labels.append(s2[a][1:-4])
				if len(vector1)!=0:
					orbitals_list.append(vector1)
				if len(vector2)!=0:
					orbitals_list.append(vector2)
				if len(vector3)!=0:
					orbitals_list.append(vector3)
				if len(vector4)!=0:
					orbitals_list.append(vector4)
				if len(vector5)!=0:
					orbitals_list.append(vector5)
				if len(vector6)!=0:
					orbitals_list.append(vector6)
				if len(vector7)!=0:
					orbitals_list.append(vector7)
				if len(vector8)!=0:
					orbitals_list.append(vector8)
				vector1=[]
				vector2=[]
				vector3=[]
				vector4=[]
				vector5=[]
				vector6=[]
				vector7=[]
				vector8=[]
		if len(s1)!=0:
			if s1[0]=='EIGENVALUES' or s1[0]=='Eigenvalues':
				for a in range(0,len(s2)):
					eigenvalues_array.append(s2[a])
			elif (re.search("DENSITY MATRIX",loglines[i])==None and re.search("Density Matrix",loglines[i])==None):
				if len(split(loglines[i][:4]))>0:
					basis_number = int(str(split(loglines[i][:4])[0]))
				if len(split(loglines[i][4:6]))>0:
					atom_number = int(str(split(loglines[i][4:6])[0]))	
				if len(split(loglines[i][6:11]))!=0:
					atom_label = str(split(loglines[i][6:11])[0])
				bf_label=str(split(loglines[i][11:20])[0]) 
				
				if basis_number>basis_number2:
					basis_number2 = basis_number
					atom_num_list.append(atom_number)
					label_list.append(atom_label)
					bf_labels.append(bf_label) 
				if numcolumns>0:
					if len(loglines[i])>=31:
						vector1.append(str(split(loglines[i][19:31])[0]))
						if str(vector1[-1])!=str(s2[0]):
							print "ERROR %s != %s" % (str(vector1[-1]),str(s2[0]))
							sys.exit(9)
					if len(loglines[i])>=41:
						vector2.append(str(split(loglines[i][31:41])[0]))
						if str(vector2[-1])!=str(s2[1]):
							print "ERROR %s != %s" % (str(vector2[-1]),str(s2[1]))
							sys.exit(9)
					if len(loglines[i])>=51:
						vector3.append(str(split(loglines[i][41:51])[0]))
						if str(vector3[-1])!=str(s2[2]):
							print "ERROR %s != %s" % (str(vector3[-1]),str(s2[2]))
							sys.exit(9)
					if len(loglines[i])>=61:
						vector4.append(str(split(loglines[i][51:61])[0]))
						if str(vector4[-1])!=str(s2[3]):
							print "ERROR %s != %s" % (str(vector4[-1]),str(s2[3]))
							sys.exit(9)
					if len(loglines[i])>=71:
						vector5.append(str(split(loglines[i][61:71])[0]))
						if str(vector5[-1])!=str(s2[4]):
							print "ERROR %s != %s" % (str(vector5[-1]),str(s2[4]))
							sys.exit(9)
				if numcolumns>5:
					if len(loglines[i])>=81:
						vector6.append(str(split(loglines[i][71:81])[0]))
						if str(vector6[-1])!=str(s2[5]):
							print "ERROR %s != %s" % (str(vector6[-1]),str(s2[5]))
							sys.exit(9)
					if len(loglines[i])>=91:
						vector7.append(str(split(loglines[i][81:91])[0]))
						if str(vector7[-1])!=str(s2[6]):
							print "ERROR %s != %s" % (str(vector7[-1]),str(s2[6]))
							sys.exit(9)
					if len(loglines[i])>=101:
						vector8.append(str(split(loglines[i][91:101])[0]))
						if str(vector8[-1])!=str(s2[7]):
							print "ERROR %s != %s" % (str(vector8[-1]),str(s2[7]))
							sys.exit(9)
					if numcolumns>8:
						print "Error: script can only handle at most 8 vector columns for SCF type orbitals"

	if len(vector1)!=0:
		orbitals_list.append(vector1)
	if len(vector2)!=0:
		orbitals_list.append(vector2)
	if len(vector3)!=0:
		orbitals_list.append(vector3)
	if len(vector4)!=0:
		orbitals_list.append(vector4)
	if len(vector5)!=0:
		orbitals_list.append(vector5)
	if len(vector6)!=0:
		orbitals_list.append(vector6)
	if len(vector7)!=0:
		orbitals_list.append(vector7)
	if len(vector8)!=0:
		orbitals_list.append(vector8)
	superstructure = []
	superstructure.append(orbitals_list)
	superstructure.append(occ_array)
	superstructure.append(bf_labels)
	superstructure.append(eigenvalues_array)
	superstructure.append(sym_labels)
	return superstructure
#end read_orbitals_scf()

#====================================================================================================================================================
def read_orbitals_nbo(start_phrase,new_orbitals_phrase,index_start_stop,loglines):
	"""The following function reads PNBO, PNLMO, etc orbitals for G09 log file and returns the orbitals.
	 
	Should work with PNBOs, PNLMOs, NBOs, etc and stores the atom numbers in the first vector and the atom label plus orbital label in second
	vector.	Need to rewrite this for classes.
	"""
	i = 0
	vector1=[]
	vector2=[]
	vector3=[]
	vector4=[]
	vector5=[]
	vector6=[]
	vector7=[]
	vector8=[]
	vector9=[]
	vector10=[]
	vector11=[]
	vector12=[]
	orbitals_list =[]
	numcolumns = 0
	atom_num_list = []
	label_list = []
	for i in range(index_start_stop[0],index_start_stop[1]):
		splitoid=split(loglines[i])
		linerange1 = loglines[i][:17]
		splitoid1=split(linerange1)
		linerange2 = loglines[i][17:]
		splitoid2=split(linerange2)
		if re.search(start_phrase,loglines[i])!=None or re.search("----------",loglines[i])!=None or str(loglines[i])==str(new_orbitals_phrase):
			if re.search(start_phrase,loglines[i])!=None:
				numcolumns = len(splitoid) - 1
			if str(loglines[i])==str(new_orbitals_phrase):
				if loglines[i]==new_orbitals_phrase:
					if len(orbitals_list)==0:
						if len(atom_num_list) > 0:
							orbitals_list.append(atom_num_list)
							orbitals_list.append(label_list)
				label_list=[]
				atom_num_list=[]
				numcolumns = 0
				if len(vector1)!=0:
					orbitals_list.append(vector1)
				if len(vector2)!=0:
					orbitals_list.append(vector2)
				if len(vector3)!=0:
					orbitals_list.append(vector3)
				if len(vector4)!=0:
					orbitals_list.append(vector4)
				if len(vector5)!=0:
					orbitals_list.append(vector5)
				if len(vector6)!=0:
					orbitals_list.append(vector6)
				if len(vector7)!=0:
					orbitals_list.append(vector7)
				if len(vector8)!=0:
					orbitals_list.append(vector8)
				if len(vector9)!=0:
					orbitals_list.append(vector9)
				if len(vector10)!=0:
					orbitals_list.append(vector10)
				if len(vector11)!=0:
					orbitals_list.append(vector11)
				if len(vector12)!=0:
					orbitals_list.append(vector12)
				vector1=[]
				vector2=[]
				vector3=[]
				vector4=[]
				vector5=[]
				vector6=[]
				vector7=[]
				vector8=[]
				vector9=[]
				vector10=[]
				vector11=[]
				vector12=[]		
		else:
			if numcolumns > 0:
				if len(loglines[i])>=15:
					if numcolumns >=0:
						temp_split=split(linerange1[11:16])
						temp_split2=split(linerange1[6:8])
						label_list.append(temp_split2+temp_split[1:-1])
						temp_split=split(linerange1[8:10])
						atom_num_list.append(temp_split)
				if len(loglines[i])>=23:	
					if numcolumns >=1:
						s=split(linerange2[:7])
						vector1.append(str(s[0]))
				if len(loglines[i])>=31:	
					if numcolumns >=2:
						s=split(linerange2[7:15])
						vector2.append(str(s[0]))
				if len(loglines[i])>=39:					
					if numcolumns >=3:
						s=split(linerange2[15:23])
						vector3.append(str(s[0]))
				if len(loglines[i])>=47:	
					if numcolumns >=4:
						s=split(linerange2[23:31])
						vector4.append(str(s[0]))
				if len(loglines[i])>=55:	
					if numcolumns >=5:
						s=split(linerange2[31:39])
						vector5.append(str(s[0]))
				if len(loglines[i])>=63:	
					if numcolumns >=6:
						s=split(linerange2[39:47])
						vector6.append(str(s[0]))
				if len(loglines[i])>=71:
					if numcolumns >=7:
						s=split(linerange2[47:55])
						vector7.append(str(s[0]))
				if len(loglines[i])>=79:
					if numcolumns >=8:
						s=split(linerange2[55:63])
						vector8.append(str(s[0]))
				if len(loglines[i])>=87:
					if numcolumns >=9:
						print "Warning 9+ orbitals detected in a column. Please double check output"
						vector9.append(linerange2[63:71])
					if numcolumns >=10:
						vector10.append(linerange2[71:79])
					if numcolumns >=11:
						vector11.append(linerange2[79:87])
					if numcolumns >=12:
						vector12.append(linerange2[87:95])
					if numcolumns >= 13:
						print "ERROR: script not designed to handle 13+ orbitals in a column!"
						sys.exit(0)				
	return orbitals_list
#end read_orbitals()

#====================================================================================================================================================
def read_nbo_overlap(start_phrase,new_orbitals_phrase,index_start_stop,loglines):
	"""Reads the NBO/NLMO/etc overlap matrices and returns info in list format.
	
	First vector corresponds to atom numbers
	Second vector corresponds to Atom Label plus Orbital Label
	Does not use fixed position but rather splitting and may be dangerous should entries in two columns be 
	sandwiched next to each other. Need to double check and possibly reimplement. Should also switch to more advanced classes.
	"""
	i = 0
	vector1=[]
	vector2=[]
	vector3=[]
	vector4=[]
	vector5=[]
	vector6=[]
	vector7=[]
	vector8=[]
	vector9=[]
	vector10=[]
	vector11=[]
	vector12=[]
	orbitals_list =[]
	splitoid =[]
	splitoid2=[]
	numcolumns = 0
	for i in range(index_start_stop[0],index_start_stop[1]):
		splitoid2=split(loglines[i])
		if re.search(start_phrase,loglines[i])!=None:
			numcolumns = len(splitoid2) - 1
		if re.search("----------",loglines[i])!=None:
			if numcolumns+1 != len(splitoid2):
				print "Parsing error in Overlap Matrix Read"
				break
		if str(loglines[i])==str(new_orbitals_phrase):
			if len(vector1)!=0:
				orbitals_list.append(vector1)
			if len(vector2)!=0:
				orbitals_list.append(vector2)
			if len(vector3)!=0:
				orbitals_list.append(vector3)
			if len(vector4)!=0:
				orbitals_list.append(vector4)
			if len(vector5)!=0:
				orbitals_list.append(vector5)
			if len(vector6)!=0:
				orbitals_list.append(vector6)
			if len(vector7)!=0:
				orbitals_list.append(vector7)
			if len(vector8)!=0:
				orbitals_list.append(vector8)
			if len(vector9)!=0:
				orbitals_list.append(vector9)
			if len(vector10)!=0:
				orbitals_list.append(vector10)
			if len(vector11)!=0:
				orbitals_list.append(vector11)
			if len(vector12)!=0:
				orbitals_list.append(vector12)
	
			vector1=[]
			vector2=[]
			vector3=[]
			vector4=[]
			vector5=[]
			vector6=[]
			vector7=[]
			vector8=[]
			vector9=[]
			vector10=[]
			vector11=[]
			vector12=[]
		if len(splitoid2)>0 and numcolumns > 0:
			if re.search(start_phrase,loglines[i])!=None or re.search("----------",loglines[i])!=None:
				abc = 1# basically does nothing!
			else:
				if len(splitoid2)>0:
					splitoid=split(loglines[i][17:])
				if loglines[i]==new_orbitals_phrase:
					numcolumns = 0
				if len(splitoid)-numcolumns == 0:
					if numcolumns >=1:
						vector1.append(splitoid[0])
					if numcolumns >=2:
						vector2.append(splitoid[1])
					if numcolumns >=3:
						vector3.append(splitoid[2])
					if numcolumns >=4:
						vector4.append(splitoid[3])
					if numcolumns >=5:
						vector5.append(splitoid[4])
					if numcolumns >=6:
						vector6.append(splitoid[5])
					if numcolumns >=7:
						vector7.append(splitoid[6])
					if numcolumns >=8:
						vector8.append(splitoid[7])
					if numcolumns >=9:
						vector9.append(splitoid[8])
					if numcolumns >=10:
						vector10.append(splitoid[9])
					if numcolumns >=11:
						vector11.append(splitoid[10])
					if numcolumns >=12:
						vector12.append(splitoid[11])
					if numcolumns >=13:
						print "Error: overlap matrix has too many columns!"
						sys.exit(0) #need to reimplement to throw an exception
			
	return orbitals_list
#end read_nbo_overlap()

#====================================================================================================================================================
def analyze_nbos(index_start_stop,loglines):
	"""Reads NBO analysis information and stores it in an appropriate list structure.
	
	Desperately needs to be reimplemented as Objects rather than messy lists.
	"""
	super_structure = []
	orbital_number = 0
	occupancy = 0.000
	read_orbital_number = 0
	prev_orb_num = 0
	bonding = []
	lone_pair = []
	core = []
	ryberg = []
	antibonding = []
	antiryberg = []
	antilone_pair = []
	anticore = []
	last_orb_type = "NON"
	b = []
	x = []
	for i in range(index_start_stop[0],index_start_stop[1]):
		#splitoid=split(loglines[i])
		#myline = str(loglines[i]).replace("("," ")
		myline = str(loglines[i]).replace(")"," ")
		myline2 = myline.replace("-","- ")
		myline = myline2.replace("1C","1 C")
		myline2 = myline.replace("2C","2 C")
		myline = myline2.replace("3C","3 C")
		myline2 = myline.replace("4C","4 C")
		myline = myline2.replace("1S","1 S")
		myline2 = myline.replace("2S","2 S")
		myline = myline2.replace("3S","3 S")
		myline2 = myline.replace("4S","4 S")
		myline = myline2.replace("5C","5 C")
		myline2 = myline.replace("1N","1 N")
		myline = myline2.replace("1M","1 M")
		myline2 = myline.replace("1B","1 B")
		myline = myline2.replace("1H","1 H")
		splitoid=split(myline)
		splitoid2=[]
		splitoid2.append(split(loglines[i][:7]))
		splitoid2.append(split(loglines[i][7:17]))
		splitoid2.append(split(loglines[i][17:27]))
		splitoid2.append(split(loglines[i][27:32]))
		splitoid2.append(split(loglines[i][32:36]))
		splitoid2.append(split(loglines[i][36:39]))
		splitoid2.append(split(loglines[i][39:42]))
		splitoid2.append(split(loglines[i][42:46]))
		splitoid2.append(split(loglines[i][46:]))
		if splitoid[0][-1:]==".":
			orbital_number = orbital_number + 1
			read_orbital_number = int(splitoid[0][:-1])
			if orbital_number != read_orbital_number:
				print "ERROR in script: read in orbital %d but thought it was %d!!!" % (read_orbital_number,orbital_number)
			occupancy = splitoid[1][1:-1]
					
			if splitoid[2][-1:]=="%":
				if len(b)>0:
					x.append(b)
				b = []
				b.append(int(splitoid[0][:-1]))
				b.append((splitoid[1][1:-1]))
				b.append((splitoid[2][:-1]))
				b.append(splitoid[3])
				if splitoid[3]=="BD*(" or splitoid[3]=="LP*(" or splitoid[3]=="RY*(" or splitoid[3]=="CR*(":
					b.append(int(splitoid[6]))
				else:
					b.append(int(splitoid[7]))
				if splitoid[3]=="BD":
					b.append(int(splitoid[9]))
			else:
				if (splitoid[2][0].isdigit())==isinstance('f',int):
					last_orb_type = splitoid[2][0:3]
			if last_orb_type=="BD":
				bonding.append(read_orbital_number)
				n = []
				n.append(occupancy)
				n.append(splitoid[6])
				if len(splitoid)>=10:
					n.append(splitoid[9])
				if len(splitoid)==9:
                                        n.append(splitoid[8])
				bonding.append(n)
			if last_orb_type=="LP":
				lone_pair.append(read_orbital_number)
				n = []
				n.append(occupancy)
				n.append(splitoid[6])
				lone_pair.append(n)
			if last_orb_type=="CR":
				core.append(read_orbital_number)
				n = []
				n.append(occupancy)
				n.append(splitoid[6])
				core.append(n)
			if last_orb_type=="RY":
				ryberg.append(read_orbital_number)
				n = []
				n.append(occupancy)
				n.append(splitoid[6])
				ryberg.append(n)
			if last_orb_type=="BD*":
				antibonding.append(read_orbital_number)
				n = []
				n.append(occupancy)
				n.append(splitoid[5])
				if len(splitoid)>=9:
					n.append(splitoid[8])
				if len(splitoid)==8:
					n.append(splitoid[7])
				antibonding.append(n)
			if last_orb_type=="CR*":
				antiryberg.append(read_orbital_number)
				n = []
				n.append(occupancy)
				n.append(splitoid[5])
				anticore.append(n)	
			if last_orb_type=="RY*":
				antiryberg.append(read_orbital_number)
				n = []
				n.append(occupancy)
				if splitoid[4].isdigit():
					n.append(splitoid[4])
				if splitoid[5].isdigit():
					n.append(splitoid[5])
				antiryberg.append(n)
			if last_orb_type=="LP*":
				antilone_pair.append(read_orbital_number)
				n = []
				n.append(occupancy)
				n.append(splitoid[5])
				n.append(splitoid[8])
				antilone_pair.append(n)
		if splitoid[0][:1]=="(":
			if splitoid[1][-1:]=="%":
				a = []
				a.append(splitoid[1][:-2])
				a.append(splitoid[2][:-1])
				a.append(splitoid[4])
				if last_orb_type == "BD":
					bonding.append(a)
					a = []
				if last_orb_type == "BD*":
					antibonding.append(a)
					a = []
				if last_orb_type == "NON":
					print "ERROR in script: incorrectly parsing percentage LAST_ORB_TYPE = NON !!!"
		if splitoid[0][-1:]=="%":
			b.append(splitoid[0][:-1])
			if len(splitoid[1])>=2:
				if splitoid[1][1:].isdigit():
					b.append(int(splitoid[1][1:]))
			else:
				if (splitoid[2].isdigit()):
					b.append(int(splitoid[2]))
	x.append(b)
	super_structure.append(x)
	super_structure.append(core)
	super_structure.append(bonding)
	super_structure.append(lone_pair)
	super_structure.append(ryberg)
	super_structure.append(antibonding)
	super_structure.append(antilone_pair)
	super_structure.append(anticore)
	super_structure.append(antiryberg)
	return super_structure
#end analyze_nbos()

#====================================================================================================================================================
def read_geometry(index_start_stop,loglines):
	'''This function reads in the molecular geometry.
	
	Geometry is read in and stored in Angstrom units.
	returns a list contain a list ordered as [atom_type, atomic_number, x, y, z, "A"]
	where atom_type is the index of first occurence of the element(atom) type.
	"A" indicates the coordinates are in Angstroms rather than in Bohr.
	'''
	print "Reading in the geometry"
	i_start = index_start_stop[0]
	i_end = index_start_stop[1]
	i = 0
	atomic_number = -9999
	molecule = []
	atoms = []
	for line in loglines[i_start:i_end]:
		atom = []
		i = i + 1
		a = split(line[0:9])
		b = split(line[9:20])
		c = split(line[30:46])
		d = split(line[46:58])
		e = split(line[58:70])
		splitoid = split(line)
		if i != int(splitoid[0]) or i != int(a[0]):
			print "ERROR: incorrectly read in geometry! \n %s --> %d != %d" % (line,i,int(a[0]))
		if int(splitoid[1]) != int(b[0]):
			print "ERROR in script: incorrect partition of fixed field. \n %s --> %d != %d" % (line,int(splitoid[1]),int(b[0]))
		if float(splitoid[3]) != float(c[0]):
			print "ERROR in script: incorrect partition of fixed field. \n %s --> %d != %d" % (line,float(splitoid[3]),float(c[0]))
		if float(splitoid[4]) != float(d[0]):
			print "ERROR in script: incorrect partition of fixed field. \n %s --> %d != %d" % (line,float(splitoid[4]),float(d[0]))
		if float(splitoid[5]) != float(e[0]):
			print "ERROR in script: incorrect partition of fixed field. \n %s --> %d != %d" % (line,float(splitoid[5]),float(e[0]))
		
		#checks if the atom type has been added to a list indicating the uniqueness of the atom based on those previously read in.		
		atomic_number = int(b[0])
		if atoms.count(atomic_number) > 0:
			atom.append(atoms.index(atomic_number)+1)
		else:
			atoms.append(atomic_number)
			atom.append(atoms.index(atomic_number)+1)
		
		#adds to the atom[] list the atomic number and xyz coordinates along with the designation "A" for Angstroms
		x = float(c[0])
		y = float(d[0])
		z = float(e[0])
		atom.append(atomic_number)
		atom.append(x)
		atom.append(y)
		atom.append(z)
		atom.append("A")
		molecule.append(atom)
	return molecule
#end read_geometry()

#====================================================================================================================================================
def read_basis_set(index_start_stop,loglines):
	'''If the gfinput option is selected then this function can read in the basis set.
	
	This function will read in general basis set input and store the information in appropriate list otherwise will give an error.
	Currently only works with S, P, D, F, G, H  basis functions.	
	'''
	#First check that the first line contains the phrase:
	#"AO basis set in the form of general basis input"	#ENDS WITH A BLANK LINE
	#print "RBS : "+str(loglines[index_start_stop[0]])
	if re.search("AO basis set in the form of general basis input",loglines[index_start_stop[0]])== None:
		print "ERROR: Are you sure that this log file has general basis printed out?"
		print "      IF NOT, please rerun the Gaussian job with the keyword gfinput"
	#different variables being read in
	bf_type = "L"
	num_contracts = 0
	contracts_linear_coeff = 0.0
	contracts_exponent_coeff = 0.0
	contract_coeff = 0.0
	contract_expo = 0.0
	contractions_exponents = []
	contractions_coefficients = []
	contractions_coefficientsP = []
	atoms = []
	atom = []
	i = 0
	cont_i = 0
	contraction_length = 0
	atom_charges = []
	contraction = []
	previous_bf_type = "None"
	number_of_sp_basis = 0
	for line in loglines[index_start_stop[0]+1:index_start_stop[1]]:
		splitoid = split(line)
		if len(splitoid)>0:
			if len(splitoid)==1:
				if splitoid[0]=="****":
					if len(contractions_exponents)>0 and len(contractions_coefficients)>0:
						if bf_type != "SP" and bf_type != "L":
							contraction.append(bf_type)
							contraction.append(contractions_exponents)
							contraction.append(contractions_coefficients)
							atom.append(contraction)
						else:
							if bf_type == "SP" or bf_type == "L":
								contraction.append("S")
								contraction.append(contractions_exponents)
								contraction.append(contractions_coefficients)
								atom.append(contraction)
								contraction = []
								contraction.append("P")
								contraction.append(contractions_exponents)
								contraction.append(contractions_coefficientsP)
								atom.append(contraction)
						contractions_exponents = []
						contractions_coefficients = []
						contractions_coefficientsP = []
						contraction = []
					atoms.append(atom)
					atom=[]
				else:
					print "potential problem not a **** in line= %s" % str(line)
			if len(splitoid)==2:
				if splitoid[0].count(".") == 0:
					i = i + 1
					if int(splitoid[0])!= len(atoms)+1:
						print "possible error in parsing of atom %d" % i
					else:
						atom_charges.append(int(splitoid[1]))
				else:
					previous_bf_type = bf_type 
					Xcontract_expo = splitoid[0]
					Xcontract_coeff = splitoid[1]
					if Xcontract_expo.count("D")==1:
						Xcontract_expo=Xcontract_expo.replace("D","E")
					if Xcontract_coeff.count("D")==1:
						Xcontract_coeff=Xcontract_coeff.replace("D","E")
					contract_expo = float(Xcontract_expo)
					contract_coeff = float(Xcontract_coeff)*contraction_norm_coeff
					contractions_exponents.append(contract_expo)
					contractions_coefficients.append(contract_coeff)
			if ((len(splitoid)==3 or len(splitoid)==4) and splitoid[0].count(".")==0):
				if len(contractions_exponents)>0 and len(contractions_coefficients)>0:
					if (previous_bf_type != "SP" and previous_bf_type != "L"):
						contraction.append(bf_type)
						contraction.append(contractions_exponents)
						contraction.append(contractions_coefficients)
					else:
						if (previous_bf_type == "SP" or previous_bf_type == "L"):
							contraction.append("S")
							contraction.append(contractions_exponents)
							contraction.append(contractions_coefficients)
							atom.append(contraction)
							contraction = []
							contraction.append("P")
							contraction.append(contractions_exponents)
							contraction.append(contractions_coefficientsP)
					if len(contractions_exponents)!=contraction_length:
						print "error contractions length!!!"
					contractions_exponents = []
					contractions_coefficients = []
					contractions_coefficientsP = []
					atom.append(contraction)
					contraction=[]
				bf_type = str(splitoid[0])
				contraction_length = int(splitoid[1])
				contraction_norm_coeff = float(splitoid[2])
				if len(splitoid)==4:
					contraction_norm_expo = float(splitoid[3])
				if contraction_norm_coeff != 1.0:
					print "error contraction_norm_coeff != 1.0"
				if bf_type == "L" or bf_type == "SP":
					number_of_sp_basis = number_of_sp_basis + 1
					if number_of_sp_basis == 1:
						print "Warning: L (SP) functions detected in basis set! Please double check every thing is OK."
			if (len(splitoid)==3 and splitoid[0].count(".")==1):
				previous_bf_type = bf_type
				if splitoid[0].count(".") == 0:
					print "ERROR: problem with parsing SP type basis set!"
				else:
					Xcontract_expo = splitoid[0]
					Xcontract_coeff = splitoid[1]
					Xcontract_coeffP = splitoid[2]
					if Xcontract_expo.count("D")==1:
						Xcontract_expo=Xcontract_expo.replace("D","E")
					if Xcontract_coeff.count("D")==1:
						Xcontract_coeff=Xcontract_coeff.replace("D","E")
					if Xcontract_coeffP.count("D")==1:
						Xcontract_coeffP=Xcontract_coeffP.replace("D","E")
					contract_expo = float(Xcontract_expo)
					contract_coeff = float(Xcontract_coeff)*contraction_norm_coeff
					contractions_exponents.append(contract_expo)
					contractions_coefficients.append(contract_coeff)
					contract_coeff = float(Xcontract_coeffP)*contraction_norm_coeff
					contractions_coefficientsP.append(contract_coeff)
	return atoms
#end read_basis_set()

#====================================================================================================================================================
def read_pseudopotentials(index_start_stop,loglines):
	'''The pseudopotentials (ECPs) are read in.
	
	If no pseudopotentials are present this will give an error.
	'''
	ecps = []
	ecp = []
	term = []
	w = 0
	name = ""
	for i in range(index_start_stop[0]+1,index_start_stop[1]-2):
		lh = loglines[i][:33]
		rh = loglines[i][34:90]
		so = loglines[i][90:]
		lhs = split(lh)
		if re.search("====",str(lhs))!=None:
			break
		rhs = split(rh)
		sos = split(so)
		if len(sos)>0 and len(lhs)==0:
			if float(str(sos[0]))!= 0.0:
				print "ERROR: read_pseudopotentials not designed to read SO-Coeffient!"
				w = w + 1
		if len(lhs)>0:
			if len(term)>0:
					ecp.append(term)
			if len(ecp)>0:
				name = name + "%d " % (len(ecp)-1)
				ecp.append(name)
				ecps.append(ecp)
			ecp = []
			term = []
			name = "%s-QMC GEN %d " % (periodic_table_lookup(int(lhs[1])),int(lhs[1])-int(lhs[2]))
		else:
			if len(rhs)>0:
				if (str(rhs[0]).isdigit()):
					term.append(rhs)
				else:
					if len(term)>0:
						ecp.append(term)
						term = []
	if len(term)>0:
		ecp.append(term)
	if len(ecp)>0:
		name = name + "%d " % (len(ecp)-1)
		ecp.append(name)	
		ecps.append(ecp)
	
	if w>0:
		print "************************************************************************************************************"
		print "WARNING: %d SO-Coefficients were detected but READ_PSEUDOPOTENTIALS is not designed to parse or process them" % w
		print "*********** GENERATED PSEUDOPOTENTIALS ARE INCORRECT ! *****************************************************"
	return ecps
#end read_pseudopotentials()

#====================================================================================================================================================
def read_energy_symmetry_eigenvalue(index_start_stop,loglines): #untested
	'''This function should read in Energy and orbital Symmetries and Eigenvalues.

	Needs thorough testing and possible reimplementation.
	'''
	start = index_start_stop[0]
	stop = index_start_stop[1]
	total_energy = 0.0
	OrbSym_start = 0
	OrbSym_end = 0
	OrbEigen = 0
	OrbEigen = 0
	orbital_symmetries = []
	orbital_eigenvalues = []
	occupied = []
	virtual = []
	occ_type = ""
	for line in loglines[start:stop]:
		if (re.search("SCF Done:",line))!=None:
			line_length=len(line)
			a = line.index("=")
			splitoid = split(line[a:line_length])
			b = str(splitoid[0])
			if b.isdigit()==True:
				total_energy = atof(b)
		if (re.search("Orbital symmetries:", line)) != None:
			OrbSym_start = OrbSym_start + 1
		if (re.search("The electronic state is", line)) != None:
			OrbEigen = OrbEigen + 1
			orbital_symmetries.append(occupied)
			orbital_symmetries.append(virtual)
			occupied=[]
			virtual=[]
		if (re.search("Condensed to atoms", line)) != None:
			OrbEigend = OrbEigend + 1
		if (OrbSym_start==1) and (OrbEigen==0):
			#parsing symmetries
			splitoid = split(line)
			if splitoid[0]=="Occupied":
				occ_type = tOCCUPIED
			if splitoid[0]=="Virtual":
				occ_type = tVIRTUAL
			for a in splitoid[1:]:
				if occ_type == tOCCUPIED:
					occupied.append(a)
				if occ_type == tVIRTUAL:
					virtual.append(a)
		if (OrbSym_start==1) and (OrbEigen==1) and (OrbEigend==0):
			splitoid = split(line)
			if splitoid[0]== "Alpha" or splitoid[0]=="Beta":
				if splitoid[1]=="occ.":
					for c in splitoid[3:]:
						occupied.append(c)
				if splitoid[1]=="virt.":
					for c in splitoid[3:]:
						virtual.append(c)
				if splitoid[0]=="Beta":
					print "WARNING: read_energy_symmetry_eigenvalue() not designed for BETA orbitals"
		if (OrbSym_start==1) and (OrbEigen==1) and (OrbEigend==1):
			orbital_eigenvalues.append(occupied)
			orbital_eigenvalues.append(virtual)
	n = []
	parsed_orbitals=ORBITALS(total_energy,orbital_eigenvalues,n,orbital_symmetries)
	return(parsed_orbitals)
#end  read_energy_symmetry_eigenvalue()

#====================================================================================================================================================
def parse_occupation_guess(index_start_stop,loglines):
	'''Parses the occupation string for older style log files.
	'''
	occupation_numbers = []
	start = index_start_stop[0]
	stop = index_start_stop[1]
	i = 0
	j = 0
	for i in range(start,stop):
		line = loglines[i]
		print "D>parse_occupation guess"
		line = line.replace("Occupied"," ")
		line = line.replace("("," ")
		line = line.replace(")"," ")
		splitoid = split(line)
		j = j + len(splitoid)
	for i in range(j):
		occupation_numbers.append("O")
	return occupation_numbers
#end 

#====================================================================================================================================================
#END OF FILE
#====================================================================================================================================================

