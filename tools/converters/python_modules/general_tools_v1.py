#!/usr/bin/env python

general_tools_version_string ="version 1 revision 5"
##############################################################################
### A collection of general functions that are useful for writing scripts. ###
### Authors: Dominik Domin                                                 ###
### Last updated: August 5, 2010					   ###
##############################################################################
import sys
import re
from string import *
import math

#===================================================================================================
def atomic_number_from_element_symbol(element_symbol):
	'''Return the atomic number of an element based on its two letter symbol

	'''
	atomic_number = -999999
	
	periodic_table_element_list = ['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al', 'Si','P','S','Cl','Ar','K','Ca','Sc', 'Ti','V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Ha']
	i = 0	
	for el in periodic_table_element_list:
		i = i + 1
		if str(element_symbol)==str(el):
			atomic_number = i
			break
	if atomic_number == -999999:
		print "ERROR: %s not found in list of elements!" % element_symbol
	return atomic_number
#end atomic_number_from_element_symbol()

#===================================================================================================
def periodic_table_lookup(atomic_number):
	'''Returns element symbol based on atomic charge
	'''
	element_name =""
	periodic_table_element_list = ['X', 'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al', 'Si','P','S','Cl','Ar','K','Ca','Sc', 'Ti','V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Ha'] #X is a ghost atom
	if atomic_number > len(periodic_table_element_list):
		print "*****************************************************************************************************************************"
		print "ERROR: atomic number %d is greater than the available Periodic Table list (currently largest element has atomic number of %d)" % (atomic_number, len(periodic_table_element_list))
		print "*****************************************************************************************************************************"
	else:	
		element_name = str(periodic_table_element_list[atomic_number])
	return element_name
#end periodic_table_lookup
	
#===================================================================================================
def check_base_file_name(base_file_name):
	'''Completely useless function.

	'''
	if base_file_name[-3:]=="log":
		outputfile = open(base_file_name[:-4]+".outputfile",'w')
	else:
		outputfile = open(base_file_name+".outputfile",'w')
	return outputfile	
#end of check_base_file_name

#===================================================================================================
def punch_orbitals_to_outputfile(outputfile, orbitals):
	"""Prints all orbitals and labels to an output file.
	
	This is mostly useful for debugging, but should be easy to modify to anything else that might
	be needed.	
	"""
        i = 0
        for a in orbitals:
                i = i + 1
                if (i==1):
                        outputfile.write("ATOM LABEL NUMBERS: ")
                if (i==2):
                        outputfile.write("Atomic Orbital Labels: ")
                for b in a:
                        outputfile.write("%s " % b)
                outputfile.write("\n")
        outputfile.write("\n")
        return outputfile
#end punch_orbitals_to_outputfile

#====================================================================================================
def script_help(program_name, arguments_list, usage_description):
	"""This general function prints a help menu based on the past arguments and quits the program.
	
	This function should simply help menus in scripts but needs to transformed into a class.
	"""
	print "*************************************************************************************"
	print "*  %s " % (usage_description)
	print "*************************************************************************************"
	print "USAGE: %s %s" % (program_name, arguments_list)
	sys.exit(0)
#end of script_help

#=====================================================================================================
def return_start_stop_index(start_phrase, end_phrase, index_start_stop,loglines):
	"""The following function returns the start and end line numbers from a search thru the loglines[] list.
	
	The function returns a list containing the start index number and stop index number based on the first instances
	of the START_PHRASE and END_PHRASE.
	"""
	index_start = 0
	index_stop = 0
	i = 0
	for line in loglines:
		i = i + 1
		if re.search(start_phrase,line)!=None:
			if index_stop == 0:
				index_start = i
		if re.search(end_phrase,line)!=None:
			if index_start > 0:
				index_stop = i
				break
	index_start_stop.append(index_start)
	index_start_stop.append(index_stop)
	return index_start_stop
#end return_start_stop_index()

#=====================================================================================================
def restrained_return_start_stop_index(initial_phrase,start_phrase,end_phrase,index_start_stop,loglines):#untested
#restrained_return_start_stop_index("-- Stationary point found.","Standard orientation:","Rotational constants (GHZ):",index_start_stop,loglines)
	"""The following function returns the start and end line numbers from a search thru the loglines[] list and returns the line numbers of the last instances or the first after the initial phrase.
	"""
	index_start = 0
	index_stop = 0
	initial_index = 0
	i = 0
	for line in loglines:
		i = i + 1
		if re.search(initial_phrase,line)!=None:
			initial_index = i
		if re.search(start_phrase,line)!=None:
			index_start = i
		if re.search(end_phrase,line)!=None:
			index_stop = i
			if initial_index > 0 and index_start > 0:
				break
	index_start_stop.append(index_start)
	index_start_stop.append(index_stop)
	return index_start_stop
#end restrained_return_start_stop_index()

#=====================================================================================================
def return_unique_list(basis):
	'''Returns a unique list when inputed a redundant list
	'''
	unique_basis=[]
	unique_basis_numbering = []
	if len(basis)==0:
		return unique_basis
	unique_basis.append(basis[0])
	i = 0
	for atom in basis:
		i = 0
		for ubasis in unique_basis:
			i = i+1
			if str(atom)==str(ubasis):
				unique_basis_numbering.append(i)
				break
			else:
				if i == len(unique_basis):
					unique_basis.append(atom)
	return unique_basis
#return_unique_list

#=====================================================================================================
def is_within(item_to_check,list_of_items):
	"""Checks if item is within a given list.

	I believe an appropriate standard python function might already do this, 
	but keeping just for compatiblity of existing code.
	"""		
	for a in list_of_items:
		if str(item_to_check) == str(a):
			return 1
	return 0
#end is_within()

#=====================================================================================================
def compute_contracted_GTO_at_point(l_angular,exponents,coefficients,r_value):
	'''Computes and returns the radial value of a contracted Gaussian Type Orbital functions.
	The input parameters exponents and coefficients must be lists in the appropriate order.
	Works only for one angular moment type contraction: i.e. does not work for SP contractions.
	Normalizations taken from gamess2qmc.pl script.
	'''
	value = 0.0
	i = -1
	pi = 3.14159265358979323846264338327950288419716939937510
	if l_angular == "S":
		for alpha in exponents:
			i = i + 1
			value = value + (coefficients[i] * math.exp(-1.0*alpha*r_value**2) * (2.0/(pi**(1.0/4.0)))*((2.0*alpha)**0.75)  )
	if l_angular == "P":
		for alpha in exponents:
                        i = i + 1
			value = value + (r_value * coefficients[i] * math.exp(-1.0*alpha*r_value**2) *(math.sqrt(8.0/3.0)/(pi**(1.0/4.0)))*((2.0*alpha)**(5.0/4.0)) )
	if l_angular == "D":
                for alpha in exponents:
                        i = i + 1
                        value = value + ((r_value**2) * coefficients[i] * math.exp(-1.0*alpha*r_value**2) *(math.sqrt(16.0/15.0)/(pi**(1.0/4.0)))*((2.0*alpha)**(7.0/4.0)) )
	if l_angular == "F":
                for alpha in exponents:
                        i = i + 1
                        value = value + ((r_value**3) * coefficients[i] * math.exp(-1.0*alpha*r_value**2) * (math.sqrt(2.0/105.0)/(pi**(1.0/4.0)))*((2.0*alpha)**(9.0/4.0)) )
	if l_angular == "G":
                for alpha in exponents:
                        i = i + 1
                        value = value + ((r_value**4) * coefficients[i] * math.exp(-1.0*alpha*r_value**2) * (math.sqrt(1.0/1155.0)/(pi**(1.0/4.0)))*((2.0*alpha)**(11.0/4.0)) )
	if l_angular == "SP":
		print "ERROR: normalization of solid spherical harmonic gaussian type function is not implemented for l_angular = SP"
		print "       Your basis set contains l_angular = %s " % l_angular
		sys.exit(0)
	return value
#end compute_contracted_GTO_at_point()

#=====================================================================================================
def compute_contracted_STO_at_point(l_angular,exponents,coefficients,r_value):
        '''Currently this is an empty function!
        '''
	print "WARNING: You have called an unwritten function: compute_contracted_STO_at_point"
	if l_angular >= 0:
                print "ERROR: normalization of solid spherical harmonic slater type function is not implemented for any l_angular "
                sys.exit(0)
#end compute_contracted_GTO_at_point()

#=====================================================================================================
def compute_radial_grids(basis, parameters):
	'''Computes and returns a radial grids for CHAMP.

	'''
	radial_grid = []
	r_values = []
	line = []

	basis_types = basis[0]
	l_angulars = basis[1]
	
	num_radial_points = parameters[0]
	gridtype = parameters[1]
	gridr0 = parameters[2]
	gridarg = parameters[3]
	
	if gridtype == 1:
		for i in range(num_radial_points):
			r = gridr0 + (gridarg*i)
			r_values.append(r)
	if gridtype == 2:
		for i in range(num_radial_points):
			r = gridr0 * (gridarg**i)
			r_values.append(r)
	if gridtype == 3:
		gridr0 = gridr0/((gridarg**(num_radial_points-1))-1.0)
		for i in range(num_radial_points):
			r = 0.0	
			r = gridr0 *(gridarg**i) - gridr0
			r_values.append(r)
	radial_grid.append(r_values)
	for i in range(len(basis_types)):
		line = []
		basis_type = basis_types[i]
		l_angular = l_angulars[i][0]
		exponents = l_angulars[i][1]
		coefficients = l_angulars[i][2]
		for j in range(len(r_values)):
			r_value = r_values[j]
			if basis_type == "GTO":
				line.append(compute_contracted_GTO_at_point(l_angular,exponents,coefficients,r_value))
			if basis_type == "STO":
				line.append(compute_contracted_STO_at_point(l_angular,exponents,coefficients,r_value))
		radial_grid.append(line)
	
	return radial_grid
#compute_radial_grid()


#=====================================================================================================
# CLASSES OF COMMONLY USED OBJECTS
#=====================================================================================================
class JASTROW:
	"""Jastrow class for CHAMP"""
	#parameters = []
	#ijast = 0
	#isc = 0
	#scalek = 0
	def __init__(self, ijast, isc, scalek, parameters):
		self.ijast=ijast
		self.isc=isc
		self.scalek=scalek
		self.parameters = parameters

#end of class JASTROW

#=====================================================================================================
class ORBITALS:
	"""Orbitals class for CHAMP

	These objects contain orbitals, symmetries and energies as would be needed for Champ input files.	
	"""
	#eigenvalues = []
	#eigenvectors = []
	#symmetries = []
	def __init__(self,total_energy,energies,orbitals,symmetries):
		self.total_energy = total_energy
		self.eigenvalues=energies
		self.eigenvectors=orbitals
		self.symmetries=symmetries
#end of class ORBITALS

#=======================================================================================================
# END OF FILE
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
