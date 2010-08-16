#!/usr/bin/env python

orbital_tools_version_string ="version 1 revision 2"
##############################################################################
### A collection of general functions for transformation of orbitals.      ###
### Authors: Dominik Domin                                                 ###
### Last updated: April 27, 2010                                           ###
##############################################################################
import sys
import re
from string import *

#====================================================================================================================================================
def standard_orbital_to_compact_orbital_format(standard_orbital):
	"""Function removes zero elements and saves just the occupied coefficents and the basis function numbers.
	
	Will need to be reworked when changing storage of orbitals.	
	"""
	compact_orbital = []
	for i in range(len(standard_orbital)):
		if standard_orbital[i] != "0.0000":
			compact_orbital.append(str(standard_orbital[i]))
			compact_orbital.append(str(i+1))
	return compact_orbital
#end standard_orbital_to_compact_orbital()

#====================================================================================================================================================
def orbitals_standard_format_to_compact_format(standard_orbitals):
	"""Functions loops over orbitals to save them in compact format

	Orbitals are saved as lists.
	"""
	compact_orbitals = []
	for standard_orbital in standard_orbitals:
		compact_orbitals.append(standard_orbital_to_compact_orbital_format(standard_orbital))
	return compact_orbitals
#end orbitals_standard_format_to_compact_format()

#====================================================================================================================================================
def standard_orbital_to_number(standard_orbital):
	"""Saves basis function numbers for non-zero elements of orbital.

	Nasty set of functions to be removed ASAP or at least modified to suport new classes.
	Function used to generate RHF/ROHF wave functions.
	"""
	compact_orbital = []
	for i in range(len(standard_orbital)):
		if standard_orbital[i] != "0.00000" and standard_orbital[i] != "0.0000":
			compact_orbital.append(str(i+1))
	return compact_orbital
#end standard_orbital_to_number()

#====================================================================================================================================================
def orbitals_standard_format_to_numbers(standard_orbitals,occupation_numbers):
	"""Loops over occupied orbitals to find the basis functions used in each of them.

	Nasty set of functions to be removed ASAP or at least modified to suport new classes.
	Function used to generate RHF/ROHF wave functions.
	"""
	compact_orbitals = []
	i = -1
	for standard_orbital in standard_orbitals:
		i = i+1
		if occupation_numbers[i]=='O':
			compact_orbitals.append(standard_orbital_to_number(standard_orbital))
	return compact_orbitals
#end orbitals_standard_format_to_numbers()

#====================================================================================================================================================
def vectors_threshold_clean_up(ivector, athreshold):
	"""Returns a list of vectors that have elements smaller than a threshold value set to zero
	
	This function calls vector_threshold_clean_up()	
	"""
	rvector = []	
	for vec in ivector:
		svec = []
		svec = vector_threshold_clean_up(vec,athreshold)
		rvector.append(svec)
	return rvector
#end vectors_clean_up()

def vector_threshold_clean_up(vec,athreshold):
	"""Function sets elements of a vector that are below a threshold value to 0.0000 """
	svec = []
	svec_elem = 0.0000
	for vec_elem in vec:
		if abs(float(vec_elem)) <= athreshold:
			svec_elem = 0.0000
		else:
			svec_elem = float(vec_elem)
		svec.append("%1.4f" % svec_elem)
	return svec
#end vector_clean_up()

#====================================================================================================================================================
def vectors_atoms_clean_up(ivectors,iatoms,Catoms):
	"""Sets all elements that belong to any atoms in Catoms list to 0.0000"""
	rvectors = []
	for avector in ivectors:
		rvectors.append(vector_atoms_clean_up(avector,iatoms,Catoms))
	return rvectors
#end vectors_atom_clean_up()

def vector_atoms_clean_up(ivector,iatoms,Catoms):
	"""Sets elements that belong to Catoms to 0.0000"""
	rvector = []
	rvector = ivector
	cleanup_index = []
	for i in range(len(ivector)):
		for a in Catoms:
			if str(a)==str(iatoms[i]):
				cleanup_index.append(i)
	for i in cleanup_index:
		rvector[i]="0.0000"
		
	return rvector
#end vector_atoms_clean_up()

#====================================================================================================================================================
def vectors_symmetry_identification(ivector):
	"""Returns symmetry vectors that state which elements are non-zero.

	Vectors returned are of the type 0 0 0 1 1 0 1 0 1 1 etc
	"""	
	rvector = []
	for vec in ivector:	
		svec=[]
		svec = vector_symmetry_identification(vec)
		rvector.append(svec)
	return rvector	
#end vectors_symmetry_identification()

def vector_symmetry_identification(vec):
	"""Returns symmetry vector that states which elements are non-zero. """
	svec=[]
	svec_elem = -1		
	for vec_elem in vec:
		if abs(float(vec_elem)) > 0.0:
			svec_elem = 1
		else:
			svec_elem = 0
		svec.append(svec_elem)
	return svec
#end vector_symmetry_identification()

#====================================================================================================================================================
def vector_compare(vector_1, vector_2):
	"""Compares two vectors if elements are 0.0000 and returns a list contaning number of compatible and incompatible elements.

	A compatible element is such that vector one is 1 and vector two is 1 or 0, but if vector one is 0 then vector two must be 0.
	"""
	if len(vector_1) != len(vector_2):
		print "Error lengths of vector_1 != vector_2"
		sys.exit(0)
	number_same_elem = 0
	number_different_elements = 0
	number_incompatible_elements = 0
	number_v1_0_v2_1 = 0
	number_v1_1_v2_0 = 0
	
	for a in range(len(vector_1)):
		if vector_1[a] == 0:
			if vector_2[a] == 0:
				number_same_elem = number_same_elem + 1	
			if vector_2[a] == 1:
				number_different_elements = number_different_elements + 1
				number_v1_0_v2_1 = number_v1_0_v2_1 + 1
				number_incompatible_elements = number_incompatible_elements + 1
		if vector_1[a] == 1:
			if vector_2[a] == 1:
				number_same_elem = number_same_elem + 1	
			if vector_2[a] == 0:
				number_different_elements = number_different_elements + 1
				number_v1_1_v2_0 = number_v1_1_v2_0 + 1
		
	if len(vector_1) != (number_same_elem + number_different_elements):
		print "Error num_s + num_d != num_elem %d %d %d" % (number_same_elem,number_different_elements,len(vector_1))
		sys.exit(0)
	comp_vec = []
	comp_vec.append(number_same_elem) 
	comp_vec.append(number_different_elements)
	comp_vec.append(number_incompatible_elements)
	comp_vec.append(number_v1_0_v2_1)
	comp_vec.append(number_v1_1_v2_0)
	return comp_vec
#end vector_compare

#====================================================================================================================================================
def vector_split_up(ivector, atoms, orbital_identification):
	"""Splits up an orbital into two based on atom contributions.

	 Looks inside a single vector for atoms[0-1] needs specification of atoms in vector list (orbital_identification) 
	and returns two vectors in a list"""
	svec = []
	svec1 = []
	svec2 = []
	i = -1
	for a in orbital_identification:
		i = i + 1
		if str(a) == str(atoms[0]):
			svec1.append(ivector[i])
			svec2.append("0.0000")
		if str(a) != str(atoms[0]):
			if str(a) != str(atoms[1]):
				svec1.append(ivector[i])
				svec2.append(ivector[i])
		if str(a) == str(atoms[1]):
			svec1.append("0.0000")
			svec2.append(ivector[i])
	svec.append(svec1)
	svec.append(svec2)
	return svec
#end vector_split_up()

#====================================================================================================================================================
def remove_redundant_orbitals_listing(redundant_pnbo_virtuals):
	"""Removes redundant listing of orbitals.

	A listing is just a list of which atomic orbital elements are occupied rathe than the coefficients inside the atomic orbitals.
	"""
	#cleaning up redundant orbitals
	pnbo_virtuals_listing = []
	for n in redundant_pnbo_virtuals:
		b = 0
		for a in pnbo_virtuals_listing:
			if n == a:
				b = 1
		if b == 0:
			pnbo_virtuals_listing.append(n)
	pnbo_virtuals_listing.sort()
	return pnbo_virtuals_listing
#end remove_redundant_orbitals()

#====================================================================================================================================================
def merge_vector(Lvector,Rvector):
	"""Merges two vectors together that have been previously split."""
	return_vector = []
	if len(Lvector)!=len(Rvector):
		print "ERROR merge_vector() failed since vectors are wrong size"
		sys.exit(0)
	for a in range(len(Lvector)):
		if Lvector[a] != "0.0000" and Rvector[a] !="0.0000":
			print "ERROR the two vectors have both the same vector element as non-zero"
			print "vector L = %s\n" % (str(Lvector))
			print "vector R = %s\n" % (str(Rvector))
			
		else:
			if Lvector[a] != "0.0000":
				return_vector.append(Lvector[a])
			if Rvector[a] != "0.0000":
				return_vector.append(Rvector[a])
			if Lvector[a] == "0.0000" and Rvector[a]== "0.0000":
				return_vector.append(Rvector[a])
	return return_vector 
#end merge_vector()

#====================================================================================================================================================
def remove_same_virtuals(inactive_pnbo_virtuals_listing, mo_virtual_listing, vb_virtual_listing):
	"""Removes MO VB virtuals from PNBO virtuals."""
	ipvl = []
	for n in inactive_pnbo_virtuals_listing:
		b=0
		for m in mo_virtual_listing:
			if n==m:
				b = 1
		for m in vb_virtual_listing:
			if n==m:
				b = 1
		if b == 0:
			ipvl.append(n)
	ipvl.sort()
	return ipvl
#end remove_same_virtuals()

#====================================================================================================================================================
# END OF FILE
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
