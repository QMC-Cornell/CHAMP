#!/usr/bin/env python

pseudovb_input_read_version_string ="version 1 revision 20"
##########################################################################################
### Functions related to parsing pseudo VB input files.                                ###
### Authors: Dominik Domin                                                             ###
### Last updated: May 12, 2010                            MODIFIED                     ###
##########################################################################################

import sys
import re
from string import *
from general_tools_v1 import *
from write_champ_v1 import *

#====================================================================================================================================================
def pvb_input_file_read(input_file,nbo_analysis):
	"""Parses a pseudo VB input file and info can be used to generate other input files.
	
	Original thought for creating CHAMP input files but could also be used for XMVB or other codes.
	"""
	inputfile = []
	inputfile = open(input_file,"r").readlines()
	active_start_stop=[]
	inactive_start_stop=[]
	structure_start_stop=[]
	jast_start_stop=[]
	options_start_stop=[]
	hl = []
	vb = []
	mo = []
	inactive = []
	lewis_structures = []
	lewis_weight = []
	istatus = 0
	cutoffpercentage = 0.0
	number_lewis_structures = 0
	return_start_stop_index("act", "end", active_start_stop,inputfile)
	return_start_stop_index("inact", "end", inactive_start_stop,inputfile)
	return_start_stop_index("struct", "end", structure_start_stop,inputfile)
	return_start_stop_index("jast", "end", jast_start_stop,inputfile)
	return_start_stop_index("options", "end", options_start_stop,inputfile)
	options = OPTIONS()
	#$OPTIONS
	if options_start_stop[0]!=0:
		if options_start_stop[1]==0:
			print "ERROR in pseudo-VB input file: missing end on options menu thus ignoring everything after options in the input file"
		else:
			for i in range(jast_start_stop[0]-1,jast_start_stop[1]-1):
				splitoid = split(inputfile[i])
				if len(splitoid)>0:
					if splitoid[0]=="cusp":
						options.menus.append("cusp")
						for a in splitoid[1:]:
							if a == "check":
								options.check="true"
							if a == "impose":
								options.impose="true"
							if a == "impose_opt":
								options.opt="true"
							if a == "impose_occ":
								options.occ="true"
					if splitoid[0]=="opt":
						options.menus.append("opt")
						for a in splitoid[1:]:
							if a == "jast":
								options.opt_jast="true"
							if a == "orb":
								options.opt_orb="true"
							if a == "csfs":
								options.opt_csfs="true"
							if a == "exp":
								options.opt_exp="true"
							
														
		
	#$JAST
	if jast_start_stop[0]==0 and jast_start_stop[1] == 0:
		jast = JASTROW(ijast,isc,scalek)
	else:
		if jast_start_stop[0]!=0 and jast_start_stop[1] != 0:
			ijast=0
			isc=0
			scalek=0
			parameters=[]
			for i in range(jast_start_stop[0]-1,jast_start_stop[1]-1):
				splitoid = split(inputfile[i])
				if len(splitoid)==1:
					if splitoid[0]=='NONE':
						jast = JASTROW(ijast,isc,scalek)
						break
					if splitoid[0].count("="):
						temp=splitoid[0][splitoid[0].index("="):]
						if splitoid[0][:splitoid[0].index("=")]=="ijast":
							ijast=temp
						if splitoid[0][:splitoid[0].index("=")]=="isc":
							isc=temp
						if splitoid[0][:splitoid[0].index("=")]=="scalek":
							scalek=temp
				if len(splitoid)==2:
					if splitoid[0].count("=")==1:
						if splitoid[0]=="ijast=":
							ijast=splitoid[1]
						if splitoid[0]=="isc=":
							isc=splitoid[1]
						if splitoid[0]=="scalek=":
							scalek=splitoid[1]
					if splitoid[1].count("=")==1:
						if splitoid[0]=="ijast":
							ijast=splitoid[1][1:]
						if splitoid[0]=="isc":
							isc=splitoid[1][1:]
						if splitoid[0]=="scalek":
							scalek=splitoid[1][1:]
				if len(splitoid)==3:
					if splitoid[1]=="=":
						if splitoid[0]=="ijast":
							ijast=splitoid[2]
						if splitoid[0]=="isc":
							isc=splitoid[2]
						if splitoid[0]=="scalek":
							scalek=splitoid[2]
				if len(splitoid)>3:
					parameters = inputfile[i]
			jast = JASTROW(ijast,isc,scalek,parameters)
	#$ACTIVE
	for i in range(active_start_stop[0],active_start_stop[1]):
		#print "ACT"
		splitoid = split(inputfile[i])
		#print str(splitoid)
		if len(splitoid)>0:
			if splitoid[0] == 'HL' or splitoid[0] == 'hl':
				for n in range(1,len(splitoid)):
					hl.append(int(splitoid[n]))
			if splitoid[0] == 'VB' or splitoid[0] == 'vb':
				for n in range(1,len(splitoid)):
					vb.append(int(splitoid[n]))
			if splitoid[0] == 'MO' or splitoid[0] == 'mo':
				for n in range(1,len(splitoid)):
					mo.append(int(splitoid[n]))	
	#$INACTIVE
	for i in range(inactive_start_stop[0],inactive_start_stop[1]):
		#print "IN"
		splitoid = split(inputfile[i])
		#print str(splitoid)
		if len(splitoid)==2 and splitoid[1][-1:]=='%':
			if splitoid[0] == 'none' or splitoid[0] == 'NONE' or splitoid[0] == '0':
				print "inactive orbitals = none"
			if splitoid[0] == 'list' or splitoid[0] == 'LIST' or splitoid[0] == '1':
				istatus = 1
			if splitoid[0] == 'all' or splitoid[0] == 'ALL' or splitoid[0] == '2':
				istatus = 2
			if splitoid[0] == 'all_except' or splitoid[0] == 'ALL_EXCEPT' or splitoid[0] == 'ALL_BUT' or splitoid[0] == '3':
				istatus = 3
			cutoffpercentage = float(splitoid[1][:-1])
		if istatus > 1:
			inactive = all_orbitals_listing(nbo_analysis)
			#print "all active orbitals are: " + str(inactive)
		if len(splitoid)> 1 and splitoid[1][-1:]!="%":
			#print "orbital designation"
			if istatus == 1:
				for n in range(len(splitoid)):
					inactive.append(int(splitoid[n]))
			if istatus > 1:
				#print "need to read in all orbitals and designate them"
				if istatus == 2:
					print "Possible input error: you have designated optimization of ALL orbitals and also are listing orbitals"
				if istatus == 3:
					print "Optimizing all orbitals except %s\n" % (str(splitoid))
					splitoid.sort()
					splitoid.reverse()
					#print str(splitoid)
					#inactive = all_orbitals_listing(nbo_orbitals)#need to create and test
					#print "inact " + str(inactive)
					for a in splitoid:
						#print "testing a = " + str(a)
						#print "inactive = " + str(inactive)
						for b in inactive:
							if str(a) == str(b):
								inactive.remove(b)
					#print "post removal " + str(inactive)
					istatus = -3

	#$STRUCTURE
	lewis_structure = []
	for i in range(structure_start_stop[0],structure_start_stop[1]):
		splitoid = split(inputfile[i])
		#print str(splitoid)
		if len(splitoid)==1 and splitoid[0].isdigit():
			#print "is digit!"
			number_lewis_structures = int(splitoid[0])
		if len(splitoid)>1 and len(splitoid) % 2 != 0:
			#print "is lewis structure!"
			if re.search(".",splitoid[0])!=None:
				lewis_structures.append(lewis_structure)
				lewis_structure = []
				lewis_weight.append(float(splitoid[0]))
				for a in range(1,len(splitoid)):
					lewis_structure.append(int(splitoid[a]))
			else:
				for a in range(0,len(splitoid)):
					lewis_structure.append(int(splitoid[a]))
			
	
	lewis_structures.append(lewis_structure)
	for b in lewis_structures:
		if len(b)==0:
			lewis_structures.remove(b)	
			
	#remove redundancy in the all vector if it contains orbitals that will be done us VB/MO treatment
	for a in vb:
		for b in inactive:
			if str(a)==str(b):
				inactive.remove(b)
	for a in mo:
		for b in inactive:
			if str(a)==str(b):
				inactive.remove(b)
	if istatus != 0:
		print "Optimizable Inactive orbitals are " + str(inactive)
	#append everything
	everything_from_input = []
	everything_from_input.append(cutoffpercentage)
	everything_from_input.append(vb)
	everything_from_input.append(mo)
	everything_from_input.append(inactive)
	everything_from_input.append(lewis_structures)
	everything_from_input.append(lewis_weight)
	everything_from_input.append(jast)
	everything_from_input.append(options)
	everything_from_input.append(hl)
	return everything_from_input
#end pvb_input_file_read()

#====================================================================================================================================================
def all_orbitals_listing(nbo_analysis):
	"""This function returns a list of all occupied orbitals based on NBO analysis."""
	optimizable_orbitals = []
	core = nbo_analysis[1]
	bonding = nbo_analysis[2]
	lone = nbo_analysis[3]
	if len(nbo_analysis[4])!= 0:
		print "POSSIBLE ERROR: RY ORBITALS DETECTED"
	n = []
	for a in range(len(core)):
		if a%2==0:
			if str(core[a]).isdigit():
				optimizable_orbitals.append(int(core[a]))
	for a in range(len(bonding)):
		if str(bonding[a]).isdigit():
			optimizable_orbitals.append(int(bonding[a]))
	for a in range(len(lone)):
		if str(lone[a]).isdigit():
			optimizable_orbitals.append(int(lone[a]))

	optimizable_orbitals.sort()

	return optimizable_orbitals
#end all_orbitals_listing

#====================================================================================================================================================

