#!/usr/bin/env python

xmvb_write_version_string = "Version 1 Revision 2"
##############################################################################
### A collection of general functions for creating XMVB input files.       ###
### Authors: Dominik Domin                                                 ###
### Last updated: April 26, 2010                                           ###
##############################################################################
import sys
import re
from string import *

#============================================================================================================
### writes compact format orbitals into a fil.gus format file
def write_fil_gus_format(orbitals,fil_gus_file):
	#print "writing orbitals to fil.gus"
	file_name=open(fil_gus_file,"w")
	for orbital in orbitals:
		file_name.write(" %3d" % (len(orbital)/2))
	file_name.write("\n")
	for orb in orbitals:
		i=0
		for elem in orb:
			if i % 8 == 0 and i != 0:
				file_name.write("\n")
			i = i + 1
			file_name.write(" %s " % (elem))
			if str(elem) == str("0.0000"):
				print "error %d" % i
		file_name.write("\n")

	file_name.close()
#end print_fil_gus_format()

#============================================================================================================
### writes a BFI file
def write_bfi_file(bfi_file,rf,ck,it,orb):
	print "writing %s bfi file" % (bfi_file)
	file_name = open(bfi_file,"w")
	file_name.write(" %s\n %s\n %s\n" % (rf,it,ck))
	if rf==".rwf":
		print "Warning: missing %rwf label"
	if it==".int":
		print "Warning: missing %int label"
	if ck==".chk":
		print "Warning: missing %chk label"
	file_name.write("           %d\n 0            %d\n\n" % (orb,orb))
	for i in range(orb):
		file_name.write(" %3d" % (i+1))
		if (i+1)%20==0:
			file_name.write("\n")
	file_name.close()
#end write_bfi_file()

#============================================================================================================
### writes a XMI file
def write_xmi_file(xmi_file,nelectron,nbasis,norb,nmul,struct,orb,ae):
	alpha=ae
	print "writing %s xmi file" % (xmi_file)
	file_name = open(xmi_file+".xmi","w")
	file_name.write(" %s\n" % (xmi_file))
	file_name.write(" $ctrl\n")
	file_name.write(" epg=7, iout=1, iop=2, itmax=500, nstruct=1, \n")
 	file_name.write(" nelectron= %d nbasis= %d norb= %d nmul= %d \n" % (nelectron,nbasis,norb,nmul))
 	file_name.write(" $end\n")
	file_name.write(" $struct\n")
	for i in range(alpha-(nmul-1)/2):
		file_name.write(" %d %d " % (i+1,i+1))
	for i in range(alpha-(nmul-1)/2,alpha):
		file_name.write(" %d   " % (i+1))
 	file_name.write("\n $end\n")

 	file_name.write(" $orb\n")
	for i in range(len(orb)):
		file_name.write(" %d " % len(orb[i]))
	file_name.write("\n")
	for a in orb:
		for b in a:
			file_name.write(" %d " % int(b))
		file_name.write("\n")
 	file_name.write(" $end\n")
	file_name.close()
#end write_xmi_file()
#============================================================================================================
# END OF FILE
#============================================================================================================
