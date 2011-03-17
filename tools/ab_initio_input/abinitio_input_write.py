#!/usr/bin/python
##########################################################################
##### Writes various ab initio chemistry input files                 #####
##### Dominik Domin  LAST MODIFIED: Feb 1, 2011                      #####
version_number = 1
subversion_number = 1
##########################################################################

from string import *

### Generating GAUSSIAN input file (g98/g03/g05/g09)
def write_gaussian(name1,name2,nmul,pseudo_yes,charge,ename,unique_atom_list_ename,unique_pseudo_list_ename,pseudo,xcor,ycor,zcor):
	'''Creates a standard hartree-fock (RHF or ROHF) single point calculation.
	Also does NBO analysis...
	'''
	outputoccfile = open(name1+".com",'w')
	outputoccfile.write("%mem=2gb\n")
	outputoccfile.write("%%chk=%s\n" %(name2[0:8]))
	outputoccfile.write("%%rwf=%s\n" %(name2[0:8]))
	outputoccfile.write("%%int=%s\n" %(name2[0:8]))
	if nmul == 1:
		if pseudo_yes == "T":
			outputoccfile.write("# rhf gen pseudo=read pop=(full,nboread) scf=(nodirect,tight) noraf nosym \n")
		else:
			outputoccfile.write("# rhf gen pop=(full,nboread) scf=(nodirect,tight) noraf nosym \n")
	if nmul > 1:
		if pseudo_yes == "T":
			outputoccfile.write("# rohf gen pseudo=read pop=(full,nboread) scf=(nodirect,tight) noraf nosym \n")
        	else:
			outputoccfile.write("# rohf gen pop=(full,nboread) scf=(nodirect,tight) noraf nosym \n")
	outputoccfile.write("# gfinput 5d 7f\n")
	outputoccfile.write("\n %s \n\n%d %d\n" % (name1,charge,nmul))
	for i in range(len(ename)):
		outputoccfile.write("%s %lf %lf %lf\n" %(ename[i],xcor[i],ycor[i],zcor[i]))
	for j in unique_atom_list_ename:
		outputoccfile.write("\n")
		gbas = open(str(j),'r').readlines()	
		for line in gbas:
			if len(split(line))>0:
				outputoccfile.write("%s" % line)			
		outputoccfile.write("****")
	outputoccfile.write("\n")
	if pseudo_yes == "T":
		outputoccfile.write("\n")
		for k in unique_pseudo_list_ename:
			pseudo = open(k,'r').readlines()
			for line in pseudo:
				if len(split(line))>0:
					outputoccfile.write("%s" % line)
		outputoccfile.write("\n")	
	else:
		outputoccfile.write("\n")
	outputoccfile.write("$NBO AONBO AOPNBO AONLMO AOPNLMO SPNBO SPNLMO $END\n\n")
	outputoccfile.close()
	print "All done writing a Gaussian formatted input file for %s" % str(name1)
#end write_gaussian()

### Generating GAMESS(US) input file 
def write_gamess(name1,name2,nmul,pseudo_yes,charge,ename,unique_atom_list_ename,unique_pseudo_list_ename,pseudo,xcor,ycor,zcor,zcharges,ecp_path,basis_path):
	'''Creates a standard Hartree-Fock (RHF or ROHF single point calculation)
	'''
	## 	Generating gamess input file
	outputoccfile = open(name1+".inp",'w')
	if nmul == 1:        
		outputoccfile.write(" $CONTRL ISPHER=1 SCFTYP=RHF PLTORB=.TRUE. UNITS=ANGS\n")
		if pseudo_yes == "T":
			outputoccfile.write("  ICHARG=%d ECP=READ LOCAL=BOYS MAXIT=50 $END\n" % (charge))
		else:
			outputoccfile.write("  ICHARG=%d LOCAL=BOYS MAXIT=50  $END\n" % (charge))
		outputoccfile.write(" $SYSTEM TIMLIM=10000 MEMORY=90000000 $END\n $SCF SHIFT=.TRUE. DIRSCF=.TRUE. $END\n $GUESS GUESS=HUCKEL $END\n $DATA\n")

	else:
		outputoccfile.write(" $CONTRL ISPHER=1 SCFTYP=ROHF MULT=%d PLTORB=.TRUE. UNITS=BOHR\n" % (nmul))
		if pseudo_yes == "T":
			outputoccfile.write("  ICHARG=%d ECP=READ LOCAL=BOYS MAXIT=75 $END\n" % (charge))
		else:
			outputoccfile.write("  ICHARG=%d LOCAL=BOYS MAXIT=75  $END\n" % (charge))
		outputoccfile.write(" $SYSTEM TIMLIM=100000 MEMORY=90000000 $END\n")
		outputoccfile.write(" $SCF SHIFT=.TRUE. DIRSCF=.TRUE. $END\n")
		outputoccfile.write(" $GUESS GUESS=HUCKEL $END\n")
		outputoccfile.write(" $DATA\n")
	outputoccfile.write("%s\nC1\n" % name2)
	for i in range(len(ename)):
		outputoccfile.write("%s %1.1f %lf %lf %lf\n" %(ename[i],zcharges[i],xcor[i],ycor[i],zcor[i]))
		gbas = []
		temp_name = str(basis_path)+ "/" + str(ename[i])
		temp2 = temp_name.replace("//","/")
		temp_name = temp2.replace("//","/")
		if unique_atom_list_ename.count(str(temp_name))>0 :
			gbas = open(str(basis_path)+str(ename[i]),'r').readlines()
		else:
			gbas = open(ename[i],'r').readlines()
		for line in gbas:
			if len(split(line))>0:
				outputoccfile.write("%s" % line)
		outputoccfile.write("\n")
	outputoccfile.write("$END\n\n")
	
	if pseudo_yes == "T":
		outputoccfile.write(" $ECP\n")
		already_used = []
		for a in ename:
			if already_used.count(str(a))>0:
				outputoccfile.write("%s-ECP\n" % str(a))
			else:	
				print str(unique_pseudo_list_ename)
				temp_name = str(ecp_path)+"/"+str(a)+".ecp"
				temp2 = temp_name.replace("//","/")
				temp_name = temp2.replace("//","/")
				if unique_pseudo_list_ename.count(str(a)+".ecp")>0 or unique_pseudo_list_ename.count(str(temp_name))>0:
					already_used.append(str(a))
					if unique_pseudo_list_ename.count(str(a)+".ecp")>0:
						ecp_index = unique_pseudo_list_ename.index(str(a)+".ecp")
					if unique_pseudo_list_ename.count(str(temp_name))>0:
						ecp_index = unique_pseudo_list_ename.index(str(temp_name))
					pseudo=open(str(unique_pseudo_list_ename[ecp_index]),'r').readlines()
					for line in pseudo:
						outputoccfile.write("%s" % line)
				else:
					already_used.append(str(a))
					outputoccfile.write("%s-ECP NONE\n" % str(a))
					
		outputoccfile.write(" $END\n")
	outputoccfile.write("\n")
	outputoccfile.close()
	print "All done writing a GAMESS (US) formatted input file for %s" % str(name1)
#end write_gamess()

### Generate MOLPRO input file
def write_molpro(mem,name1,name2,nmul,pseudo_yes,charge,ename,unique_atom_list_ename,unique_pseudo_list_ename,pseudo,xcor,ycor,zcor,zcharges,ecp_path,basis_path):
	''' Creates a simple molpro Hartree-Fock calculation... or that is what it thinks its doing...
	
	'''
	## 	Generating gamess input file
	outputoccfile = open(name1+"-mp.in",'w')
	outputoccfile.write("***,%s\n" % name1)
	outputoccfile.write("memory,%d,m\ngeomtyp=xyz\n" % mem)
	outputoccfile.write("set,CHARGE=%d\n" % charge)
	outputoccfile.write("geometry={nosym; NOORIENT;\n  %d\n" % len(ename))
	outputoccfile.write("%s\n" % name2)
	for i in range(len(ename)):
		outputoccfile.write("%s   %lf   %lf   %lf\n" % (ename[i],xcor[i],ycor[i],zcor[i]) )
	outputoccfile.write("}\n")
	outputoccfile.write("basis={\n")
	if pseudo_yes == "T":
		for k in unique_pseudo_list_ename:
			pseudo = open(k,'r').readlines()
			for line in pseudo:
				if len(split(line))>0:
					outputoccfile.write("%s" % line)
		outputoccfile.write("\n")
	for j in unique_atom_list_ename:
		outputoccfile.write("\n")
		gbas = open(str(j),'r').readlines()	
		for line in gbas:
			if len(split(line))>0:
				outputoccfile.write("%s" % line)
	outputoccfile.write("}\n")
	outputoccfile.write("HF\n")		
	outputoccfile.write("--\n")
	print "All done writing a Molpro formatted input file for %s" % str(name1)
	print "**warning untested** input files might not run"
#end write_molpro()

### Generate NWCHEM input file
def write_nwchem(mem,name1,name2,nmul,pseudo_yes,charge,ename,unique_atom_list_ename,unique_pseudo_list_ename,pseudo,xcor,ycor,zcor,zcharges,ecp_path,basis_path):
	'''Another untested input file writing mechanism
	'''
	outf = open(name1+"-nw.in",'w')
	outf.write("start %s\nmemory total %d mb\necho\n"% (name1,mem))
	outf.write("geometry noautosym noautoz units angstrom\n")
	for i in range(len(ename)):
		outf.write("%s %s %s %s\n" % (ename[i],xcor[i],ycor[i],zcor[i]))
	outf.write("end\n\n")
	outf.write("basis\n")
	for iatom in unique_atom_list_ename:
		lines = open(iatom,'r').readlines()
		for line in lines:
			outf.write("%s" % line)
	outf.write("end\n\n")
	if pseudo_yes =="T":
		outf.write("ecp\n")
		for iatom in unique_pseudo_list_ename:
			lines = open(iatom,'r').readlines()
			for line in lines:
				outf.write("%s" % line)
		outf.write("end\n\n")
	
	outf.write("task hf \n\n")
	
	print "All done writing a NWChem formatted input file for %s" % str(name1)
	print "**warning untested** input files might not run"
#end write_nwchem()
