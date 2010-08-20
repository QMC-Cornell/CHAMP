#!/usr/bin/env python

write_champ_version_string = "version 1 revision 5"
############################################################
#
# Functions useful for writing to CHAMP format
# 
#  Author: Dominik Domin
#  Last Modified on August 5, 2010
############################################################

#===========================================================
# DEFAULTS for when nothing is specified
#===========================================================
class champ_default:#control
	mc =0
	seed = 1837465927472523
	etrial = -0.0
	nstep = 100
	nblk=10
	nblkeq=1
	nconf=100
	nconf_new=0
	tau=0.01
	isite=0
#end of class champ_defaults
class OPTIONS:
	"""Options class for CHAMP for printing run options in CHAMP input files used by pseudovb_input_read

	"""
	def __init__(self):
		self.menus=[]
		self.check="false"
		self.impose="false"
		self.occ="fase"
		self.opt="false"
		self.opt_jast="false"
		self.opt_orb="false"
		self.opt_csfs="false"
		self.opt_exp="false"
#end of class OPTIONS



#===========================================================
# FUNCTIONS
#===========================================================
#====================================================================================================================================================
def write_champ_pseudopotentials(filename,cecps):
	'''Writes a pseudopotential file for champ.
	'''
	print "creating pseudopotential file %s" % filename
	file_name = open(filename,"w")
	ecps = cecps[:-1]
	Lname = cecps[-1:]
	name = str(Lname[0])
	file_name.write("%s\n" % name)
	file_name.write("%d\n" % (len(ecps)))
	for ecp in ecps:
		file_name.write("%d\n" % len(ecp))
		for a in ecp:
			file_name.write("%s %s %s\n" % (str(a[2]), str(a[0]),str(a[1]))) 
	file_name.close()
#end write_champ_pseudopotentials()

#====================================================================================================================================================
def write_champ_orbitals(filename,orbitals,energies,symmetry):
	'''Prints orbitals in a format acceptable to CHAMP.

	'''
	print "writing CHAMP orbitals: Untested function please double check!"
	file_name = open(filename,"a")
	file_name.write("orbitals\n")

	file_name.write(" coefficients\n")
	for orb in orbitals:
		for elem in orb: 
			file_name.write(" %3.8f" % float(elem))
		file_name.write("\n")		
	file_name.write(" end\n")
	if len(energies)!=0:
		file_name.write(" energies\n")
		for elem in energies:
			file_name.write(" %4.9f" % float(elem))
		file_name.write("\n")
		file_name.write(" end\n")

	if len(symmetry)!=0:
		file_name.write(" symmetry\n")
		for elem in symmetry:
			file_name.write(" %s" % str(elem))
		file_name.write("\n")
		file_name.write(" end\n")

	file_name.write("end\n\n")
	file_name.close()
#end write_orbitals_champ()

#====================================================================================================================================================
def write_champ_control(filename,mc,seed,etrial,nstep,nblk,nblkeq,nconf,nconf_new,tau,isite):#need to possibly add additional settings for VMC and DMC
	'''Prints the control menu in modern champ format.'''
	print "writing champ control menu"
	file_name = open(filename,"a")
	file_name.write("control\n")
	file_name.write(" seed=%d \n etrial=%lf \n" % (seed,etrial))
	file_name.write(" nstep=%d \n nblk=%d \n nblkeq=%d \n" % (nstep,nblk,nblkeq))
	if mc=="VMC":
		print "VMC run"
		if nconf>1:
			print "warning CHAMP VMC designed to work only with nconf=1 thus forcing nconf=1 rather than %d" % (nconf)
		file_name.write(" nconf=1 nconf_new=0 tau=%lf isite=%d\n" % (tau,isite))
		#need more VMC specific features
	if mc=="DMC":
		print "DMC run"
		file_name.write(" nconf=%f nconf_new=%d tau=%lf isite=%d\n" % (nconf,nconf_new,tau,isite))
		#need more DMC specific features
	file_name.write("end\n\n")
	file_name.close()
#end write_champ_control()
 
#====================================================================================================================================================
def write_champ_geometry(filename,geometry,ecps):
	'''Prints the geometry menu in modern champ format.

	For the moment this only contains the nuclear coordinates.
	'''
	print "writing champ geometry menu"
	file_name = open(filename,"a")
	file_name.write("nuclei\n")
	ase = []
	if len(ecps)!=0:
		nloc = ecps[0]
		nquad = ecps[1]
		lpotp1 = ecps[2]
		ase = ecps[3]
		file_name.write(" nloc=%d\n" % nloc)
		file_name.write(" nquad=%d\n" % nquad)
		file_name.write(" lpotp1 ")
		for lp in lpotp1:
			file_name.write(" %d " % lp)
		file_name.write("end\n")
	file_name.write(" geometry\n")
	i = -1
	for a in geometry:
		i = i + 1
		se = 0
		if len(ase)>0:
			se = float(ase[i])
		if a[5]=="A":
			file_name.write("  %d %2.1f %3.8f  %3.8f  %3.8f\n" % (a[0],(a[1]-se),a[2]*1.889726878,a[3]*1.889726878,a[4]*1.889726878))
		else:
			file_name.write("  %d %2.1f %3.8f  %3.8f  %3.8f\n" % (a[0],(a[1]-se),a[2],a[3],a[4]))
	file_name.write(" end\n")
	file_name.write("end\n\n")
	file_name.close()
#end write_champ_geometry()

#====================================================================================================================================================
def write_champ_wavefunction(filename,alpha_electrons,beta_electrons):
	'''Prints the wavefunction menu in modern champ format.'''
	print "writing champ wavefunction menu"
	file_name = open(filename,"a")
	file_name.write("wavefunction\n")
	file_name.write("nelec=%d nup=%d \n" % (alpha_electrons+beta_electrons,alpha_electrons))
	file_name.write("end\n\n")
	file_name.close()
#end write_champ_wavefunction()

#====================================================================================================================================================
def write_champ_basis_grid(filename,grid_list,options):#untested
	'''Creates a numerical grid that should be readible by champ.

	'''
	nrbas=options[0]
	igrid=options[1]
	nr=options[2]
	exp_h_bas=options[3]
	r0_bas=options[4]
	icusp=options[5]
	file_name = open(filename,"w")
	print "writing to grid file: %s" % filename
	print "untested function please double check"
	file_name.write("%d %d %d %1.6f %1.6f %d   nrbas,igrid,nr,exp_h_bas,r0_bas,icusp\n" % (nrbas,igrid,nr,exp_h_bas,r0_bas,icusp))
	
	for i in range(len(grid_list[0])):
		for j in range(len(grid_list)):
			file_name.write("  %.12E" % (grid_list[j][i]))
		file_name.write("\n")
	file_name.close()
#end write_champ_basis_grid()

#====================================================================================================================================================
def write_champ_basis(filename,basis_type,which_analytical_basis,basis): #unfinished
	'''Prints the basis menu in modern champ format.

	For the moment it only prints analytical Slaters.
	'''
	#print "writing CHAMP basis set section: untested"
	# none|slater|gaussian|gauss-slater
	file_name = open(filename,"a")
	file_name.write("basis\n")
	#file_name.write(" numr=%d \n" % (numr) )
	#file_name.write(" ibasis=%d \n" % (ibasis))
	file_name.write(" which_analytical_basis=%s\n" % (which_analytical_basis))
	if basis_type ==  "numerical":
		file_name.write(" radial_basis_type = numerical\n")
		file_name.write(" basis_functions\n")
		#print "D> numerical basis selected"
		file_name.write("%s" % (grid_basis(basis)))
		
	if basis_type == "analytical":
		file_name.write(" radial_basis_type = analytical\n")
		file_name.write(" basis_functions\n")
		#print "D> analytical basis selected"
		file_name.write("%s" % (analytical_basis(basis)))
	if basis_type == "mixed_analytical_numerical":
		file_name.write(" radial_basis_type = numerical-analytical\n")
		file_name.write(" basis_functions\n")
		#print "D> mixed analytical numerical basis selected"
		file_name.write("%s" % (mixed_analytical_grid_basis(basis)))
		
		#write_champ_basis_grid("basis",grid_list,options)
	#must double check with latest source code
	file_name.write(" end\n")
	file_name.write("end\n\n")
	file_name.close()
#end write_champ_basis()

#====================================================================================================================================================
def mixed_analytical_grid_basis(basis):
	'''Kernel for writing mixed analytical and grid basis sets.

	   returns a string containing everything that is needed to print basis sets info

	'''
	basis_string = ""
	i = 0
	atoms = basis[0]
	labels = [["1S"],["2PX","2PY","2PZ"],["3D0","3D+1","3D-1","3D+2","3D-2"],["4F0","4F+1","4F-1","4F+2","4F-2","4F+3","4F-3"],["5G0","5G+1","5G-1","5G+2","5G-2","5G+3","5G-3","5G+4","5G-4"]]
	for iatom in atoms:
		c = 0
		i = i + 1
		basis_string = basis_string + "%d\n" % (i)
		for contraction in iatom:
			if len(contraction[1])>=1:
				if contraction[0]=="S":
					lmax = 1
				if contraction[0]=="P":
					lmax =3	
				if contraction[0]=="D":
					lmax = 5
				if contraction[0]=="F":
					lmax = 7
				if contraction[0]=="G":
					lmax = 9
				if (contraction[0]!="S" or contraction[0]!="D" or contraction[0]!="P" or contraction[0]!="F" or contraction[0]!="G")==0:
					print "ERROR: %s is not supported" % contraction[0]
				c = c + 1
				for n in range(lmax):
					if len(contraction[1])==1:
						basis_string = basis_string + "%s %5.6f %d\n" % (labels[(lmax-1)/2][n],float(str(contraction[1][0])),c)
					else:
						basis_string = basis_string + "%s 0.0 %d\n" % (labels[(lmax-1)/2][n],c)
		
	return basis_string
#end mixed_analytical_grid_basis()

#====================================================================================================================================================
def analytical_basis(basis):
	'''Kernel of write_champ_basis for analytical basis sets (Slater or Gaussian).

	'''
	print "D: using analytical_basis()"
	print "ABORTING SCRIPT UNFINISHED FUNCTION :  analytical_basis()"
	sys.exit(0)
#end

#====================================================================================================================================================
def grid_basis(basis):
	'''Kernel of write_champ_basis for numeric (grid) basis sets.

	'''
	#print "D: using mixed_analytical_grid_basis()"
	basis_string = ""
	i = 0
	atoms = basis[0]
	labels = [["1S"],["2PX","2PY","2PZ"],["3D0","3D+1","3D-1","3D+2","3D-2"],["4F0","4F+1","4F-1","4F+2","4F-2","4F+3","4F-3"],["5G0","5G+1","5G-1","5G+2","5G-2","5G+3","5G-3","5G+4","5G-4"]]
	for iatom in atoms:
		c = 0
		i = i + 1
		basis_string = basis_string + "%d\n" % (i)
		for contraction in iatom:
			if len(contraction[1])>=1:
				if contraction[0]=="S":
					lmax = 1
				if contraction[0]=="P":
					lmax =3	
				if contraction[0]=="D":
					lmax = 5
				if contraction[0]=="F":
					lmax = 7
				if contraction[0]=="G":
					lmax = 9
				if (contraction[0]!="S" or contraction[0]!="D" or contraction[0]!="P" or contraction[0]!="F" or contraction[0]!="G")==0:
					print "ERROR: %s is not supported" % contraction[0]
				c = c + 1
				for n in range(lmax):
					basis_string = basis_string + "%s %d\n" % (labels[(lmax-1)/2][n],c)
	return basis_string
#end grid_basis()

#====================================================================================================================================================
def write_champ_csfs(filename,csfs):
	'''Writes the CSF information to a CHAMP input file.

	'''
	determinants = csfs[0]
	csf_coef = csfs[1]
	dets_in_csfs = csfs[2]
	vb_weights = csfs[3]
	
	file_name = open(filename,"a")
	file_name.write("csfs\n")
	file_name.write(" determinants\n")
	file_name.write(" %s" % determinants)
	file_name.write(" end\n")
	file_name.write(" csf_coef ")
	file_name.write(" %s" %  csf_coef)
	file_name.write(" end\n")
	file_name.write(" dets_in_csfs\n")
	file_name.write(" %s" % dets_in_csfs)
	file_name.write(" end\n")
	file_name.write(" vb_weights = %s\n" % vb_weights)
	file_name.write("end\n\n")
	file_name.close()
#end write_champ_csfs()

#====================================================================================================================================================
def write_champ_jastrow(filename,jastrow):
	'''A method for printing out the jastrow section 
	
	'''
	ijas = jastrow.ijast
	isc = jastrow.isc
	scalek = jastrow.scalek
	parameters = jastrow.parameters
	file_name = open(filename,"a")
	file_name.write("jastrow\n")
	file_name.write(" ijas=%d isc=%d scalek=%2.2f \n" % (ijas,isc,scalek))
	file_name.write(" parameters \n") 
	file_name.write("%s\n" % (parameters))
	file_name.write(" end\nend\n\n")
	file_name.close()
#end

def champ_jastrow_parameters(unique_basis):
	'''Returns standard initial parameters for ijas=4 isc=4
	
	'''
	jastrow_parameters = ""
	for i in unique_basis:
		jastrow_parameters = jastrow_parameters + "0. 0. 0. 0. 0. 0. (a(iparmj),iparmj=1,nparma)\n"
	jastrow_parameters = jastrow_parameters + "0.5 1. 0. 0. 0. 0. (b(iparmj),iparmj=1,nparmb)\n"
	for i in unique_basis:
		jastrow_parameters = jastrow_parameters + "0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. (c(iparmj),iparmj=1,nparmc)\n"
	
	return jastrow_parameters
#end champ_jastrow_parameters()

#====================================================================================================================================================
def write_champ_optimization(filename,pvar,ethresh,params):
	'''Prints an optimization menu in modern Champ format.

	Need to potentially add other settings after discussing with Benoit and Julien.
	'''
	print "writing a champ optimization menu: "
	file_name = open(filename,"a")
	file_name.write("optimization\n")
	file_name.write(" parameters %s end\n" % params)
	#must double check with latest source code to see if to add more things!
	file_name.write(" p_var = %lf \n energy_threshold = %lf \n" % (pvar,ethresh))
	file_name.write("end\n\n")
	file_name.close()
#end write_champ_optimization()

#====================================================================================================================================================
def write_champ_correlated_sampling(filename,correlated_sampling_vector,geometry,istrech,alfstr):#untested
	'''Prints the correlated sampling menu in modern champ format.

	Takes in a vector of the displacements for correlated sampling.	
	'''
	print "writing champ correlated sampling menu: UNTESTED FUNCTION"
	file_name = open(filename,"a")
	file_name.write("correlated_sampling\n")
	file_name.write(" nforce = %d \n istrech = %d \n alfstr = %d \n geometry_displacements\n" % (len(correlated_sampling_vector)+1,istrech,alfstr))
 	for a in geometry:
		file_name.write(" %3.8f %3.8f %3.8f\n" % (a[2],a[3],a[4]))
	for b in correlated_sampling_vector:
		i = -1
		for a in geometry:
			i = i + 1
			file_name.write(" %3.8f %3.8f %3.8f\n" % (a[2]+b[i][0],a[3]+b[i][1],a[4]+b[i][2]))
	file_name.write(" end\n")
	file_name.write("end\n\n")
	file_name.close()
#end write_champ_correlated_sampling()

#====================================================================================================================================================
def write_champ_cusps(filename,options):
	'''Writes a custom cusp menu for CHAMP

	''' 
	print "writing champ cusps menu: UNTESTED FUNCTION"
	check=options.check
	impose=options.impose
	opt=options.opt
	occ=options.occ
	file_name = open(filename,"a")
	file_name.write("cusp\n")
	file_name.write(" check_cusp_en = %s\n impose_cusp_en = %s\n" % (check,impose))
	if opt==True or occ==True:
		file_name.write(" impose_cusp_en_opt = %s\n impose_cusp_en_occ = %s\n" % (opt,occ))
 	file_name.write("end\n\n")
	file_name.close()
#end write_champ_cusps()

#====================================================================================================================================================
def write_champ_forces(filename,force_types, directions):
	'''Writes a custom forces menu for CHAMP

	'''
	print "Writing FORCES menu: UNTESTED FUNCTION"
	print "PLEASE DOUBLE CHECK THE FUNCTION write_champ_forces() BEFORE USING"
	file_name = open(filename,"a")
        file_name.write("forces\n")
	if len(directions) > 0:
		file_name.write(" components ")
		for a in directions:
			file_name.write(" %s " % a)
		file_name.write(" end\n")
	if len(force_types)>0:
		for b in force_types:
			file_name.write(" estimator = %s \n" % b)
        file_name.write("end \n\n")
#end write_champ_forces()

#====================================================================================================================================================
def write_champ_stop(filename):
	'''Prints a stop command within a CHAMP input file

	'''
	file_name = open(filename,"a")
        file_name.write("exit\n")
	file_name.close()
#end write_champ_stop()

#====================================================================================================================================================
# END OF FILE
#====================================================================================================================================================
 
