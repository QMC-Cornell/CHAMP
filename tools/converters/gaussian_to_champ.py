#!/usr/bin/python
#!/usr/local/bin/python2.4

###################################################################################
#
# GAUSSIAN to CHAMP convertion script
# Author: J. Toulouse
# Created: 23 Jan 2008
# Last modified: 28 Apr 2009
#
###################################################################################

# modules
import sys
import string
import re
import operator
import getopt
from Numeric import *

# ========================================================
# Global variables and initializations
# ========================================================

l_sort_basis = False

ncent = 0
nuclear_coordinates = []
nuclear_charges = []
cent = []
cent_to_cent_type = []
cent_type_to_first_cent = []
nctype = 0
nuclear_charges_unique = []
iwctype = []
angstrom_to_bohr = 1.889716165

total_charge = 0
nelec = 0
nup = 0
ndn = 0

nbasis = 0
basis_atom_indexes_allsubsection = []
basis_atom_indexes = []
basis_labels_allsubsection = []
basis_labels = []
basis_labels_without_number = []
basis_labels_by_cent_type = []
basis_labels_without_number_by_cent_type = []

orbital_symmetry = True
norb = 0
orbital_indexes = []
orbital_occupations = []
orbital_symmetries = []
orbital_eigenvalues = []
orbital_coefficients = []
det_up_orb_lab = []
det_dn_orb_lab = []

ndet = 1
ncsf = 1
csf_coef = [1.]
ndet_in_csf = [1]
iwdet_in_csf = [[1]]
cdet_csf = [[1.]]

energy = 0.

nparma = 4
nparmb = 5
nparmc = 15
nparmcsf = 0

# ========================================================
def help_menu():
# ========================================================
# print help menu
# --------------------------------------------------------
  print "USAGE: %s" % sys.argv[0], "-i gaussian_output_file -o champ_input_file"
  print "options: -h: print this menu"
  print "         -s: sort basis functions in standard CHAMP order"
  
  sys.exit(0)

# ========================================================
def not_empty_string (s):
  "filter not empty strings"
  return s != ''
  sys.exit(0)

# ========================================================
def read_total_charge ():
  "read the total charge of the molecule."
  global total_charge
  found = False

  for line in lines:
    object_matched = re.search ("\s+Charge\s*\=\s*(\S+)", line)
    if object_matched:
      object = object_matched.group()
      total_charge = int(float((string.split (object, '='))[1]))
      found = True
      break
      
  if not found:
    print "\nERROR: total charge not found"
    sys.exit(0)

#  print "total charge = ", total_charge

# ========================================================
def read_electron_number ():
  "read the number of electrons."
  global nelec, nup, ndn
  found = False

  for line in lines:
    object_matched = re.search ("\d+\s+alpha electrons\s+\d+\s+beta electrons", line)
    if object_matched:
      object = object_matched.group()
      object_splitted = string.split (object)
      nup = int(object_splitted[0])
      ndn = int(object_splitted[3])
      found = True
      break

  if not found:
    print "\nERROR: alpha and beta electron numbers not found"
    sys.exit(0)

#  print "number of electrons read."
  nelec = nup + ndn
#  print "nelec = ", nelec
#  print "nup = ", nup
#  print "ndn = ", ndn

# ========================================================
def read_geometry ():
  "read geometry parameters of the molecule."
  global ncent, nctype, cent
  global nuclear_charges_unique
  global cent_type_to_first_cent

  found = False
  j = 0
  for i in range(len(lines)):
    if re.search ("Input orientation:", lines[i]):
      if re.search ("Coordinates", lines[i+2]):
        if re.search ("Number", lines[i+3]):
          if re.search ("-------", lines[i+4]):
            j = i + 5
            geometry_list = []
            while (j < len(lines)):
              if re.search ("-------", lines[j]):
                found = True
                break
              geometry_list.append (string.split (lines[j]))
              j = j + 1
            if found:
              break

#  if Standard orientation (ie calculations with symmetry) is found, take it
  j = 0
  for i in range(len(lines)):
    if re.search ("Standard orientation:", lines[i]):
      found = False
      if re.search ("Coordinates", lines[i+2]):
        if re.search ("Number", lines[i+3]):
          if re.search ("-------", lines[i+4]):
            j = i + 5
            geometry_list = []
            while (j < len(lines)):
              if re.search ("-------", lines[j]):
                found = True
                break
              geometry_list.append (string.split (lines[j]))
              j = j + 1
            if found:
              break

  if not found:
    print "\nERROR: geometry not found"
    sys.exit(0)

#  print "geometry read"
#  print "geometry_list=", geometry_list

  ncent = len(geometry_list)
  for i in range(ncent):
    nuclear_charges.append (int(geometry_list[i][1]))
    nuclear_coordinates.append (map(float, geometry_list[i][3:6]))

  nuclear_charges_unique = list(set(nuclear_charges))
  nctype = len(nuclear_charges_unique)

  for i in range(ncent):
    cent_to_cent_type.append (nuclear_charges_unique.index (nuclear_charges[i])+1)
#  print "cent_to_cent_type=",cent_to_cent_type

  for i in range(nctype):
    cent_type_to_first_cent.append (nuclear_charges.index (nuclear_charges_unique[i])+1)
#  print "cent_type_to_first_cent=",cent_type_to_first_cent

  cent = (angstrom_to_bohr * array(nuclear_coordinates)).tolist()

  for i in range(ncent):
    iwctype.append(nuclear_charges_unique.index (nuclear_charges[i])+1)

#  print "nuclear_charges =", nuclear_charges
#  print "nuclear_charges_unique =", nuclear_charges_unique
#  print "iwctype =", iwctype
#  print "cent =",cent

# ========================================================
def read_energy ():
  """read the calculated energy."""
  global energy
  found = False

  for line in lines:
    if re.search ("SCF Done:", line):
      object_splitted = string.split (line)
      energy = float(object_splitted [4])
      found = True
      break

  if not found:
    print "\nERROR: energy not found"
    sys.exit(0)

#  print "SCF energy read"

# ========================================================
def read_basis ():
  "read one-electron basis."
  global nbasis
  found = False

  for line in lines:
    if re.search ("^\s+\d+\s+basis functions", line):
      object_splitted = string.split (line)
      nbasis = int(object_splitted [0])
      found = True
      break

  if not found:
    print "\nERROR: number of basis functions not found"
    sys.exit(0)

#  print "nbasis = ", nbasis

# ========================================================
def read_orbitals ():
  "read molecular orbitals."
  global norb
  global orbital_coefficients
  global basis_atom_indexes, basis_labels
  global basis_labels_without_number
  global basis_labels_by_cent_type
  global basis_labels_without_number_by_cent_type

  basis_atom_indexes_subsection = []
  basis_labels_subsection = []
  orbital_coefficients_subsection = []
  found = False
  j = 0
  subsection = 0
  begin_subsection = False
  for i in range(len(lines)):

#   begin of orbital section
    if re.search ("Molecular Orbital Coefficients", lines[i]):
      found = True
      j = i + 1
      while j <len(lines): 

#       read orbital indexes, occupations and eigenvalues
#       example:                       1         2         3         4         5
#                                      O         O         O         O         O
#                EIGENVALUES --   -26.43129 -26.43113  -1.76868  -1.50082  -0.81349
        if re.search ("^(\s+\d+)+$", lines[j]):
          
          orbital_indexes.extend (string.split(lines[j]))
#          print "line ",j,", orbital indexes:",lines[j]
          j = j + 1
          if not re.search ("([OV]\s+)+", lines[j]):
#            print "\nline: ", lines[j]
            print "\nERROR: the orbital section is not in the expected format"
            sys.exit(0)
          orbital_occupations.extend (string.split(lines[j]))
#          print "line ",j,", orbital occupation:",lines[j]
          j = j + 1

#         read eigenvalues if correct format is found
#          if not re.search ("EIGENVALUES \-\-\s+(\*)+", lines[j]):
          if not re.search ("EIGENVALUES \-\-\s+([0-9\-\.]+\s+)+", lines[j]):
#              print "\nline: ", lines[j]
              print "\nWARNING: the orbital eigenvalues are not in the expected format (but they are not used anyway)"
#              sys.exit(0)
          else:
            orbital_eigenvalues.extend ((string.split(lines[j]))[2:])
#            print "line ",j,", orbital eigenvalues:",lines[j]

          subsection = subsection + 1
          if (subsection >= 2):
            basis_atom_indexes_allsubsection.append (basis_atom_indexes_subsection)
            basis_labels_allsubsection.append (basis_labels_subsection)
            basis_atom_indexes_subsection = []
            basis_labels_subsection = []
            orbital_coefficients.extend ((transpose(orbital_coefficients_subsection)).tolist())
            orbital_coefficients_subsection = []

#       read the line of orbital coefficients on the first basis function on each atomic center
#       example: 1 1   F  1S          0.20743   0.20641   0.01669   0.01512   0.00000
        elif re.search ("^\s+\d+\s+\d+\s+[A-Za-z]+\s+\d+[SPDFGH][A-Z\+\-\d\ ]+(\s+[\d\-]+)",lines[j]):
#          print "line ",j,", orbital coefficient new center:",lines[j]
          line_splitted = string.split (lines[j])
          atom_index = line_splitted[1]
          basis_atom_indexes_subsection.append (atom_index)
          basis_labels_subsection.append (line_splitted[3])
          orbital_coefficients_subsection.append (map(float,line_splitted[4:]))

#       read normal line of orbital coefficients
#       example: 2        2S          0.50758   0.50920  -0.21505  -0.22492   0.00000
        elif re.search ("^\s+\d+\s+\d+[SPDFGH][A-Z\+\-\d\ ]+(\s+[\d\-]+)",lines[j]):
#          print "line ",j,", orbital coefficient:",lines[j]
          lines[j] = string.replace(lines[j],'D 0','D0')
          line_splitted = filter (not_empty_string,string.split(lines[j],' '))
#          print "line_splitted=",line_splitted
          basis_atom_indexes_subsection.append (atom_index)
          basis_labels_subsection.append (line_splitted[1])
#          print "line_splitted[2:]=",line_splitted[2:]
          orbital_coefficients_subsection.append (map(float,map(string.strip,line_splitted[2:])))

#       DENSITY MATRIX marks the end of the orbital section
        elif re.search ("DENSITY MATRIX.", lines[j]):
          basis_atom_indexes_allsubsection.append (basis_atom_indexes_subsection)
          basis_labels_allsubsection.append (basis_labels_subsection)
          orbital_coefficients.extend ((transpose(orbital_coefficients_subsection)).tolist())
          break

        j = j + 1
      break
       
  if not found:
    print "\nERROR: orbital coefficients not found"
    sys.exit(0)

#  print "orbitals read"

  norb = len(orbital_indexes)
#  print "basis_labels_allsubsection=",basis_labels_allsubsection
#  print "basis_atom_indexes_allsubsection=",basis_atom_indexes_allsubsection
#  print "orbital_coefficients=",orbital_coefficients
#  print "orbital_indexes=", orbital_indexes
#  print "orbital_occupations=",orbital_occupations
#  print "orbital_eigenvalues=",orbital_eigenvalues

# get orbital symmetries
# example: (SGG)--O  (SGU)--O  (SGG)--O  (SGU)--O  (PIU)--O
  for i in range(norb):
   object_matched = re.search ("\(\S+\)", orbital_occupations[i])
   if object_matched:
      object = object_matched.group()
      object = string.replace(object,'(','')
      object = string.replace(object,')','')
#      print object
      orbital_symmetries.append (object)
   else:
      orbital_symmetries.append ('?')

#  print "orbital_symmetries=",orbital_symmetries

# keep only the basis atom indexes and basis labels of the first subsection
  basis_atom_indexes = basis_atom_indexes_allsubsection[0]
  basis_labels = basis_labels_allsubsection[0]
#  print "basis_atom_indexes=",basis_atom_indexes
#  print "basis_labels=",basis_labels

# check dimensions
  if (len(basis_atom_indexes) != nbasis):
    print "\nERROR: dimension of basis_atom_indexes =",len(basis_atom_indexes)," is different from nbasis =",nbasis
    sys.exit(0)
  if (len(basis_labels) != nbasis):
    print "\nERROR: dimension of basis_labels =",len(basis_labels)," is different from nbasis =",nbasis
    sys.exit(0)

# check if basis atom indexes of other subsections are the same
  for i in range(1,len(basis_atom_indexes_allsubsection)):
    if (basis_atom_indexes_allsubsection[i] != basis_atom_indexes):
      print "\nERROR: sequences of basis atom indexes are not identical for each orbital subsection!"
      print "basis atom indexes for subsection # 1 :",basis_atom_indexes
      print "basis atom indexes for subsection #",i+1,":",basis_atom_indexes_allsubsection[i]
      sys.exit(0)

# check if basis labels of other subsections are the same
  for i in range(1,len(basis_labels_allsubsection)):
    if (basis_labels_allsubsection[i] != basis_labels):
      print "\nERROR: sequences of basis labels are not identical for each orbital subsection!"
      print "basis labels for subsection # 1 :",basis_labels
      print "basis labels for subsection #",i+1,":",basis_labels_allsubsection[i]
      sys.exit(0)

# remove numbers of basis labels, e.g.: 1S -> S, 10PX -> PX, etc...
  for i in range(nbasis):
    object_matched = re.search ("[SPDFGH].*", basis_labels[i])
    if object_matched:
       object = object_matched.group()
       basis_labels_without_number.append (object)

# check dimensions
  if (len(basis_labels_without_number) != nbasis):
    print "\nERROR: dimension of basis_labels_without_number =",len(basis_labels_without_number)," is different from nbasis =",nbasis
    sys.exit(0)

#  print "\nbasis_labels=",basis_labels
#  print "basis_labels_without_number=",basis_labels_without_number

# basis lables for each unique center type
  for cent_type in range(nctype):
    for i in range(nbasis):
      if (int(basis_atom_indexes[i]) == int(cent_type_to_first_cent [cent_type])):
        basis_labels_by_cent_type.append (basis_labels[i])
        basis_labels_without_number_by_cent_type.append (basis_labels_without_number[i])

#  print "basis_labels_by_cent_type=",basis_labels_by_cent_type

# ========================================================
def sort_basis ():
  "put basis function in champ order."
#  global orbital_coefficients

  print "sorting basis functions in standard CHAMP order..."

# define champ basis order
  basis_labels_to_champ_order = {'S':1,'PX':2,'PY':3,'PZ':4,'D0':5,'D+2':6,'D-2':7,'D+1':8,'D-1':9}

#  for i in range(nbasis):
#    print "basis function # ",i+1
#    print "basis_atom_index=",basis_atom_indexes[i]
#    print "basis_label=",basis_labels[i]
#    print "basis_labels_to_champ_order=",basis_labels_to_champ_order [basis_labels[i][1:]]

# sort basis functions
  for i in range(nbasis):
    for j in range(i+1,nbasis):
      if (basis_atom_indexes[i] == basis_atom_indexes[j]): 
        if (basis_labels_to_champ_order [basis_labels_without_number[j]] < basis_labels_to_champ_order[basis_labels_without_number[i]]):
#         reorder arrays
          temp = basis_labels [j]
          del basis_labels [j]
          basis_labels.insert (i, temp)

          temp = basis_labels_without_number [j]
          del basis_labels_without_number [j]
          basis_labels_without_number.insert (i, temp)

          temp = basis_atom_indexes [j]
          del basis_atom_indexes [j]
          basis_atom_indexes.insert (i, temp)

#          basis_labels [i], basis_labels [j] = basis_labels [j], basis_labels [i]
#          basis_labels_without_number [i], basis_labels_without_number [j] = basis_labels_without_number [j], basis_labels_without_number [i]
#          basis_atom_indexes [i], basis_atom_indexes [j] = basis_atom_indexes [j], basis_atom_indexes [i]

          for k in range(norb):
            temp = orbital_coefficients[k][j]
            del orbital_coefficients[k][j]
            orbital_coefficients[k].insert (i, temp)
#            orbital_coefficients[k][i], orbital_coefficients[k][j] = orbital_coefficients[k][j], orbital_coefficients[k][i]

#  print "orbital_coefficients[0]=",orbital_coefficients[0]
#  print "basis_atom_indexes=",basis_atom_indexes
#  print "basis_labels=",basis_labels

# ========================================================
def calculate_parameter_number ():
  "calculate the number of parameters."
  global nparm
  nparm = nctype*nparma + nparmb + nctype*nparmc + nparmcsf

# ========================================================
def check_electron_number ():
  """check number of electrons."""
  global nelec
  global total_charge
  global total_nuclear_charge
  total_nuclear_charge = sum(nuclear_charges)
  nelec_check = total_nuclear_charge - total_charge
  if (nelec_check != nelec):
    print "\nERROR: electron number inconsistent with nuclear charges and total charge!"
    sys.exit(0)

# ========================================================
# main code
# ========================================================

print "GAUSSIAN -> CHAMP convertion script"
print "WARNING: for now, it works only for single-determinant calculations!"

# check argument number
if len(sys.argv) <= 1:
  print "\nERROR: missing arguments"
  help_menu()

# command line arguments
#file_input_string = sys.argv[1]
#file_output_string = sys.argv[2]

# command line arguments
options, extra = getopt.getopt(sys.argv[1:],'i:o:hs')

for opt, val in options:
  if opt == "-h":
     help_menu()
  elif opt == "-i":
   file_input_string = val
  elif opt == "-o":
   file_output_string = val
  elif opt == "-s":
   l_sort_basis = True
  else:
     print "\nERROR: unknown option:", opt, val

# open and read input file
file_input = open(file_input_string,'r')
lines = file_input.readlines()

# read
sys.stdout.write('reading GAUSSIAN file >'+ str(file_input_string) +'< ... ')
read_total_charge ()
read_geometry ()
read_electron_number ()
read_energy ()
read_basis ()
read_orbitals ()

sys.stdout.write('done\n')

# sort basis in standard CHAMP order?
if (l_sort_basis):
  sort_basis ()


# close input file
file_input.close()

# calculate
calculate_parameter_number ()

# check
check_electron_number()

# default orbital occupation in determinant
for idet in range(ndet):
  det_up_orb_lab.append (range(1,nup+1))
  det_dn_orb_lab.append (range(1,ndn+1))

# write CHAMP input file in new format
sys.stdout.write('writing CHAMP file >'+ str(file_output_string) +'< ... ')
file_output = open(file_output_string,'w')
file_output.write('control\n')
file_output.write(' seed=1837465927472523\n')
file_output.write(' etrial=%.1f' %energy+ '\n')
file_output.write(' nstep=100 nblk=10 nblkeq=1 nconf=100\n')
file_output.write(' tau=0.01\n')
file_output.write('end\n\n')

file_output.write('nuclei\n')
file_output.write(' geometry\n')
for i in range(ncent):
  file_output.write('  ' + str(iwctype[i]) + ' %.1f' %nuclear_charges[i])
  for j in range(3):
    file_output.write('  %0.8f' % cent[i][j] + '  ')
  file_output.write('\n')
file_output.write(' end\n')
file_output.write('end\n\n')

file_output.write('wavefunction\n')
file_output.write(' nelec=' + str(nelec)+' nup='+str(nup)+'\n')
file_output.write('end\n\n')

file_output.write('basis\n')
file_output.write(' basis_functions\n')
for i in range(nctype):
  file_output.write('%d' %(i+1) + '\n')
  for j in range(len(basis_labels_without_number_by_cent_type)):
    file_output.write('?'+basis_labels_without_number_by_cent_type[j] + '   %0.8f' % 0. + '\n' )
file_output.write(' end\n')
file_output.write('end\n\n')

file_output.write('orbitals\n')
file_output.write(' coefficients')
for i in range(norb):
  file_output.write('\n')
  for j in range(nbasis):
    file_output.write(' %0.8e' % orbital_coefficients[i][j] + '  ')
file_output.write('\n end\n')
if (orbital_symmetry):
  file_output.write('symmetry')
  for i in range(norb):
    file_output.write(' '+orbital_symmetries[i])
  file_output.write(' end\n')
file_output.write('end\n\n')

file_output.write('csfs\n')
file_output.write(' determinants\n')
for idet in range(ndet):
  for i in range(nup):
    file_output.write('%d' % det_up_orb_lab[idet][i] + ' ')
  file_output.write(' ')
  for i in range(ndn):
    file_output.write('%d' % det_dn_orb_lab[idet][i] + ' ')
file_output.write('\n')
file_output.write(' end\n')
file_output.write('end\n\n')

file_output.write('jastrow\n')
file_output.write(' ijas=4 isc=4\n')
file_output.write(' scalek=0.8\n')
file_output.write(' parameters\n')
for i in range(nctype):
  file_output.write('0. 0. 0. 0. 0. 0. (a(iparmj),iparmj=1,nparma)\n')
file_output.write('0.5 1. 0. 0. 0. 0. (b(iparmj),iparmj=1,nparmb)\n')
for i in range(nctype):
  file_output.write('0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. (c(iparmj),iparmj=1,nparmc)\n')
file_output.write(' end\n')
file_output.write('end\n\n')

file_output.write('optimization\n')
file_output.write(' parameters jastrow end\n')
file_output.write('end\n')

# close output file
file_output.close()
sys.stdout.write('done\n')
sys.stdout.write('WARNING: basis information needs to be completed in CHAMP file!\n')



#########################################################################################"
# write CHAMP input file (old format)
sys.exit(0)

sys.stdout.write('writing CHAMP file >'+ str(file_output_string) +'< ... ')
file_output = open(file_output_string,'w')
file_output.write('your title                               title\n')
file_output.write('1837465927472523                         irn\n')
file_output.write('0 1                                      iperiodic,ibasis\n')
file_output.write('0.5  %.1f' %energy + '  \'  Hartrees\'                hb,etrial,eunit\n')
file_output.write('1000  10   1   100   0                   nstep,nblk,nblkeq,nconf,nconf_new\n')
file_output.write('0    0    1    -2                        idump,irstar,isite,ipr\n')
file_output.write('6  1.  5.  1.  1.                        imetro delta,deltar,deltat fbias\n')
file_output.write('2 1 1 1 1 0 0 0 0                        idmc,ipq,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e\n')
file_output.write('50  .01                                  nfprod,tau\n')
file_output.write('0  -3   1  0                             nloc,numr,nforce,nefp\n')
file_output.write(str(nelec)+' '+str(nup)+'                                     nelec,nup\n\n')

file_output.write('\'* Geometry section\'\n')
file_output.write('3                                        ndim\n')
file_output.write(str(nctype)+' '+str(ncent)+'                                      nctype,ncent\n')
for i in range(ncent):
  file_output.write(str(iwctype[i])+' ')
file_output.write(' (iwctype(i),i=1,ncent)\n')
for i in range(nctype):
  file_output.write('%.1f' %nuclear_charges_unique[i] + ' ')
file_output.write('(znuc(i),i=1,nctype)')
for i in range(ncent):
  file_output.write('\n')
  for j in range(3):
    file_output.write('  %0.8f' % cent[i][j] + '  ')
  file_output.write(str(i+1))
file_output.write('  ((cent(k,i),k=1,3),i=1,ncent)\n')

file_output.write('\n')
file_output.write('\'* Determinantal section\'\n')
file_output.write('0 0                                      inum_orb,iorb_used\n')
file_output.write(str(ndet)+' '+str(nbasis)+' '+str(norb)+'                                    ndet,nbasis,norb\n')
for i in range(nctype):
  file_output.write('0   0  0 0 0   0  0 0 0  0 0 0 0 0   0 0 0 0  0 0 0 0 0  0 0 0 0 0 0 0   0 0 0 0  0 0 0 0 0  0 0 0 0 0 0 0   0 0 0 0 0 0 0 0 0  0  0 0 0  0 0 0 0 0 n1s,n2s,n2px,n2py,np2z,...,n5g,sa,pa,da')
for i in range(norb):
  file_output.write('\n')
  for j in range(nbasis):
    file_output.write('%0.8e' % orbital_coefficients[i][j] + '  ')
  if i == 0:
    file_output.write('((coef(ibasis,iorb),ibasis=1,nbasis),iorb=1,norb)')
file_output.write('\n')
for i in range(nbasis):
  file_output.write('%0.8f' % 0. + ' ')
file_output.write('(zex(i),i=1,nbasis)\n')
for idet in range(ndet):
  for i in range(nup):
    file_output.write('%d' % det_up_orb_lab[idet][i] + ' ')
  file_output.write(' ')
  for i in range(ndn):
    file_output.write('%d' % det_dn_orb_lab[idet][i] + ' ')
file_output.write('(iworbd(j,idet),j=1,nelec)\n')
file_output.write(str(ncsf)+' ncsf\n')
for i in range(ncsf):
  file_output.write('%0.8f' % csf_coef[i] + ' ')
file_output.write(' (csf_coef(icsf),icsf=1,ncsf)\n')
for i in range(ncsf):
  file_output.write('%d' % ndet_in_csf[i] + ' ')
file_output.write('(ndet_in_csf(icsf),icsf=1,ncsf)\n')
for i in range(ncsf):
  for j in range(ndet_in_csf[i]):
    file_output.write('%d' %iwdet_in_csf[i][j] + ' ')
  file_output.write('(iwdet_in_csf(idet_in_csf,%d'%(i+1)+'),idet_in_csf=1,ndet_in_csf(%d'%(i+1)+'))\n')
  for j in range(ndet_in_csf[i]):
    file_output.write('%f' %cdet_csf[i][j] + ' ')
  file_output.write('(cdet_csf(idet_in_csf,%d'%(i+1)+'),idet_in_csf=1,ndet_in_csf(%d'%(i+1)+'))\n')
  
file_output.write('\n')
file_output.write('\'* Jastrow section\'\n')
file_output.write('1             ianalyt_lap\n')
file_output.write('4 4 1 1 5 0   ijas,isc,nspin1,nspin2,nord,ifock\n')
file_output.write('5 5 5         norda,nordb,nordc\n')
file_output.write('0.8 0. scalek,a21\n')
for i in range(nctype):
  file_output.write('0. 0. 0. 0. 0. 0. (a(iparmj),iparmj=1,nparma)\n')
file_output.write('0.5 1. 0. 0. 0. 0. (b(iparmj),iparmj=1,nparmb)\n')
for i in range(nctype):
  file_output.write('0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. (c(iparmj),iparmj=1,nparmc)\n')


file_output.write('\n')
file_output.write('\'* Optimization section\'\n')
file_output.write('100 1000 1.d-6 0. 1.d-3     nopt_iter,nblk_max,add_diag(1),p_var,tol_energy\n')
file_output.write('1000 '+str(nparm)+' 1 1 5 1000 21101 1 NDATA,NPARM,icusp,icusp2,NSIG,NCALLS,iopt,ipr\n')
file_output.write('0 0 0 0 i3body,irewgt,iaver,istrech\n')
file_output.write('0 0 0 0 0 0 0 0 0 0 ipos,idcds,idcdr,idcdt,id2cds,id2cdr,id2cdt,idbds,idbdr,idbdt\n')
for i in range(norb):
  file_output.write('%d' % 1 + ' ')
file_output.write('(lo(iorb),iorb=1,norb)\n')
file_output.write('0 ')
for i in range(nctype):
  file_output.write(str(nparma) + ' ')
file_output.write(str(nparmb) + ' ')
for i in range(nctype):
  file_output.write(str(nparmc) + ' ')
file_output.write('0 '+str(nparmcsf)+' 0 0  nparml,nparma,nparmb,nparmc,nparmf,nparmcsf,nparms,nparmg\n')
file_output.write('    (iworb(iparm),iwbasi(iparm),iparm=1,nlarml)\n')
file_output.write('    (iwbase(iparm),iparm=1,nparm-nparml)\n')
file_output.write('  ')
for i in range(2,nparmcsf+2):
  file_output.write(str(i) + ' ')
file_output.write('(iwcsf(iparm),iparm=1,nparmcsf)\n')
for i in range(nctype):
  file_output.write('    3 4 5 6 (iwjasa(iparm),iparm=1,nparma)\n')
file_output.write('2 3 4 5 6 (iwjasb(iparm),iparm=1,nparmb)\n')
for i in range(nctype):
  file_output.write('    3   5   7 8 9    11    13 14 15 16 17 18    20 21    23 (iwjasc(iparm),iparm=1,nparmc)\n')
file_output.write('0 0       necn,nebase\n')
file_output.write('          ((ieorb(j,i),iebasi(j,i),j=1,2),i=1,necn)\n')
file_output.write('          ((iebase(j,i),j=1,2),i=1,nebase)\n')
for i in range(norb):
  file_output.write('%d' % 0 + ' ')
file_output.write('(ipivot(j),j=1,norb)\n')
file_output.write('%.1f' % energy + ' eave\n')
file_output.write('1.d-6 5. 1 15 4 pmarquardt,tau,noutput,nstep,ibold\n')
file_output.write('T F analytic,cholesky\n')
file_output.write('end\n')

if (not l_sort_basis):
  file_output.write('\n')
  file_output.write('basis\n')
  file_output.write(' basis_functions\n')
  for i in range(nctype):
    file_output.write('  %d' %(i+1) + ' ')
    for j in range(len(basis_labels_without_number_by_cent_type)):
      file_output.write(basis_labels_without_number_by_cent_type[j] + ' ' )
    file_output.write('\n')
  file_output.write(' end\n')
  file_output.write('end\n')

file_output.write('\n')
file_output.write('optimization\n')
file_output.write(' parameters jastrow end\n')
file_output.write(' method = linear\n')
file_output.write('end\n')

# close output file
file_output.close()
sys.stdout.write('done\n')
sys.stdout.write('WARNING: basis information needs to be entered in CHAMP file!\n')
