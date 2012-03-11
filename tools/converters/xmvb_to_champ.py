#!/usr/bin/python
#!/usr/local/bin/python2.4

###################################################################################
#
# XMVB to CHAMP convertion script
# Author: J. Toulouse
# Created: 02 Mar 2009
# Last modified: 30 Apr 2011
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

angstrom_to_bohr = 1.889716165

nelec = 0
nup = 0
ndn = 0

nbasis = 0

norb = 0
orbital_indexes = []
orbital_occupations = []
orbital_eigenvalues = []
orbital_coefficients = []
orbital_symmetries = []
det_up_orb_lab = []
det_dn_orb_lab = []

ndet = 0
ncsf = 0
csf_coef = []
ndet_in_csf = [1]
structures = []
structures_initial = []
determinants = []
det_coef = []
cdet_csf = [[1.]]
det_label_in_csf = []
det_coef_in_csf = []

energy = 0.

# ========================================================
def help_menu():
# ========================================================
# print help menu
# --------------------------------------------------------
  print "USAGE: %s" % sys.argv[0], "-i xmvb_output_file -o champ_input_file"
  print "the script overwrites orbital and CSF information in champ_input_file"
  print "options: -h: print this menu"
  print "options: -d file.xdat: get orbitals (including virtual orbitals) from xdat file"
  print "options: -c: reconstruct list of determinants from the VB structures instead of reading them from xmvb_output_file"
  
  sys.exit(0)

# ========================================================
def not_empty_string (s):
  "filter not empty strings"
  return s != ''
  sys.exit(0)

# ========================================================
def read_energy ():
  """read the calculated energy."""
  global energy
  found = False

  for line in lines:
    if re.search ("Total Energy:", line):
      object_splitted = string.split (line)
      energy = float(object_splitted [2])
#      print "\nenergy=", energy,
      found = True
      break

  if not found:
    print "\nERROR: energy not found"
    sys.exit(0)

# ========================================================
def read_structure_number ():
  "read the number of structures."
  global ncsf
  found = False

  for line in lines:
    if re.search ("Number of Structures:", line):
      object_splitted = string.split (line)
      ncsf = int(object_splitted [3])
#      print "\nNumber of Structures=",ncsf
      found = True
      break

  if not found:
    print "\nERROR: Number of Structures not found"
    sys.exit(0)

# ========================================================
def read_structures_initial ():
  "read VB structures printed initially in the output without the coefficients."
  global ncsf
  current_structure = []
  found = False
  structures_initial_nb = 0

  for i in range(len(lines)):

#   begin of structure section
    if re.search ("The following structures are used in calculation:", lines[i]):

      found = True
      j = i + 1
      while j <len(lines): 

#        print "\nline=",  lines[j]

#       read first line of new structure
#       example:   1 *****    1  1  3  3
        if re.search ("^\s+\d+\s+[\*]+\s+(\s+\d+)+$", lines[j]):

          structures_initial_nb = structures_initial_nb + 1
          if (structures_initial_nb >= 2):
#            print "current_structure=",current_structure
            structures_initial.append (current_structure)
            current_structure = []
          
#          print "\nnew structure matched"
          line_splitted = filter (not_empty_string,string.split(lines[j],' '))
          for o in range(2,len(line_splitted)):
            current_structure.append(string.replace(line_splitted[o],'\n',''))


#       read continuation line of structure
#       example:                          8  8  9 10
        elif re.search ("^(\s+\d+)+$",lines[j]):
          line_splitted = filter (not_empty_string,string.split(lines[j],' '))
          for o in range(0,len(line_splitted)):
            current_structure.append(string.replace(line_splitted[o],'\n',''))

#       "Nuclear Repulsion Energy" marks the end of the initial structure section
        elif re.search ("Nuclear Repulsion Energy", lines[j]):
          structures_initial.append (current_structure)
          break

        j = j + 1
      break
       
  if not found:
    print "\nWARNING: no initial VB structures found"
    sys.exit(0)

#  print "\nstructures_initial_nb=",structures_initial_nb
#  print "structures_initial=",structures_initial

  if (structures_initial_nb != ncsf):
    print "\nERROR: mismatch in number of initial structures =",structures_initial_nb
    sys.exit(0)

# ========================================================
def read_structures ():
  "read VB structures with coefficients."
  global ncsf
  current_structure = []
  found = False
  structures_nb = 0

  for i in range(len(lines)):

#   begin of structure section
    if re.search ("COEFFICIENTS OF STRUCTURES", lines[i]):

      found = True
      j = i + 1
      while j <len(lines): 

#        print "\nline=",  lines[j]

#       read first line of new structure
#       example: 1     0.73801  ******    1  1  2  2  3  3  4  4  5  5  6  6  7  7
        if re.search ("^\s+\d+\s+[\d\.\-]+\s+[\*]+\s+(\s+\d+)+$", lines[j]):

          structures_nb = structures_nb + 1
          if (structures_nb >= 2):
#            print "current_structure=",current_structure
            structures.append (current_structure)
            current_structure = []
          
#          print "\nnew structure matched"
          line_splitted = filter (not_empty_string,string.split(lines[j],' '))
          csf_coef.append (float(line_splitted[1]))
          for o in range(3,len(line_splitted)):
            current_structure.append(string.replace(line_splitted[o],'\n',''))


#       read continuation line of structure
#       example:                          8  8  9 10
        elif re.search ("^(\s+\d+)+$",lines[j]):
          line_splitted = filter (not_empty_string,string.split(lines[j],' '))
          for o in range(0,len(line_splitted)):
            current_structure.append(string.replace(line_splitted[o],'\n',''))

#       COEFFICIENTS OF DETERMINANTS marks the end of the structure section
        elif re.search ("COEFFICIENTS OF DETERMINANTS", lines[j]):
          structures.append (current_structure)
          break

        j = j + 1
      break
       
  if not found:
#   if only one structure, take it for the initial section and assume coef to be 1
    if ncsf == 1:
       structures_nb = 1
       structures.append (structures_initial[0])
       csf_coef.append (1.)
    else:
      print "\nWARNING: VB structures not found"
      sys.exit(0)

#  print "csf_coef=",csf_coef
#  print "structures=",structures

  if (structures_nb != ncsf):
    print "\nERROR: mismatch in number of structures =",structures_nb
    sys.exit(0)

# check that the structures are all different
  for csf_i in range(ncsf):
    for csf_j in range(csf_i+1,ncsf):
      if structures[csf_i] == structures[csf_j]:
        print "\nERROR: structures ",csf_i+1, " and ",csf_j+1," are identical."
        sys.exit(0)
    

# ========================================================
def read_determinants ():
  "read determinants."
  global nelec, nup, ndn
  global ndet
  spin_up_spin_down = []
  current_determinant = []
  found = False

#              ******  COEFFICIENTS OF DETERMINANTS ****** 
#
#                            a  a  a  a  a  a  a  a  a  b  b  b  b  b
#                            b  b  b  b
# 
#   1     0.50791  ******    1  2  3  4  5  6  7  8  9  1  2  3  4  5
#                            6  7  8 10
#   2     0.50791  ******    1  2  3  4  5  6  7  8 10  1  2  3  4  5
#                            6  7  8  9
#   3     0.27754  ******    1  2 11 12 13 14 15 16 17  1  2 11 12 13
#                           14 15 16 17
#   4     0.27754  ******    1  2 18 19 20 21 22 23 24  1  2 18 19 20
#                           21 22 23 24

  for i in range(len(lines)):

#   begin of structure section
    if re.search ("COEFFICIENTS OF DETERMINANTS", lines[i]):

      found = True
      j = i + 1
      while j <len(lines): 

#        print "\nline=",  lines[j]

#       read spin information
#       example: a  a  a  a  a  a  a  a  a  b  b  b  b  b
        if re.search ("^(\s+[ab])+$", lines[j]):
          line_splitted = filter (not_empty_string,string.split(lines[j],' '))
          for o in range(len(line_splitted)):
            spin_up_spin_down.append(string.replace(line_splitted[o],'\n',''))

#       read first line of new determinant
#       example: 1     0.50791  ******    1  2  3  4  5  6  7  8  9  1  2  3  4  5
        elif re.search ("^\s+\d+\s+[\d\.\-]+\s+[\*]+\s+(\s+\d+)+$",lines[j]):

          ndet = ndet + 1
          if (ndet >= 2):
#            print "current_determinant=",current_determinant
            determinants.append (current_determinant)
            current_determinant = []

          line_splitted = filter (not_empty_string,string.split(lines[j],' '))
          det_coef.append (float(line_splitted[1]))
          for o in range(3,len(line_splitted)):
            current_determinant.append(string.replace(line_splitted[o],'\n',''))

#       read continuation line of determinant
#       example:                          6  7  8 10
        elif re.search ("^(\s+\d+)+$",lines[j]):
          line_splitted = filter (not_empty_string,string.split(lines[j],' '))
          for o in range(0,len(line_splitted)):
            current_determinant.append(string.replace(line_splitted[o],'\n',''))

#       WEIGHTS OF STRUCTURES marks the end of the determinant section
        elif re.search ("WEIGHTS OF STRUCTURES", lines[j]):
          determinants.append (current_determinant)
          break

        j = j + 1
      break
       
  if not found:
#   if only one structure, construct determinants
    if ncsf == 1:
      print "\nWARNING: no determinants found, constructing determinants for single structure by singlet couplings (to be checked)"
      determinant_up =[]
      determinant_dn =[]
      coupled_orbital_1=[]
      coupled_orbital_2=[]
      nup = 0
      ndn = 0
      o=0
      while o < len(structures[0]):
        orb_up = structures[0][o]
        determinant_up.append (orb_up)
        nup = nup + 1
        o = o + 1
        if o >= len(structures[0]):
          break
        orb_dn = structures[0][o]
        if (orb_dn != orb_up):
          coupled_orbital_1.append(orb_up)
          coupled_orbital_2.append(orb_dn)
        determinant_dn.append (orb_dn)
        ndn = ndn + 1
        o = o + 1
      determinants.append (determinant_up+determinant_dn)
      det_coef.append (1.)
      nelec = nup+ndn
#    construct all determinants by spin-singlet couplings of unpaired electron
      for pair in range(len(coupled_orbital_1)):
        ndet=len(determinants)
        for det in range(ndet):
          new_determinant=[]
          for orb in determinants[det]:
             if (orb == coupled_orbital_1[pair]):
               new_determinant.append (coupled_orbital_2[pair])
             elif (orb == coupled_orbital_2[pair]):
               new_determinant.append (coupled_orbital_1[pair])
             else:
               new_determinant.append (str(orb))
          determinants.append (new_determinant)
          det_coef.append (1.)
    else:
      print "\nERROR: determinants not found"
      sys.exit(0)

# calculate numbers of spin-up and spin-down electrons
  if found:
#    print "spin_up_spin_down=",spin_up_spin_down
    nup = 0
    ndn = 0
    for s in spin_up_spin_down:
      if (s == 'a'):
        nup = nup + 1
      elif (s == 'b'):
        ndn = ndn + 1
    nelec = nup+ndn

#  print "nup=",nup
#  print "ndn=",ndn
#  print "nelec=",nelec
#  print "determinants=",determinants
#  print "det_coef=",det_coef

# ========================================================
def read_orbitals ():
  "read molecular orbitals."
  global norb
  global nbasis
  global orbital_coefficients
  orbital_coefficients_subsection = []
  found = False
  j = 0
  subsection = 0
  begin_subsection = False
  for i in range(len(lines)):

#   begin of orbital section
    if re.search ("OPTIMIZED ORBITALS", lines[i]):
      found = True
      j = i + 1
      while j <len(lines): 

#       read orbital indexes, occupations and eigenvalues
#       example: 1          2          3          4          5
        if re.search ("^(\s+\d+)+$", lines[j]):
          
          orbital_indexes.extend (string.split(lines[j]))
#          print "line ",j,", orbital indexes:",lines[j]

          subsection = subsection + 1
          if (subsection >= 2):
            orbital_coefficients.extend ((transpose(orbital_coefficients_subsection)).tolist())
            orbital_coefficients_subsection = []

#       read normal line of orbital coefficients
#       example: 1     0.207568   0.208271   0.018863   0.000000   0.000000
        elif re.search ("^\s+\d+\s+(\s+[\d\-]+)+",lines[j]):
#          print "line ",j,", orbital coefficient:",lines[j]
          line_splitted = filter (not_empty_string,string.split(lines[j],' '))
#          print "line_splitted=",line_splitted
#          print "line_splitted[1:]=",line_splitted[1:]
          orbital_coefficients_subsection.append (map(float,map(string.strip,line_splitted[1:])))

#       Dipole moment marks the end of the orbital section
        elif re.search ("Dipole moment", lines[j]):
          orbital_coefficients.extend ((transpose(orbital_coefficients_subsection)).tolist())
          break

        j = j + 1
      break
       
  if not found:
    print "\nERROR: orbital coefficients not found in xmvb output file"
    sys.exit(0)

#  print "orbitals read"

  norb = len(orbital_indexes)
  nbasis = len(orbital_coefficients[0])
#  print "norb=",norb
#  print "orbital_indexes=", orbital_indexes
#  print "orbital_coefficients=",orbital_coefficients

# ========================================================
def read_orbitals_from_xdat ():
  "read molecular orbitals from xdat file."
  global norb
  global nbasis
  global orbital_coefficients
  orbital_coefficients_subsection = []
  found = False
  j = 0
  subsection = 0
  begin_subsection = False
  orbital_indexes = []
  orbital_coefficients = []

  for i in range(len(lines)):

#   begin of orbital section
    if re.search ("ORBITALS IN TERMS OF BASIS FUNCTIONS", lines[i]):
      found = True
      j = i + 1
      while j <len(lines): 

#       read orbital indexes
#       example: 1          2          3          4          5
        if re.search ("^(\s+\d+)+$", lines[j]):
          
          orbital_indexes.extend (string.split(lines[j]))
#          print "line ",j,", orbital indexes:",lines[j]

          subsection = subsection + 1
          if (subsection >= 2):
            orbital_coefficients.extend ((transpose(orbital_coefficients_subsection)).tolist())
            orbital_coefficients_subsection = []

#       read normal line of orbital coefficients
#       example: 1     0.207568   0.208271   0.018863   0.000000   0.000000
        elif re.search ("^\s+\d+\s+(\s+[\d\-]+)+",lines[j]):
#          print "line ",j,", orbital coefficient:",lines[j]
          line_splitted = filter (not_empty_string,string.split(lines[j],' '))
#          print "line_splitted=",line_splitted
#          print "line_splitted[1:]=",line_splitted[1:]
          orbital_coefficients_subsection.append (map(float,map(string.strip,line_splitted[1:])))

#       Dipole moment marks the end of the orbital section
#        elif re.search ("Dipole moment", lines[j]):
#          orbital_coefficients.extend ((transpose(orbital_coefficients_subsection)).tolist())
#          break

        j = j + 1
      break
  orbital_coefficients.extend ((transpose(orbital_coefficients_subsection)).tolist())
       
  if not found:
    print "\nERROR: orbital coefficients not found in xdat file"
    sys.exit(0)

#  print "orbitals read"

  norb = len(orbital_indexes)
  nbasis = len(orbital_coefficients[0])
#  print "norb=",norb
#  print "nbasis=",nbasis
#  print "orbital_indexes=", orbital_indexes
#  print "orbital_coefficients=",orbital_coefficients

# ========================================================
def determine_orbitals_sym ():
  "determine orbitals symmetry classes."
  global orbital_symmetries

  
# determine occupied orbitals
  orbital_occupied = zeros([norb], Int)
  for det_i in range(ndet):
    for elec_i in range(nelec):
       orbital_occupied[int(determinants[det_i][elec_i])-1]=1
  #print "orbital_occupied=",orbital_occupied

# determine orbital "symmetry" classes
  orbital_sym_classes = []
  orb_sym_nb = 0
  for orb_i in range(norb):
    orbital_sym_class = zeros([nbasis], Int)
    for bas_i in range(nbasis):
       if (not orbital_coefficients[orb_i][bas_i] == 0):
         orbital_sym_class[bas_i]=1
    #print "orbital_sym_class=",orbital_sym_class

    found = False
    for orb_sym in range(orb_sym_nb):
      found = True
      for bas_i in range(nbasis):
        if (orbital_occupied[orb_i] == 1):
          if (not orbital_sym_class[bas_i] == orbital_sym_classes [orb_sym][bas_i]):
            found = False
            break
#       for virtual orbital, check only that orbital_sym_class is "included" in one the orbital_sym_classes
        elif (orbital_occupied[orb_i] == 0):
          if (orbital_sym_class[bas_i] == 1 and not orbital_sym_class[bas_i] == orbital_sym_classes [orb_sym][bas_i]):
            found = False
            break
      if (found):
        break
    if (found):
#      print "existing class ",orb_sym+1
      orbital_symmetries.append (orb_sym+1)
    else:
      orb_sym_nb = orb_sym_nb + 1
#      print "new class ",orb_sym_nb
      orbital_symmetries.append (orb_sym_nb)
      orbital_sym_classes.append (orbital_sym_class)


# filter for symmetry classes that are alike (different only by one basis function) and suppress one                 #
  for orb_sym_i in range(orb_sym_nb-1):                                                                              #
    suppressed = False                                                                                               #
    for orb_sym_j in range(orb_sym_i+1,orb_sym_nb):                                                                  #
# if a class has been suppressed, the indexes has to be moved up                                                     #
       if (suppressed):                                                                                              #
          orb_sym_j = orb_sym_j - 1                                                                                  #
       tester = 0                                                                                                    #
       for bas_i in range(nbasis):                                                                                   #
          if (not orbital_sym_classes [orb_sym_i][bas_i] == orbital_sym_classes [orb_sym_j][bas_i]):                 #
             tester = tester+1                                                                                       #
# remember which orbital is to be suppressed (the one that has lesser basis functions)                               #
             if (orbital_sym_classes [orb_sym_i][bas_i] > orbital_sym_classes [orb_sym_j][bas_i]):                   #
                orb_sym_target = orb_sym_j                                                                           #
                orb_sym_new = orb_sym_i                                                                              #
             else:                                                                                                   #
                orb_sym_target = orb_sym_i                                                                           #
                orb_sym_new = orb_sym_j                                                                              #
# if classes are different by more than one basis : skip                                                             #
             if tester > 1:                                                                                          #
               break                                                                                                 #
          if tester > 1:                                                                                             #
            break                                                                                                    #
# if indeed a class is to be suppressed : modify -orbital_sym_classes (suppress the class)                           #
#                                                -orb_sym_nb (one less)                                              #
#                                                -orbital_symmetries (everyone beyond is shifted down)               #
       if tester == 1:                                                                                               #
          suppressed = True                                                                                          #
          orbital_sym_classes.pop(orb_sym_target)                                                                    #
          orb_sym_nb = orb_sym_nb - 1                                                                                #
          for orb_i in range(norb):                                                                                  #
             if orbital_symmetries[orb_i] == orb_sym_target+1:                                                       #
                orbital_symmetries[orb_i] = orb_sym_new+1                                                            #
             if orbital_symmetries[orb_i] > orb_sym_target+1:                                                        #
                orbital_symmetries[orb_i] = orbital_symmetries[orb_i] - 1                                            #
          print "WARNING: a symmetry class had to be removed (it was similar to another one)"                         #

 
# check that two symmetry classes do not share the same basis functions
# [optionnal] and that any basis function is in a symmetry class
  orbital_sym_class_sum = zeros([nbasis], Int)
  for orb_sym in range(orb_sym_nb):
#    print orb_sym+1,orbital_sym_classes[orb_sym]
    for bas_i in range(nbasis):
      orbital_sym_class_sum[bas_i] = orbital_sym_class_sum[bas_i] + orbital_sym_classes [orb_sym][bas_i]
  for bas_i in range(nbasis):
#    comment this out because some basis functions may not be used (e.g., some d functions)
#    if (orbital_sym_class_sum[bas_i] < 1):
#      print "\nWARNING: basis function #", bas_i+1," is not in any symmetry class!\n"
#      sys.exit(0)
    if (orbital_sym_class_sum[bas_i] > 1):
      print "\nERROR: basis function #", bas_i+1," is in more than one symmetry class!\n"
      sys.exit(0)

#  print "orbital_sym_class_sum=",orbital_sym_class_sum
       

# ========================================================
def dets_in_csfs ():
  "determine determinants in csfs."
  current_det_label_in_csf = []
  structures_orb_occ = zeros([ncsf,norb], Int)
  determinants_orb_occ = zeros([ndet,norb], Int)

# determine orbital occupation in structures
  for csf_i in range(ncsf):
    for elec_i in range(nelec):
      orb_i = int(structures[csf_i][elec_i])-1
      structures_orb_occ[csf_i][orb_i] = structures_orb_occ[csf_i][orb_i] + 1
#  print "structures_orb_occ=\n",structures_orb_occ

# determine orbital occupation in determinants
  for det_i in range(ndet):
    for elec_i in range(nelec):
      orb_i = int(determinants[det_i][elec_i])-1
      determinants_orb_occ[det_i][orb_i] = determinants_orb_occ[det_i][orb_i] + 1
#  print "determinants_orb_occ=\n",determinants_orb_occ


# determine determinants in each structure
  for csf_i in range(ncsf):
    if (csf_i >= 1): 
      det_label_in_csf.append (current_det_label_in_csf)
      current_det_label_in_csf = []
#    print "csf_i=",csf_i
#    print "structures_orb_occ[csf_i]=",structures_orb_occ[csf_i]
    for det_i in range(ndet):
      det_in_csf = True
      for orb_i in range(norb):
         if (not determinants_orb_occ[det_i][orb_i]  == structures_orb_occ[csf_i][orb_i]):
           det_in_csf = False
           break
      if (det_in_csf):
        current_det_label_in_csf.append (det_i+1)
#        print "determinant is in CSF"
  det_label_in_csf.append (current_det_label_in_csf)
#  print "det_label_in_csf=", det_label_in_csf

  det_used = []
  for i in range(ndet):
    det_used.append (0)

  current_det_coef_in_csf = []
  for csf_i in range(ncsf):
    if (csf_i >= 1):
      det_coef_in_csf.append (current_det_coef_in_csf)
      current_det_coef_in_csf = []
    for det_i in det_label_in_csf [csf_i]:
      det_used[det_i-1] = det_used[det_i-1] + 1
      current_det_coef_in_csf.append (det_coef[det_i-1]/csf_coef[csf_i])
  det_coef_in_csf.append (current_det_coef_in_csf)
#  print "det_coef_in_csf=",det_coef_in_csf

#  print "det_used=",det_used
  for det_i in range(ndet):
    if (det_used[det_i] == 0):
      print "WARNING: determinant # ",det_i+1," is not used\n"
    if (det_used[det_i] > 1):
      print "\nERROR: determinant # ",det_i+1," is used in more than one structure. The coefficients of determinants in each structure cannot be determined.\n"
      print "\nUse option '-c' to reconstruct the list of determinants from the structure.\n"
      sys.exit(0)

# ========================================================
def construct_determinants ():
  "construct determinants (and dets_in_csfs) if determinants not read"
  global nelec, nup, ndn
  global ndet

  print "\nWARNING: reconstruct list of determinants from the structures, the structures will only be approximately normalized (to be checked)"

# construct the list of non-sorted and non-unique composing determinants from the structures
  ndet = 0
  determinants_unsorted=[]
  for csf_i in range(ncsf):
    coupled_orbital_1=[]
    coupled_orbital_2=[]
    determinant_up =[]
    determinant_dn =[]
    determinants_in_csf=[]
    current_det_label_in_csf = []
    current_det_coef_in_csf=[]
    ndet = ndet + 1
    nup = 0
    ndn = 0
    o=0
    while o < len(structures[csf_i]):
      orb_up = structures[csf_i][o]
      determinant_up.append (orb_up)
      nup = nup + 1
      o = o + 1
      if o >= len(structures[csf_i]):
        break
      orb_dn = structures[csf_i][o]
      if (orb_dn != orb_up):
        coupled_orbital_1.append(orb_up)
        coupled_orbital_2.append(orb_dn)
      determinant_dn.append (orb_dn)
      ndn = ndn + 1
      o = o + 1
    determinants_in_csf.append (determinant_up+determinant_dn)
    current_det_label_in_csf.append (ndet)
    current_det_coef_in_csf.append (1.)
    
    nelec = nup+ndn
#  construct all determinants by spin-singlet couplings of unpaired electron
    for pair in range(len(coupled_orbital_1)):
      ndet_csf=len(determinants_in_csf)
      for det in range(ndet_csf):
        new_determinant=[]
        for orb in determinants_in_csf[det]:
           if (orb == coupled_orbital_1[pair]):
             new_determinant.append (coupled_orbital_2[pair])
           elif (orb == coupled_orbital_2[pair]):
             new_determinant.append (coupled_orbital_1[pair])
           else:
             new_determinant.append (str(orb))
        ndet = ndet + 1
        determinants_in_csf.append (new_determinant)
        current_det_label_in_csf.append (ndet)
        current_det_coef_in_csf.append (1.) # coefficients are always positive before sorting orbitals in determinants
    for det in range(len(determinants_in_csf)):
      determinants_unsorted.append (determinants_in_csf[det])
    det_label_in_csf.append (current_det_label_in_csf)
    det_coef_in_csf.append (current_det_coef_in_csf)


# sort orbitals in spin-up and spin-down determinants and keep track of sign changes on the coefficients
  for det_i in range(ndet):
#    print "determinant       : ",determinants_unsorted[det_i]
    det_up_and_sign = sort_and_sign(determinants_unsorted[det_i][:nup])
    det_dn_and_sign = sort_and_sign(determinants_unsorted[det_i][nup:])
    determinants.append (det_up_and_sign[0]+det_dn_and_sign[0])
#    print "determinant sorted: ",determinants[det_i]
#   change sign of the coefficients involving the current sorted determinant
    for csf_i in range(ncsf):
      for det_j in range(len(det_label_in_csf[csf_i])):
        det = det_label_in_csf[csf_i][det_j]
        if det-1 == det_i:
          det_coef_in_csf[csf_i][det_j]=det_coef_in_csf[csf_i][det_j]*det_up_and_sign[1]*det_dn_and_sign[1]
  
# keep only list of unique determinants
  det_i=1
  while det_i < len(determinants):
   det_unique = True
   for det_j in range(det_i):
     if determinants[det_i] == determinants[det_j]:
#       print "identical determinants:", det_i+1, det_j+1
       det_unique = False
       break
   if det_unique:
     det_i = det_i +1
   else:
     del determinants[det_i]
     for csf_i in range(ncsf):
       for det_k in range(len(det_label_in_csf[csf_i])):
         det = det_label_in_csf[csf_i][det_k]
         if det-1 == det_i:
           det_label_in_csf[csf_i][det_k]=det_j+1
         if det-1 > det_i:
           det_label_in_csf[csf_i][det_k]=det_label_in_csf[csf_i][det_k]-1
     
  ndet=len(determinants)

# approximate normalization of the determinant coefficients (assuming the determinants are orthonormal, which they are not)
  for csf_i in range(ncsf):
    norm = 0
    for det_i in range(len(det_coef_in_csf[csf_i])):
      norm = norm + det_coef_in_csf[csf_i][det_i]**2
    for det_i in range(len(det_coef_in_csf[csf_i])):
      det_coef_in_csf[csf_i][det_i] = det_coef_in_csf[csf_i][det_i]/sqrt(norm)
      
  
#  print "nup=",nup
#  print "ndn=",ndn
#  print "nelec=",nelec
#  print "determinants=",determinants
#  print "det_label_in_csf=",det_label_in_csf
#  print "det_coef_in_csf=",det_coef_in_csf
      
# ========================================================
def sort_and_sign (array):
  "sort array and determine the signature of the sort permutation."

  sign=1
#  print "original array=",array
  for i in range(len(array)):
    for j in range(i+1,len(array)):
      if float(array[j]) < float(array[i]):
        tempi=array[i]
        tempj=array[j]
        del array[i]
        array.insert(i,tempj)
        del array[j]
        array.insert(j,tempi)
        sign=sign*(-1)
  return(array,sign)
#  print "sorted array=",array
#  print "sign=",sign

# ========================================================
# main code
# ========================================================

print "XMVB -> CHAMP convertion script"
print "WARNING: not fully tested!"

# check argument number
if len(sys.argv) <= 1:
  print "\nERROR: missing arguments"
  help_menu()

# command line arguments
#file_input_string = sys.argv[1]
#file_output_string = sys.argv[2]

# command line arguments
options, extra = getopt.getopt(sys.argv[1:],'i:o:d:hsc')

l_xdat = False
l_construct_det = False

for opt, val in options:
  if opt == "-h":
     help_menu()
  elif opt == "-i":
   file_input_string = val
  elif opt == "-o":
   file_output_string = val
  elif opt == "-d":
   file_xdat_string = val
   l_xdat = True
  elif opt == "-c":
   l_construct_det = True
  else:
     print "\nERROR: unknown option:", opt, val

# read xmvb input file
file_input = open(file_input_string,'r')
sys.stdout.write('reading XMVB file >'+ str(file_input_string) +'< ... ')
lines = file_input.readlines()
read_energy ()
read_structure_number ()
read_structures_initial ()
read_structures ()
if not l_construct_det:
  read_determinants ()
read_orbitals ()
file_input.close()
sys.stdout.write('done\n')

# read xdat file
if l_xdat:
  file_xdat = open(file_xdat_string,'r')
  lines = file_xdat.readlines()
  sys.stdout.write('reading xdat file >'+ str(file_xdat_string) +'< ... ')
  read_orbitals_from_xdat ()
  sys.stdout.write('done\n')
  file_xdat.close()

determine_orbitals_sym();

# calculate
if not l_construct_det:
  dets_in_csfs ()
else:
  construct_determinants ()

# check

# default orbital occupation in determinant
for idet in range(ndet):
  det_up_orb_lab.append (range(1,nup+1))
  det_dn_orb_lab.append (range(1,ndn+1))

# read CHAMP input file
file_output = open(file_output_string,'r')
lines = file_output.readlines()
file_output.close()

# overwrite CHAMP input file
sys.stdout.write('overwriting CHAMP file >'+ str(file_output_string) +'< ... ')
file_output = open(file_output_string,'w')

skip_line = False

for l in range(len(lines)):
  line = lines[l]
  if (skip_line):
    if re.search ("^end", line):
      skip_line = False
    continue
    
  if re.search ("etrial", line):
    file_output.write(' etrial=%.1f' %energy+ '\n')
    continue

  if re.search ("\d+\s+ncsf", line):
    file_output.write(str(ncsf)+' ncsf\n')
    continue

  if re.search ("^\s*orbitals\s*$", line):
    skip_line = True
    file_output.write('orbitals\n')
    file_output.write(' coefficients')
    for i in range(norb):
      file_output.write('\n')
      for j in range(nbasis):
        file_output.write(' %0.8e' % orbital_coefficients[i][j] + '  ')
    file_output.write('\n end\n')
    file_output.write(' symmetry')
    for i in range(norb):
      file_output.write(' '+str(orbital_symmetries[i]))
    file_output.write(' end\n')
    file_output.write('end\n')
    continue

  if re.search ("^\s*csfs\s*$", line):
    skip_line = True
    file_output.write('csfs\n')
    file_output.write(' determinants\n')
    for idet in range(ndet):
      for i in range(nelec):
        file_output.write('%s' % determinants[idet][i] + ' ')
        if (i == nup -1):
          file_output.write('  ')
      file_output.write('\n')
    file_output.write(' end\n')
    file_output.write(' csf_coef')
    for i in range(ncsf):
      file_output.write(' %.8f' %csf_coef[i])
    file_output.write(' end\n')
    file_output.write(' dets_in_csfs\n')
    for csf_i in range(ncsf):
      for det_i in range(len(det_label_in_csf [csf_i])):
        file_output.write('%s ' %det_label_in_csf[csf_i][det_i])
      file_output.write('\n')
      for det_i in range(len(det_coef_in_csf [csf_i])):
        file_output.write('%.8f ' %det_coef_in_csf[csf_i][det_i])
      file_output.write('\n')
    file_output.write(' end\n')
    file_output.write(' vb_weights = true\n')
    file_output.write('end\n')
    continue

  file_output.write(line)

   
# close output file
file_output.close()
sys.stdout.write('done\n')

