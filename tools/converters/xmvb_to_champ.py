#!/usr/local/bin/python2.4
#!/usr/bin/python

###################################################################################
#
# XMVB to CHAMP convertion script
# Author: J. Toulouse
# Created: 02 Mar 2009
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
det_up_orb_lab = []
det_dn_orb_lab = []

ndet = 0
ncsf = 0
csf_coef = []
ndet_in_csf = [1]
structures = []
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
  print "options: -h: print this menu"
  
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
    print "ERROR: energy not found"
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
    print "ERROR: Number of Structures not found"
    sys.exit(0)

# ========================================================
def read_structures ():
  "read VB structures."
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
    print "ERROR: VB structures not found"
    sys.exit(0)

#  print "csf_coef=",csf_coef
#  print "structures=",structures

  if (structures_nb != ncsf):
    print "ERROR: mismatch in number of structures =",structures_nb
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
    print "ERROR: determinants not found"
    sys.exit(0)

#  print "spin_up_spin_down=",spin_up_spin_down

# calculate numbers of spin-up and spin-down electrons
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
    print "ERROR: orbital coefficients not found"
    sys.exit(0)

#  print "orbitals read"

  norb = len(orbital_indexes)
  nbasis = len(orbital_coefficients[0])
#  print "orbital_indexes=", orbital_indexes
#  print "orbital_coefficients=",orbital_coefficients

# ========================================================
def dets_in_csfs ():
  "determine determinants in csfs."
  current_det_label_in_csf = []

  for csf_i in range(ncsf):
    if (csf_i >= 1): 
      det_label_in_csf.append (current_det_label_in_csf)
      current_det_label_in_csf = []
#    print "csf_i=",csf_i
#    print "structure=",structures[csf_i]
    for det_i in range(ndet):
#      print "determinants=",determinants[det_i]
      det_in_csf = True
      for elec_i in range(nelec):
         orb_found = False
         for elec_j in range(nelec):
           if (determinants[det_i][elec_i]  == structures[csf_i][elec_j]):
             orb_found = True 
             break
         if (not orb_found):
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
      print "ERROR: determinant # ",det_i+1," is used in more than one structure. This case is not handled yet.\n"
      sys.exit(0)
      

# ========================================================
# main code
# ========================================================

print "XMVB -> CHAMP convertion script"
print "WARNING: not fully tested!"

# check argument number
if len(sys.argv) <= 1:
  print "ERROR: missing arguments"
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
  else:
     print "ERROR: unknown option:", opt, val

# open and read input file
file_input = open(file_input_string,'r')
lines = file_input.readlines()

# read
sys.stdout.write('reading XMVB file >'+ str(file_input_string) +'< ... ')
read_energy ()
read_structure_number ()
read_structures ()
read_determinants ()
read_orbitals ()
sys.stdout.write('done\n')

# close input file
file_input.close()

# calculate
dets_in_csfs ()

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

  if re.search ("orbitals", line):
    skip_line = True
    file_output.write('orbitals\n')
    file_output.write(' coefficients')
    for i in range(norb):
      file_output.write('\n')
      for j in range(nbasis):
        file_output.write(' %0.8e' % orbital_coefficients[i][j] + '  ')
    file_output.write('\n end\n')
    file_output.write('end\n')
    continue

  if re.search ("csfs", line):
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
    file_output.write('end\n')
    continue

  file_output.write(line)

   
# close output file
file_output.close()
sys.stdout.write('done\n')

