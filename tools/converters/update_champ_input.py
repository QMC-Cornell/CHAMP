#!/usr/bin/python
#!/usr/local/bin/python2.4

###################################################################################
#
# update_champ_input.py script for creating a new input with the optimized parameters from a previous run
# Author: J. Toulouse
# Created: 02 May 2011
# Last modified: 03 May 2011
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


# ========================================================
def help_menu():
# ========================================================
# print help menu
# --------------------------------------------------------
  print "USAGE: %s" % sys.argv[0], "-i champ_output_file -o champ_input_file"
  print "this script extracts the wave function parameters from champ_output_file and overwrites them in champ_input_file"
  print "options: -h: print this menu"
  print "options: -n 4: use the wave function at iteration 4"
  sys.exit(0)

# ========================================================
# main code
# ========================================================

#print "updating CHAMP input:"

# check argument number
if len(sys.argv) <= 1:
  print "\nERROR: missing arguments"
  help_menu()

# command line arguments
#file_input_string = sys.argv[1]
#file_output_string = sys.argv[2]

# command line arguments
options, extra = getopt.getopt(sys.argv[1:],'i:o:n:')

l_iteration_wanted = False

for opt, val in options:
  if opt == "-h":
     help_menu()
  elif opt == "-i":
   file_input_string = val
  elif opt == "-o":
   file_output_string = val
  elif opt == "-n":
   l_iteration_wanted = True
   iteration_wanted = int(val)
  else:
     print "\nERROR: unknown option:", opt, val

# reading champ output file
file_input = open(file_input_string,'r')
sys.stdout.write('reading CHAMP output file >'+ str(file_input_string) +'<\n')
lines = file_input.readlines()
found = False
file_input.close()
#sys.stdout.write('done\n')


# search for wave function
found_wf = False
found_wf_iter = False
found_wf_best = False

for i in range(len(lines)):

  if re.search ("SVN Revision:", lines[i]):
    revision_number = int(string.split((string.split (lines[i])[2]),',')[0])
    if (revision_number < 542):
      print "\nERROR: the output file has SVN Revision number = ",revision_number, ". This script requires SVN Revision number >= 542"
      sys.exit(0)

  if re.search ("Wave function for next iteration #", lines[i]):
    found_wf = True
    line_wf = i
    iter = int(string.split (lines[i])[6])
#    print "found wave function of iter=",iter
    if l_iteration_wanted:
      if iteration_wanted == iter:
        found_wf_iter = True
        line_wf_iter = i

  if re.search ("Best wave function:", lines[i]):
    found_wf_best = True
    line_wf_best = i
#    print "found best wave function"


if not l_iteration_wanted: 
  if found_wf_best:
    print "Use best wave function."
    line_wf = line_wf_best
  elif found_wf:
    print "Best wave function not found. Use the last found wave function of iteration # ", iter
  else:
    print "\nERROR: no wave function found."
    sys.exit(0)

if l_iteration_wanted: 
  if found_wf_iter:
    print "Use wave function of iteration # ",iteration_wanted
    line_wf = line_wf_iter
  else:
    print "\nERROR: wave function of iteration # ",iteration_wanted, " not found"
    sys.exit(0)

jastrow_read = False
csfs_read = False
orbitals_read = False
exponents_read = False
lines_jastrow=[]
lines_csfs=[]
lines_orbitals=[]
lines_exponents=[]

for i in range(line_wf,len(lines)):
  if (not jastrow_read) and re.search ("^jastrow$", lines[i]):
    j=i+1
    while j<len(lines): 
      if re.search (" parameters$", lines[j]):
        j=j+1
      if re.search ("\s*end\s*", lines[j]):
        break
      lines_jastrow.append(lines[j])
      jastrow_read = True
      j=j+1

  if (not csfs_read) and re.search ("^csfs$", lines[i]):
    j=i+1
    while j<len(lines): 
      if re.search (" csf_coef$", lines[j]):
        j=j+1
      if re.search ("\s*end\s*", lines[j]):
        break
      lines_csfs.append(lines[j])
      csfs_read = True
      j=j+1

  if (not orbitals_read) and re.search ("^orbitals$", lines[i]):
    j=i+1
    while j<len(lines): 
      if re.search (" coefficients$", lines[j]):
        j=j+1
      if re.search ("\s*end\s*", lines[j]):
        break
      lines_orbitals.append(lines[j])
      orbitals_read = True
      j=j+1

  if (not exponents_read) and re.search ("^basis$", lines[i]):
    j=i+1
    while j<len(lines): 
      if re.search (" basis_functions$", lines[j]):
        j=j+1
      if re.search ("\s*end\s*", lines[j]):
        break
      lines_exponents.append(lines[j])
      exponents_read = True
      j=j+1


if (not jastrow_read and not csfs_read and not orbitals_read and not exponents_read):
  print "\nERROR: no optimized parameters found"
  sys.exit(0)
 

# opening CHAMP input file
file_output = open(file_output_string,'r')
lines = file_output.readlines()
file_output.close()

# overwrite CHAMP input file
sys.stdout.write('overwriting CHAMP input file >'+ str(file_output_string) +'<\n')
file_output = open(file_output_string,'w')

skip_line=False
in_jastrow_menu=False

for i in range(len(lines)):
  if (skip_line):
    if re.search ("\s*end\s*", lines[i]):
      skip_line = False
    continue

  if jastrow_read:
    if re.search ("^\s*jastrow\s*$", lines[i]):
      in_jastrow_menu=True
    
    if in_jastrow_menu and re.search ("\s*parameters\s*", lines[i]):
      skip_line = True
      print "updating Jastrow parameters"
      file_output.write(' parameters\n')
      for k in range(len(lines_jastrow)):
        file_output.write(lines_jastrow[k])
      file_output.write(' end\n')
      in_jastrow_menu=False
      continue

  if csfs_read:
    if re.search ("\s*csf_coef\s*", lines[i]):
      if not re.search ("\s*end\s*", lines[i]):
        skip_line = True
      print "updating CSF parameters"
      file_output.write(' csf_coef\n')
      for k in range(len(lines_csfs)):
        file_output.write(lines_csfs[k])
      file_output.write(' end\n')
      continue

  if orbitals_read:
    if re.search ("\s*coefficients\s*", lines[i]):
      skip_line = True
      print "updating orbital parameters"
      file_output.write(' coefficients\n')
      for k in range(len(lines_orbitals)):
        file_output.write(lines_orbitals[k])
      file_output.write(' end\n')
      continue

  if exponents_read:
    if re.search ("\s*basis_functions\s*", lines[i]):
      skip_line = True
      print "updating exponent parameters"
      file_output.write(' basis_functions\n')
      for k in range(len(lines_exponents)):
        file_output.write(lines_exponents[k])
      file_output.write(' end\n')
      continue
     
  file_output.write(lines[i])

   
# close output file
file_output.close()
#sys.stdout.write('done\n')

