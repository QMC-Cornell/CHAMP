#!/usr/bin/python
#takes basis functions from GAMESS input and puts them on a radial grid.
#Work with both Gaussian and Gauss-Slater basis functions.
#Written by Frank Petruzielo 10/18/11. Last modified 10/20/11
import getopt
import sys
import re
import numpy
import math

def usage():
    print "./basis_grid.py gamess_input_file"

################################################################################################################################## 
class BasisGrid(object):
    '''Orbitals from a gamess calculation

    Input (stored as attributes):
    input_file:                      name of input file containing gamess calculation

    Other attributes:
    basis_functions:                 Store all basis function information in array. This info is set with
                                     collect_basis_funcions
    n_basis:                         Dictionary storing number of basis functions for each atom type
    atom_names:                      List of unique atoms

    Methods:
    collect_basis_functions:         Get all basis function information: atom type, function type,
                                     angular momentum, and exponents and coefficients if contraction
                                     or exponent if GS primitive
    create_grid:                     Put basis_basis functions on a radial (exponential) grid
    '''
    
    def __init__(self, input_file):
        self.input_file = input_file

    def create_grid(self, n_grid = 1600, r_max = 7.000000, a_density = 1.003000, grid_type = 3, icusp = 0):
        '''
        Put basis_basis functions on a radial (exponential) grid
        Input
        n_grid:                      Number of grid points                  
        r_max:                       Maximum radial value
        a_density:                   The larger the parameter is, the more dense the grid is near the origin
        grid_type:                   Not used. Can be implemented simply
        icusp:                       Not used.
        '''
        #create radial exponential grid
        r_scaled = r_max / (a_density**(n_grid-1) - 1)
        grid = []
        for k in range(n_grid):
            grid.append(r_scaled * (a_density ** k - 1))

        #evaluate basis functions on the grid. Store them in a dictionary indexed by atom name.
        #initialize
        values = {}
        for i in self.atom_names:
            values[i] = []        
            
        #add values for each grid point
        for r in grid:
            for i in self.atom_names:
                values[i].append([])        
            for i in self.basis_functions:
                values[i.atom_name][-1].append(i.get_value(r))
                
        #write to files
        for i in range(len(self.atom_names)):
            basis_file = re.sub("\.inp", ".basis." + self.atom_names[i].upper() + "." + str(i+1) + "_new", self.input_file)
            f = open(basis_file, 'w')
            f.write("%d %d %d %f %f %d   nrbas,igrid,nr,exp_h_bas,r0_bas,icusp. Note that the exact form of GS functions was used.\n" % (self.n_basis[self.atom_names[i]], grid_type, n_grid, a_density, r_max, icusp))
            format_string = "  %.12E" + len(values[self.atom_names[i]][0]) * " %.12E"
            for k in range(len(grid)):
                temp_values = [grid[k]] + values[self.atom_names[i]][k]
                f.write(format_string % tuple(temp_values) + "\n")
            f.close()

                ################################################################################################################################## 
    def collect_basis_functions(self):
        '''
        Get all basis function information: atom type, function type,
        angular momentum, and exponents and coefficients if contraction or exponent if GS primitive
        '''
        f = open(self.input_file, 'r')
        in_basis = False
        in_atom = False
        self.basis_functions = []

        for line in f:

            #looking for start of basis section in input file
            match = re.search("^\s*\$DATA", line, re.IGNORECASE)
            if match:        
                in_basis = True
                continue

            if in_basis:
                #looking for end of basis section in input file
                match = re.search("^\s*\$END", line, re.IGNORECASE)
                if match:
                    break
                        
                #find the line with the atom name and store this
                match = re.search("^\s*(\w+).*\.?.*\..*\..*\.", line, re.IGNORECASE)
                if match:
                    atom_name  = match.group(1).lower()
                    in_atom = True
                    continue

                if in_atom:
                    #search for the angular momentum and number of expansion coefficients if applicable
                    match = re.search("^\s*([a-z])\s+(\d+)\s+", line, re.IGNORECASE)
                    if match:
                        ang_mom = match.group(1).lower()
                        n_expand = int(match.group(2))
                        coef = []
                        exp = []
                        continue
                    #collect expansion
                    match = re.search("^\s*\d+\s+(\d+\.\d+)\s+(-?\d+\.\d+)\s*$", line, re.IGNORECASE)
                    if match:
                        exp.append(float(match.group(1)))
                        coef.append(float(match.group(2)))
                        if len(exp) == n_expand:
                            func_type = "contraction"
                            self.basis_functions.append(BasisFunction(atom_name, func_type, ang_mom, exp, coef))
                        continue
                    #search for the angular momentum and exponent for GS function
                    match = re.search("^\s*GS-([a-z])\s+\d+\s+(\d+\.\d+)\s*$", line, re.IGNORECASE)
                    if match:
                        ang_mom = match.group(1).lower()
                        exp = float(match.group(2))
                        func_type = "gs"
                        self.basis_functions.append(BasisFunction(atom_name, func_type, ang_mom, exp))
                        continue

        f.close()

        #determine the unique atoms in the calculation and the number of basis functions for each atom
        temp_atom_names = [self.basis_functions[i].atom_name for i in range(len(self.basis_functions))]
        unique_atom_names = []
        for element in temp_atom_names:
            if element not in unique_atom_names:
                unique_atom_names.append(element)

        self.atom_names = unique_atom_names
        self.n_basis = {}
        for i in unique_atom_names:
            self.n_basis[i] = temp_atom_names.count(i)
                      
################################################################################################################################## 
class BasisFunction(object):
    '''Basis function for a quantum chemistry calculation

    Input (stored as attributes):
    atom_name:                    name of the atom the basis function is for
    func_type:                    contraction or Gauss-Slater
    ang_mom:                      angular momentum of the function

    Other attributes:
    exp:                          exponents in contraction
    coef:                         coefficients in contraction

    Methods:
    get_value:                    Get value of radial part of basis function at r

    '''    
    def __init__(self, atom_name, func_type, ang_mom, exp, coef = 1.):
        self.atom_name = atom_name
        self.func_type = func_type
        self.ang_mom = ang_mom
        self.exp = exp
        self.coef = coef

        accepted_types = ["contraction", "gs"]
        accepted_angular_mom = ["s", "p", "d", "f"]

        if self.func_type not in accepted_types:
            print "Basis function is not of an acceptable type"
            sys.exit(2)

        if self.ang_mom not in accepted_angular_mom:
            print "Basis function is not of an acceptable angular momentum"
            sys.exit(2)
        
    def get_value(self, r):
        '''
        Get value of radial part of basis function at r
        ''' 
        ang_mom_to_n = {"s": 1, "p": 2, "d": 3, "f": 4}
        n = ang_mom_to_n[self.ang_mom]
        if self.func_type == "gs":
            # gauss-slater function
            normalization_dict = {1: 1.12646742161049, 2: 0.5766099503612371, 3: 0.1965811411216233, 4: 0.05027565586963441}
            normalization = self.exp ** (n + 0.5) * normalization_dict[n]
            value = normalization * r ** (n-1) * numpy.exp( - (r * self.exp)**2 / (1. +  r * self.exp))
        else:
            # contraction of gaussians
            value = 0.
            for i in range(len(self.exp)):
                normalization = (2.0 * (2.0 * self.exp[i]) ** (n + 0.5) / math.gamma(n + 0.5)) ** 0.5
                value += self.coef[i] * normalization * r **(n-1) * numpy.exp(-r**2 * self.exp[i]) 
        return value

#############################################################################################
def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "ht:", ["help"])
    except getopt.GetoptError, err:
        print >> sys.stderr, str(err) 
        usage()
        sys.exit(2)        
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        else:
            assert False, "unhandled option"
    if args:
        input_file = args[0].strip()
    else:
        usage()
        sys.exit(2)

    basis_grid = BasisGrid(input_file)
    basis_grid.collect_basis_functions()
    basis_grid.create_grid()

if __name__ == '__main__':
    main()
