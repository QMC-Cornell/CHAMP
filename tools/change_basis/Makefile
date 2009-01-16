# This Makefile is to be used with gmake.  On Linux systems, gmake is
# the default make.

  FC = ifort
  FFLAGS = -zero -save -extend_source -w -r8 -O3 -pad
# LDFLAGS = -static

change_basis: change_basis.o ../../qmc/vmc/matinv.o
	$(FC) $(FFLAGS) $(LDFLAGS) change_basis.o ../../qmc/vmc/matinv.o -o $@
