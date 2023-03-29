# This directory contains the subset of the general purpose math routines written by
# Cyrus Umrigar that are used in CHAMP.  A more complete set of routines is:
# OBJS = my_second.o ranf.o shellnn.o intpol.o intpol_monot.o grater.o \
# hartre.o deriv.o deriv2.o divergence.o lapla3.o laplac.o fderiv.o simson.o bode.o \
# hartre_grid2.o deriv_grid2.o deriv2_grid2.o divergence_grid2.o laplac_grid2.o fderiv_grid2.o \
# bode_grid2.o

FC = gfortran
FFLAGS=-finit-local-zero -O3 -ffixed-line-length-132 -x f77-cpp-input -Wall -fbounds-check

OBJS = my_second.o intpol.o simson.o shell.o

libcyrus.a: $(OBJS)
	ar r $@ $(OBJS)
	ranlib $@

clean:
	-rm -f *.o

clean_all:
	-rm -f *.o *.a
