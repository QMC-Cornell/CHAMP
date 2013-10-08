# This Makefile is to be used with gmake.  On Linux systems, gmake is the default make.

LIBHOME=/home/cyrus/lib
LIBS=-L$(LIBHOME) -lcyrus

#FC=ifort
#FFLAGS=-zero -r8 -pad -align all -fpp -ip -O3 -mtune=core2 -march=core2 -xT -extend_source
#F90FLAGS= -zero -r8 -pad -align all -O3 -fpp -ip
#F90FLAGS= -zero -r8 -pad -align all -O3 -fpp -pg
#F90FLAGS=-zero -r8 -pad -align all -O0 -fpp -ip -check all -check nooutput_conversion -debug extended -debug-parameters all -g -traceback -error_limit 1 -DDEBUG

 FC=gfortran
#FFLAGS=-finit-local-zero -O0  -ffixed-line-length-132 -Wall -g -fbounds-check
#F90FLAGS = -finit-local-zero -O0  -ffree-line-length-none -x f95-cpp-input -Wall -g -fbounds-check
FFLAGS=-g -finit-local-zero -O3 -ffixed-line-length-132 -Wall -pg -fbacktrace -fbounds-check
F90FLAGS =-g -finit-local-zero -O3 -ffree-line-length-none -x f95-cpp-input -Wall -pg -fbacktrace -fbounds-check
#FFLAGS=-finit-local-zero -O3 -ffixed-line-length-132 -Wall
#F90FLAGS = -finit-local-zero -O3 -ffree-line-length-none -x f95-cpp-input -Wall #-DNUM_ORBITALS_GT_127
# LDFLAGS = -L/usr/local/tools/mpip/lib -lmpiP -lm#-static
 MFC = mpif90
 MPI = -cpp #-DMPI 
#LDFLAGS = -static

.SUFFIXES:
.SUFFIXES: .f90 .f95 .o .f .c

default: symmetrize_den

.c.o:
	$(CC) $(CC_FLAGS) -o $@ -c $<

.f.o:
	$(FC) $(FFLAGS) -o $@ -c $<

.f90.o:
	$(FC) $(F90FLAGS) $(MPI) -o $@ -c $<

symmetrize_den.o: types.o

symmetrize_den: symmetrize_den.o types.o
	$(FC) $(FFLAGS) $(LDFLAGS) symmetrize_den.o types.o $(LIBS) -o $@

clean:
	rm -f $(OBJS) *.mod
