FC = gfortran
F90FLAGS = -g -Wall -fbounds-check -fbacktrace -ffree-line-length-none

.SUFFIXES:
.SUFFIXES: .f90 .f95 .o .f .c

.f90.o:
	$(FC) $(F90FLAGS) $(MPI) -o $@ -c $<

select_csf: select_csf.o
	$(FC) $(F90FLAGS) $(LDFLAGS) select_csf.o -o $@
