# This Makefile is to be used with gmake.  On Linux systems, gmake is
# the default make.
# This Makefile should not be used on its own but should instead be called
# from the Makefile in lib2 (to set FC and FFLAGS).

FC = gfortran
# FFLAGS=-finit-local-zero -O3 -ffixed-line-length-132 -x f77-cpp-input -Wall -fbounds-check
FFLAGS=-finit-local-zero -O3 -ffixed-line-length-132 -x f77-cpp-input -Wall

OBJS = dasum.o daxpy.o dcopy.o ddot.o dnrm2.o drot.o drotg.o dscal.o dswap.o idamax.o \
       dgemm.o dgemv.o dger.o dsymv.o dsyr2.o dsyr2k.o dtrmm.o dtrmv.o \
       dsymm.o dsyr.o dsyrk.o dtrsm.o dtrsv.o

libblas.a: $(OBJS)
	$(AR) $(ARFLAGS) $@ $(OBJS)
	ranlib $@

clean:
	-rm -f *.o

clean_all:
	-rm -f *.o *.a
