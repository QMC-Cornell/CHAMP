# This Makefile is to be used with gmake.  On Linux systems, gmake is
# the default make.
# This Makefile should not be used on its own but should instead be called
# from the Makefile in lib2 (to set FC and FFLAGS).

OBJS = dasum.o daxpy.o dcopy.o ddot.o dnrm2.o drot.o drotg.o dscal.o dswap.o idamax.o \
       dgemm.o dgemv.o dger.o dsymv.o dsyr2.o dsyr2k.o dtrmm.o dtrmv.o \
       dsymm.o dsyr.o dsyrk.o dtrsm.o dtrsv.o

libblas.a: $(OBJS)
	$(AR) $(ARFLAGS) $@ $(OBJS)
	ranlib $@
