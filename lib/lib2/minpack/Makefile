# This Makefile is to be used with gmake.  On Linux systems, gmake is
# the default make.
# This Makefile should not be used on its own but should instead be called
# from the Makefile in lib2 (to set FC and FFLAGS).

ifndef FFLAGS
 FC     = gfortran
 FFLAGS = -finit-local-zero -O3 -ffixed-line-length-132 -Wall
endif

OBJS = chkder.o  dpmpar.o  fdjac2.o  hybrd1.o   lmder.o   lmdipt.o  objfcn.o  r1mpyq.o  ssqjac.o \
covar.o   enorm.o   grdfcn.o  hybrj.o    lmder1.o  lmpar.o   qform.o   r1updt.o  vecfcn.o \
dmchar.o  errjac.o  hesfcn.o  hybrj1.o   lmdif.o   lmstr.o   qrfac.o   rwupdt.o  vecjac.o \
dogleg.o  fdjac1.o  hybrd.o   lhesfcn.o  lmdif1.o  lmstr1.o  qrsolv.o  ssqfcn.o

libminpack.a: $(OBJS)
	ar r $@ $(OBJS)
	ranlib $@

clean:
	rm -f *.o

clean_all:
	rm -f *.o *.a
