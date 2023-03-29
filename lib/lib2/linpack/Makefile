# This Makefile is to be used with gmake.  On Linux systems, gmake is
# the default make.
# This Makefile should not be used on its own but should instead be called
# from the Makefile in lib2 (to set FC and FFLAGS).

FC = gfortran
# FFLAGS=-finit-local-zero -O3 -ffixed-line-length-132 -x f77-cpp-input -Wall -fbounds-check
FFLAGS=-finit-local-zero -O3 -ffixed-line-length-132 -x f77-cpp-input -Wall

OBJS = dcabs1.o dchdc.o dchdd.o dchex.o dchud.o dgbco.o dgbdi.o dgbfa.o dgbsl.o dgeco.o \
dgedi.o dgefa.o dgesl.o dgtsl.o dmach.o dpbco.o dpbdi.o dpbfa.o dpbsl.o dpoco.o \
dpodi.o dpofa.o dposl.o dppco.o dppdi.o dppfa.o dppsl.o dptsl.o dqrdc.o dqrsl.o \
dsico.o dsidi.o dsifa.o dsisl.o dspco.o dspdi.o dspfa.o dspsl.o dsvdc.o dtrco.o \
dtrdi.o dtrsl.o
# zchdc.o zchdd.o zchex.o zchud.o zgbco.o zgbdi.o zgbfa.o zgbsl.o \
# zgeco.o zgedi.o zgefa.o zgesl.o zgtsl.o zhico.o zhidi.o zhifa.o zhisl.o zhpco.o \
# zhpdi.o zhpfa.o zhpsl.o zmach.o zpbco.o zpbdi.o zpbfa.o zpbsl.o zpoco.o zpodi.o \
# zpofa.o zposl.o zppco.o zppdi.o zppfa.o zppsl.o zptsl.o zqrdc.o zqrsl.o zsico.o \
# zsidi.o zsifa.o zsisl.o zspco.o zspdi.o zspfa.o zspsl.o zsvdc.o ztrco.o ztrdi.o \
# ztrsl.o
 
liblinpack.a: $(OBJS)
	ar r $@ $(OBJS)
	ranlib $@

clean:
	-rm -f *.o

clean_all:
	-rm -f *.o *.a
