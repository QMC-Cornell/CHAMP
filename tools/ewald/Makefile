#SHELL  = /bin/csh

 FC     = g77
 FFLAGS = -g -C -O2 -malign-double -ffixed-line-length-none -fno-automatic -Wall -fbounds-check
#FFLAGS = -O3 -malign-double -ffixed-line-length-none -fno-automatic

LIBHOME = $(HOME)
# LIBHOME = $(QMCHOME)
# # LIBHOME = /data/bega/filippi
# QUENCH  = -L$(LIBHOME)/SimulatedAnnealing/quench_anneal/lib -lquench -lquench_seq
# QUENCHMPI=-L$(LIBHOME)/SimulatedAnnealing/quench_anneal/lib -lquench_mpi
CYRUSLIB= -L$(LIBHOME)/lib -lcyrus
BLAS    = -L$(LIBHOME)/lib2/blas -lblas
LINPACK = -L$(LIBHOME)/lib2/linpack -llinpack

V = $(HOME)/qmc/vmc

# LIB = $(HOME)/lib/lib.a $(HOME)/lib2/linpack/liblinpack.a $(HOME)/lib2/blas/libblas.a

INC    = -I$V
FFLAGS := $(FFLAGS) $(INC)

SRC_EWALD = test_ewald.f $V/ewald.f $V/distances.f $V/find_image.f $V/matinv.f \
            $V/readps_tm.f $V/dgelg.f $V/spline2.f $V/rotqua.f $V/rannyu.f

OBJ_EWALD = test_ewald.o $V/ewald.o $V/distances.o $V/find_image.o $V/matinv.o \
# OBJ_EWALD = test_ewald.o ewald.o distances.o $V/find_image.o $V/matinv.o \
            $V/readps_tm.o $V/dgelg.o $V/spline2.o $V/rotqua.o $V/rannyu.o

OBJ_EWALD_new = test_ewald.o ewald.o $V/find_image.o $V/matinv.o \
            $V/readps_tm.o $V/dgelg.o $V/spline2.o $V/rotqua.o $V/rannyu.o

# OBJ_PSP = test_psp.o ewald.o $V/find_image.o $V/matinv.o \

OBJ_PSP = test_psp.o $V/find_image.o $V/matinv.o \
            readps_tm.o $V/dgelg.o spline2.o $V/rotqua.o $V/rannyu.o

OBJ_IMG = test_find_image.o find_image.o $V/matinv.o

ewald.debug: $(SRC_EWALD)
	$(HOME)/ftnchek-3.1.2*/ftnchek -columns=130 -pretty=no-long-line -array=no-dimensions \
          $(SRC_EWALD) > oewald.ftnchek

ewald: $(OBJ_EWALD)
	$(FC) $(FFLAGS) $(OBJ_EWALD) $(CYRUSLIB) $(LINPACK) $(BLAS) -o $@

ewald_new: $(OBJ_EWALD_new)
	$(FC) $(FFLAGS) $(OBJ_EWALD_new) $(CYRUSLIB) $(LINPACK) $(BLAS) -o $@

psp: $(OBJ_PSP)
	$(FC) $(FFLAGS) $(OBJ_PSP) $(CYRUSLIB) $(LINPACK) $(BLAS) -o $@

find_image: $(OBJ_IMG)
	$(HOME)/ftnchek-3.1.2*/ftnchek -columns=130 -pretty=no-long-line -array=no-dimensions test_find_image.f find_image.f $V/matinv.f
	$(FC) $(FFLAGS) $(OBJ_IMG) -o $@
