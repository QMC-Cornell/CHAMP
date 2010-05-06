#SHELL  = /bin/csh

# FC     = g77
# FFLAGS = -g -C -O2 -malign-double -ffixed-line-length-none -fno-automatic -Wall -fbounds-check
##FFLAGS = -O3 -malign-double -ffixed-line-length-none -fno-automatic

 FC     = gfortran
 FFLAGS = -finit-local-zero -O3 -fdefault-real-8 -ffixed-line-length-132

LIBHOME = $(HOME)
# LIBHOME = $(QMCHOME)
# # LIBHOME = /data/bega/filippi
# QUENCH  = -L$(LIBHOME)/SimulatedAnnealing/quench_anneal/lib -lquench -lquench_seq
# QUENCHMPI=-L$(LIBHOME)/SimulatedAnnealing/quench_anneal/lib -lquench_mpi
CYRUSLIB= -L$(LIBHOME)/lib -lcyrus
BLAS    = -L$(LIBHOME)/lib2/blas -lblas
LINPACK = -L$(LIBHOME)/lib2/linpack -llinpack

psp_needs_tm: psp_needs_tm.o
	$(FC) $(FFLAGS) psp_needs_tm.o -o $@
