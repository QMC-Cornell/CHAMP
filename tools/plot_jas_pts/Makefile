QMCHOME = ..
LIBHOME = $(HOME)
QUENCH  = -L$(LIBHOME)/SimulatedAnnealing/quench_anneal/lib -lquench -lquench_seq
# QUENCHMPI=-L$(LIBHOME)/SimulatedAnnealing/quench_anneal/lib -lquench_mpi
CYRUSLIB= -L$(LIBHOME)/lib -lcyrus
BLAS    = -L$(LIBHOME)/lib2/blas -lblas
LINPACK = -L$(LIBHOME)/lib2/linpack -llinpack

#FC     = pgf77
#FFLAGS = -fast -r8 -Mextend -Msave -Mdalign -Mbounds
#FFLAGS = -fast -r8 -Mextend -Msave

 FC     = g77
#FFLAGS = -O2 -malign-double -ffixed-line-length-none -fno-automatic -Wall
 FFLAGS = -O3 -malign-double -ffixed-line-length-none -fno-automatic

#FC     = xlf
#FFLAGS = -gC -qextchk -qsave -qspnans -qdpc -bloadmap:map
#FFLAGS =  -O2 -qextchk -qsave -qspnans -qdpc -bloadmap:map

plot_jas_pts: plot_jas_pts.o
	$(FC) $(FFLAGS) plot_jas_pts.o -o $@
