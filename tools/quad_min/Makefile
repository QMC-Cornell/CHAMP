# This Makefile is to be used with gmake.  On Linux systems, gmake is
# the default make.

system := $(shell uname)
hostname := $(shell uname -n)

planck_compiler=g77
# planck_compiler=ifort

LAPACK = -llapack

ifeq ($(system),Linux)

ifeq ($(hostname),planck)
# AMD Athlon cluster of Arias group
ifeq ($(planck_compiler),g77)
# to link g77 mpi libraries, make sure /usr/bin is ahead of /usr/local/bin
# FC = g77-3.2
  FC = g77-3.3
# FC = g77
  FFLAGS = -g -C -malign-double -ffixed-line-length-none -fno-automatic -Wall -fbounds-check
# FFLAGS = -O2 -malign-double -ffixed-line-length-none -fno-automatic -Wall -fbounds-check
# FFLAGS = -O3 -malign-double -ffixed-line-length-none -fno-automatic
# FFLAGS = -pg -O3 -malign-double -ffixed-line-length-none -fno-automatic
  LDFLAGS = -static
  FC_MPI = mpif77 -fc=g77
# FFLAGS_MPI = -O2 -malign-double -ffixed-line-length-none -fno-automatic -Wall -fbounds-check
  FFLAGS_MPI = -O3 -malign-double -ffixed-line-length-none -fno-automatic
  LDFLAGS_MPI = -static
endif

ifeq ($(planck_compiler),ifort)
# to link ifort mpi libraries, make sure /usr/local/bin is ahead of /usr/bin
# In the mpif77 line, the fc=ifc tells it to use ifc rather than the default g77
# fpp2 is the preprocessor
# -ax{i|M|K|W} i=Pentium Pro, M-Pentium MMX, K=Pentium3, W=Pentium4,Xeon
# On Pentium4 -axW results in 3% faster execution than not using it and 2% faster execution than using -axi, but uncertainties are equally big
# -O3 seems to be about the same a O2 for ifc (v7).  Have not tested ifort.
# -w turns off warnings
# -fpp2 for preprocessor (not needed for this code)
# Richard Hennig says to use -ip for interprocess (helps a lot with VASP)
# He says VASP did not run correctly with -fast or -xW
# He says to use :ib on mpirun to get infiniband network
# He found -O3 -axWN -ipo -pad runs fastest, but there were several choices within the timing noise.
  FC = ifort
# FFLAGS = -g -C -d4 -extend_source -w -r8 -fpp2 -O0 -axi
# FFLAGS = -extend_source -w -r8 -fpp2 -O2 -axi
# FFLAGS = -extend_source -w -r8 -O3 -axiN -ipo -pad
  FFLAGS = -extend_source -w -r8 -O3 -axi -pad
# LDFLAGS = -Vaxlib -static
  LDFLAGS = -static
  FC_MPI = mpif77 -fc=ifort
# FFLAGS_MPI = -g -C -extend_source -w -r8 -O2 -axi
# FFLAGS_MPI = -extend_source -w -r8 -O2 -axi
  FFLAGS_MPI = -extend_source -w -r8 -O3 -axi -pad
# LDFLAGS_MPI = -Vaxlib -static
  LDFLAGS_MPI = -static -lmpi
endif
endif

ifeq ($(hostname),nanolab.cnf.cornell.edu)
  FC = ifort
  LAPACK  = -L/opt/intel/mkl72/lib/32 -lmkl_lapack -lmkl_ia32 -lguide -lpthread
# FFLAGS = -g -CB -d4 -zero -save -extend_source -w -r8 -fpp2 -O0 -ax
# FFLAGS = -zero -save -extend_source -w -r8 -fpp2 -O2 -ax
# FFLAGS = -zero -save -extend_source -w -r8 -O3 -axN -ipo -pad
  FFLAGS = -zero -save -extend_source -w -r8 -O3 -pad
  LDFLAGS = -Vaxlib -static -nothread
# LDFLAGS = -static -nothread
  FC_MPI = /usr/local/bin/mpif77 -fc=ifort
# FFLAGS_MPI = -g -CB -zero -save -extend_source -w -r8 -O0 -ax
# FFLAGS_MPI = -g -CB -zero -save -extend_source -w -r8 -O2 -ax
# FFLAGS_MPI = -zero -save -extend_source -w -r8 -O2 -ax
  FFLAGS_MPI = -zero -save -extend_source -w -r8 -O3 -ax -pad
# LDFLAGS_MPI = -Vaxlib -static
  LDFLAGS_MPI = -static -nothread -lmpi
endif

endif

quad_min: quad_min.o
	$(FC) $(FFLAG) $(LDFLAGS) quad_min.o -o $@
 
newton_chol_cyrus: newton_chol_cyrus.o
	$(FC) $(FFLAGS) newton_chol_cyrus.o $(LAPACK) -o $@
 
linear: linear.o
	$(FC) $(FFLAGS) linear.o $(LAPACK) -o $@
 
nigh: nigh.o
	$(FC) $(FFLAGS) nigh.o $(LAPACK) -o $@
