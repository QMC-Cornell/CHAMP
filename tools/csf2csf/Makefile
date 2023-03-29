# This Makefile is to be used with gmake.  On Linux systems, gmake is
# the default make.

system := $(shell uname)
hostname := $(shell uname -n)

# planck_compiler=g77
 planck_compiler=ifort


ifeq ($(system),Linux)

 FC = ifort
 FFLAGS = -zero -save -extend_source -w -r8 -O3 -pad
#LDFLAGS = -nothread

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
#{i|M|K|W} i=Pentium Pro, M-Pentium MMX, K=Pentium3, W=Pentium4,Xeon
# On Pentium4W results in 3% faster execution than not using it and 2% faster execution than using -axi, but uncertainties are equally big
# -O3 seems to be about the same a O2 for ifc (v7).  Have not tested ifort.
# -w turns off warnings
# -fpp2 for preprocessor (not needed for this code)
# Richard Hennig says to use -ip for interprocess (helps a lot with VASP)
# He says VASP did not run correctly with -fast or -xW
# He says to use :ib on mpirun to get infiniband network
# He found -O3WN -ipo -pad runs fastest, but there were several choices within the timing noise.
  FC = ifort
  FFLAGS = -g -CB -d4 -extend_source -w -r8 -fpp2 -O0
# FFLAGS = -extend_source -w -r8 -fpp2 -O2
# FFLAGS = -extend_source -w -r8 -O3N -ipo -pad
# FFLAGS = -extend_source -w -r8 -O3 -pad
# LDFLAGS = -Vaxlib -static
  LDFLAGS = -static
  FC_MPI = mpif77 -fc=ifort
# FFLAGS_MPI = -g -CB -extend_source -w -r8 -O2
# FFLAGS_MPI = -extend_source -w -r8 -O2
  FFLAGS_MPI = -extend_source -w -r8 -O3 -pad
# LDFLAGS_MPI = -Vaxlib -static
  LDFLAGS_MPI = -static -lmpi
endif
endif

# dft cluster at Cornell
ifeq ($(hostname),dft)
  FC=gfortran
  FFLAGS=-finit-local-zero -O3 -ffixed-line-length-132 -Wall
  FFLAGS_DEBUG=$(FFLAGS) -g -fbounds-check
endif

#qmc1 at cornell
ifeq ($(hostname),qmc1)
  FC=ifort
  FFLAGS=-zero -r8 -pad -align all -fpp -ip -O3 -mtune=core2 -march=core2 -xT -extend_source                 
  FFLAGS_DEBUG=$(FFLAGS) -check all -check nooutput_conversion -debug extended -debug-parameters all -g -traceback -error_limit 1 -DDEBUG
endif

#qmc2 at cornell
ifeq ($(hostname),qmc2)
  FC=gfortran
  FFLAGS=-finit-local-zero -O3 -fdefault-real-8 -ffixed-line-length-132 -Wall
  FFLAGS_DEBUG=$(FFLAGS) -g -fbounds-check
endif

ifeq ($(hostname),qmc4)
  FC=gfortran
  FFLAGS=-finit-local-zero -O3 -fdefault-real-8 -ffixed-line-length-132 -Wall
  FFLAGS_DEBUG=$(FFLAGS) -g -fbounds-check
endif

ifeq ($(hostname),nanolab.cnf.cornell.edu)
  FC = ifort
# FFLAGS = -g -CB -d4 -zero -save -extend_source -w -r8 -fpp2 -O0
# FFLAGS = -zero -save -extend_source -w -r8 -fpp2 -O2
# FFLAGS = -zero -save -extend_source -w -r8 -O3 -ipo -pad
  FFLAGS = -zero -save -extend_source -w -r8 -O3 -pad
# LDFLAGS = -Vaxlib -static
# LDFLAGS = -static -nothread
  LDFLAGS = -nothread
  FC_MPI = /usr/lam_intel/bin/mpif77
# FFLAGS_MPI = -g -CB -zero -save -extend_source -w -r8 -O0
# FFLAGS_MPI = -g -CB -zero -save -extend_source -w -r8 -O2
# FFLAGS_MPI = -zero -save -extend_source -w -r8 -O2
  FFLAGS_MPI = -zero -save -extend_source -w -r8 -O3 -pad
# LDFLAGS_MPI = -Vaxlib -static
# LDFLAGS_MPI = -static -nothread -lmpi
  LDFLAGS_MPI = -nothread -lmpi -lpthread
endif

endif

# det_to_csf and det2csf are outdated
all: ie csf2csf

ie: ie.o
	$(FC) $(FFLAGS) $(LDFLAGS) ie.o -o $@

det_to_csf: det_to_csf.o
	$(FC) $(FFLAGS) $(LDFLAGS) det_to_csf.o -o $@

det2csf: det2csf.o
	$(FC) $(FFLAGS) $(LDFLAGS) det2csf.o -o $@

csf2csf: csf2csf.o
	$(FC) $(FFLAGS) $(LDFLAGS) csf2csf.o -o $@
