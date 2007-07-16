# This Makefile is to be used with gmake.  On Linux systems, gmake is
# the default make.

system := $(shell uname)
hostname := $(shell uname -n)
planck_compiler=g77
QMCHOME = $(shell pwd)
LIBHOME = $(HOME)
# LIBHOME = $(QMCHOME)
QUENCH  = -L$(LIBHOME)/SimulatedAnnealing/quench_anneal/lib -lquench_g77 -lquench_seq_g77
# QUENCH  = -L$(LIBHOME)/SimulatedAnnealing/quench_anneal/lib -lquench -lquench_seq -llapack -lblas
CYRUSLIB= -L$(LIBHOME)/lib -lcyrus_g77
BLAS    = -L$(LIBHOME)/lib2/blas -lblas_g77
LINPACK = -L$(LIBHOME)/lib2/linpack -llinpack_g77
# LAPACK  = -L$(LIBHOME)/SimulatedAnnealing/quench_anneal.new2/lib -llapack_g77
PSPLINE = -L$(LIBHOME)/lib2/pspline/pspline -lpspline_g77

PERIODIC = T

export FC FFLAGS LDFLAGS FC_MPI FFLAGS_MPI LDFLAGS_MPI LAMHF77 \
       PERIODIC QMCHOME LIBHOME CYRUSLIB QUENCH QUENCHMPI \
       BLAS LINPACK LAPACK PSPLINE

ifeq ($(findstring IRIX,$(system)),IRIX)
  FC = f77
  FFLAGS = -O3 -n32 -mips4 -DIRIX
# FFLAGS = -g -C -n32 -mips4 -DIRIX
  FC_MPI = mpif77
  FFLAGS_MPI = -O3 -n32 -mips4 -DIRIX -DMPI -DNON_STANDARD
# FFLAGS_MPI = -g -C -n32 -mips4 -DIRIX -DMPI -DNON_STANDARD
endif

# This was originally for the Cornell SP2, later modified for the NERSC SP3, compile node Seaborg.
# Wolfgang needed to replace flush with flush_.  -bmaxdata:0x7FFFFFFF to get 2GB of RAM.
ifeq ($(system),AIX)
  FC = xlf
  FFLAGS = -O2 -bmaxdata:0x70000000 -qnosave -qsigtrap -WF,-DAIX -qfixed=132
  FC_MPI = mpxlf
  FFLAGS_MPI = -O2 -qfixed=132 -bmaxdata:0x70000000 -qnosave -qsigtrap -WF,-DAIX,-DMPI,-DNON_STANDARD
endif

ifeq ($(system),Linux)

ifeq ($(hostname),planck)
# AMD Athlon cluster of Arias group
ifeq ($(planck_compiler),g77)
# to link g77 mpi libraries, make sure /usr/bin is ahead of /usr/local/bin
# FC = g77-3.2
  FC = g77-3.0
# FC = g77
# FFLAGS = -g -C -malign-double -ffixed-line-length-none -fno-automatic -Wall -fbounds-check
# FFLAGS = -O2 -malign-double -ffixed-line-length-none -fno-automatic -Wall -fbounds-check
  FFLAGS = -O3 -malign-double -ffixed-line-length-none -fno-automatic
# FFLAGS = -pg -O3 -malign-double -ffixed-line-length-none -fno-automatic
# LDFLAGS = -static
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

ifeq ($(hostname),charming)
# If you have the Portland Group pgf compiler, it is probably better.
  FC    = pgf77
# FFLAGS = -fast -r8 -Mextend -Mdalign -Msave -Mbounds
  FFLAGS = -fast -r8 -Mextend -Mdalign -Msave
  FC_MPI = mpif77.pgf
  FFLAGS_MPI = -fast -r8 -DCLOBBER -DMPI
# FFLAGS_MPI = -O -DMPI -DAIX -YEXT_NAMES=LCS -I/usr/local/lam/include
endif

ifeq ($(hostname),memling.lorentz.leidenuniv.nl)
# -Msecond_underscore and -I/usr/include required for LAM-MPI
# -Msecond_underscore puts 2 underscores rather than one
  FC    = pgf77
# FFLAGS = -O0 -g -r8 -Mextend -Mdalign -Msave -Msecond_underscore
  FFLAGS = -fast -r8 -Mextend -Mdalign -Msave -Msecond_underscore
# FFLAGS = -O0 -g -r8 -Mextend -Mdalign -Msave
# FC_MPI = /opt/abinitio/bin/mpif77
  FC_MPI =  hf77
# FFLAGS_MPI = -g -r8 -Mextend -Mdalign -Msave -I/usr/include -Msecond_underscore
  FFLAGS_MPI = -fast -r8 -Mextend -Mdalign -Msave -I/usr/include -Msecond_underscore
# FFLAGS_MPI = -fast -r8 -Mextend -Mdalign -Msave -I/usr/include
endif

# ifeq ($(hostname),oscbw01)
#ifeq ($(hostname),amd-login1)
## AMD Cluster at OSC using Portland Group Compiler
#  FC = pgf77
## FFLAGS = -fast -r8 -Mextend -Mdalign -Msave -Mbounds
#  FFLAGS = -fast -r8 -Mextend -Mdalign -Msave
#  FC_MPI = mpif77
## FFLAGS_MPI = -fast -r8 -Mextend -Mdalign -Msave -Mbounds
#  FFLAGS_MPI = -fast -r8 -Mextend -Mdalign -Msave
#endif

ifeq ($(hostname),amd-login1)
# AMD Cluster (128 1.4 GHz dual proc 2GB) at OSC using Intel ifc Compiler
# In the mpif77 line, the fc=ifc tells it to use ifc rather than the default pgf
# fpp2 is the preprocessor
# -ax{i|M|K|W} i=Pentium Pro, M-Pentium MMX, K=Pentium3, W=Pentium4,Xeon
# On AMD, -axi results in 2% faster execution than not using it or using -axW, but uncertainties are equally big
# -O3 seems to about the same a O2.
# -Vaxlib is needed for flush to work under ifc.  Not needed for ifort.
  FC = ifc
# FFLAGS = -g -C -d4 -FI -extend_source -O0 -w -r8 -fpp2 -axi
# FFLAGS = -FI -extend_source -O0 -w -r8 -fpp2 -axi
  FFLAGS = -FI -extend_source -O2 -w -r8 -fpp2 -axi -ipo -pad
  LDFLAGS = -Vaxlib -static
  FC_MPI = mpif77 -fc=ifc
# FFLAGS_MPI = -g -C -FI -extend_source -O2 -w -r8 -axi
  FFLAGS_MPI = -FI -extend_source -O2 -w -r8 -axi -ipo -pad
  LDFLAGS_MPI = -Vaxlib -static
endif

ifeq ($(hostname),piv-login1)
# Pentium4 cluster (256 2.4 GHz dual proc 4GB) at OSC using Intel ifort Compiler
# In the mpif77 line, the fc=ifc tells it to use ifc rather than the default pgf
# fpp2 is the preprocessor
# -ax{i|M|K|W} i=Pentium Pro, M-Pentium MMX, K=Pentium3, W=Pentium4,Xeon
# On Pentium4 -axW results in 3% faster execution than not using it and 2% faster execution than using -axi, but uncertainties are equally big
# -O3 seems to be about the same a O2 for ifc (v7).  Have not tested ifort.
# -w turns off warnings
# -fpp2 for preprocessor (not needed for this code)
# Richard Hennig says to use -ipo for interprocess (helps a lot with VASP) (but adds to compilation time)
# He says VASP did not run correctly with -fast or -xW
# He says to use :ib on mpirun to get infiniband network
# He found -O3 -axWN -ipo -pad runs fastest, but there were several choices within the timing noise.
  FC = ifort
# FFLAGS = -g -C -d4 -extend_source -w -r8 -fpp2 -O0 -axW
# FFLAGS = -extend_source -w -r8 -fpp2 -O2 -axW
  FFLAGS = -extend_source -w -r8 -O3 -axWN -ip -pad
# LDFLAGS = -Vaxlib -static
  LDFLAGS = -static
  FC_MPI = mpif77
# FFLAGS_MPI = -g -C -extend_source -w -r8 -O2 -axW
# FFLAGS_MPI = -extend_source -w -r8 -O2 -axW
  FFLAGS_MPI = -extend_source -w -r8 -O3 -axWN -ip -pad
# LDFLAGS_MPI = -Vaxlib -static
  LDFLAGS_MPI = -static
endif

ifeq ($(hostname),ia64)
# Intel Itanium cluster at OSC using MPI and Intel Compiler efc
# -Vaxlib is needed while linking for flush(6) to work.
# Richard says the built-in error fn., erf, in Intel Fortran compiler 7.1 is broken
# so use the mkl library 6.0 error functions, derf, derfc by linking to mkl_vml library.
# -FI -extend_source or -FI -4L132 do fixed format with 132 columns
# LIBHOME = $(QMCHOME)
  FC = efc
  FFLAGS = -FI -extend_source -O2 -w -r8 -Dderf=INTELERF -Dderfc=INTELERFC -fpp2
  LDFLAGS = -Vaxlib -L/usr/local/mkl-6.0-beta/lib -lmkl_vml_itp
  FC_MPI = mpif77
  FFLAGS_MPI = -FI -extend_source -O2 -w -r8
  LDFLAGS_MPI = -Vaxlib -L/usr/local/mkl-6.0-beta/lib -lmkl_vml_itp
# BLAS = $(MKL) -lderf
  BLAS = $(MKL) -L$(HOME)/QMC/qmc/lib -laux
endif

# I compile on mck-login1, Richard on java.
ifeq ($(hostname),mck-login1)
# ifeq ($(hostname),java)
# Intel Itanium cluster at OSC using MPI and Intel Compiler efc
# LIBHOME = $(QMCHOME)
  FC = efc
  FFLAGS = -FI -extend_source -O2 -w -r8
  LDFLAGS = -Vaxlib -L/usr/local/mkl-6.0-beta/lib -lmkl_vml_itp
  FC_MPI = mpif77
  FFLAGS_MPI = -FI -extend_source -O2 -w -r8
  LDFLAGS_MPI = -Vaxlib -L/usr/local/mkl-6.0-beta/lib -lmkl_vml_itp
  BLAS = $(MKL) -L$(HOME)/QMC/qmc/lib -laux
endif

URI_CLUSTER = false
ifeq ($(hostname),mi)
  URI_CLUSTER = true
endif
ifeq ($(hostname),do)
  URI_CLUSTER = true
endif
ifeq ($(hostname),marx)
  URI_CLUSTER = true
endif
# -ax{i|M|K|W} i=Pentium Pro, M-Pentium MMX, K=Pentium3, W=Pentium4,Xeon
ifeq ($(URI_CLUSTER),true)
  FC = ifc
  FFLAGS = -extend_source -fpp -O2 -w -pad -axi
# -CA -CB -CU -g
# LDFLAGS = -Vaxlib -L/usr/local/mkl-6.0-beta/lib -lmkl_vml_itp
  LDFLAGS = -Vaxlib
  FFLAGS_NOOPT = -c -extend_source -fpp -O0 -w -pad -axi
  FC_MPI = mpif77
  FLINKER_MPI = mpif77
  FFLAGS_MPI = -extend_source -fpp -O2 -w -pad -axi
# LDFLAGS_MPI = -Vaxlib -L/usr/local/mkl-6.0-beta/lib -lmkl_vml_itp
  LDFLAGS_MPI = -Vaxlib
endif

endif

# Not sure if these are correct for SunOS
#ifeq ($(system),SunOS)
#ifeq ($(hostname),sasn100)
#  FC    = pgf77
#  FFLAGS = -fast -r8 -Mextend -Mdalign -Msave -Mbounds
## FFLAGS = -fast -r8 -Mextend -Mdalign -Msave
#  FC_MPI = pgf90
#  FFLAGS_MPI = -O -DMPI -DAIX -YEXT_NAMES=LCS -I/usr/local/lam/include
#else
#  FC = f77
#  FFLAGS = -O2 -stackvar
#endif
#endif

# If hostname is sasn100 then cross-compile at Sandia on SunOS for Linux
# for pgf compiler at Sandia, use cif90 rather than cif77 because cif77
# does not understand automatic arrays
ifeq ($(system),SunOS)
ifeq ($(hostname),sasn100)
# FC    = cif90
  FC    = cif77
# FFLAGS = -O2
# FFLAGS = -g -C
  FFLAGS = -fast -r8 -Mextend -Mdalign -Msave -Mbounds
# FFLAGS = -fast -r8 -Mextend -Mdalign -Msave
# FC_MPI = cif90
  FC_MPI = cif77
  FFLAGS_MPI = -fast -r8 -Mextend -Mdalign -Msave -Mbounds -lmpi
else
  FC = f77
  FFLAGS = -O2 -stackvar
endif
endif

ifeq ($(system),OSF1)
# The -fpe1 tells the compiler how to deal with floating point exceptions
ifeq ($(hostname),eigernorth.sandia.gov)
# FC = f77
# FC = $(CPLANTHOME)/bin/f90
  FC = /usr/local/cplant/west/current/bin/f77
# FFLAGS = -g -extend_source -align dcommons -fpe1
  FFLAGS = -O -extend_source -align dcommons -fpe1
  FC_MPI = /usr/local/cplant/west/current/bin/f77
# FFLAGS_MPI = -g -extend_source -align dcommons -fpe1
  FFLAGS_MPI = -O -extend_source -align dcommons -fpe1
endif
ifeq ($(hostname),eigersouth.sandia.gov)
# FC = f77
# FC = $(CPLANTHOME)/bin/f90
  FC = /usr/local/cplant/west/current/bin/f77
# FFLAGS = -g -extend_source -align dcommons -fpe1
  FFLAGS = -O -extend_source -align dcommons -fpe1
  FC_MPI = /usr/local/cplant/west/current/bin/f77
# FFLAGS_MPI = -g -extend_source -align dcommons -fpe1
  FFLAGS_MPI = -O -extend_source -align dcommons -fpe1
endif
endif

#ifeq ($(system),CYGWIN_NT-5.0)
#  FC = f77.exe
#  FFLAGS = -unix -fpp:"/n"
#  FC_MPI = f77.exe
#  FLINKER_MPI = link.exe
#  FFLAGS_MPI = -unix -fpp:"/n" -DCLOBBER -DMPI -Ic:/Progra~1/MPIPro/include
#  AR = ar
#endif
#ifeq ($(system),CYGWIN_NT-5.0)
#  FC = g77
#  FFLAGS =  -fcase-upper -fno-underscoring -enable-stdcall-fixup -O3 -malign-double -ffixed-line-length-none -fno-automatic -L/cygdrive/c/progra~1/mifd68~1/Lib
#  FC_MPI = g77
#  FFLAGS_MPI = -fcase-upper -fno-underscoring -enable-stdcall-fixup -O2 -malign-double -ffixed-line-length-none -fno-automatic -Wall -fbounds-check -I/cygdrive/c/PROGRA~1/MPIPro/include
#  LDFLAGS_MPI = -L/cygdrive/c/PROGRA~1/MPIPro/lib -lMPIPro_stdc -lMPIPro -L/cygdrive/c/progra~1/mifd68~1/Lib
#endif
# -4Yportlib is needed while compiling and linking for flush(6) to work.  It still does not
#  work unless I have "use iflport" in the program because we use -Qlowercase and the library
#  has capital FLUSH.
ifeq ($(system),CYGWIN_NT-5.0)
  FC = win32fe ifl
  FFLAGS = -FI -4L132 -Qlowercase -w -4Yportlib -O2
# FFLAGS = -Qlowercase -w -O2 -4L132 /Qoption,link,/verbose
# BLAS = $(MKL) -lderf -Vaxlib
# LDFLAGS = -L/cygdrive/c/progra~1/mifd68~1/Lib /Qoption,link,/verbose mkl_c.lib mkl_p3.lib
  FC_MPI = win32fe ifl
# Used to use MPIPro164 for v2 and MPIPro for v1,vplus and change library names.  Now MPIPro for all.
# Another way to deal with spaces in names is "Program~s" instead of "Program\ Files"
  FFLAGS_MPI = -FI -4L132 -Qlowercase -w -O2 -I/cygdrive/c/Program\ Files/MPIPro/include
  LDFLAGS_MPI = -Vaxlib -L/cygdrive/c/Program\ Files/MPIPro/lib MPIPro.lib MPIPro_cdec.lib
#-L/cygdrive/c/progra~1/mifd68~1/Lib mkl_c.lib mkl_p3.lib
endif

# environment variable for LAM-MPI
LAMHF77 = $(FC)
