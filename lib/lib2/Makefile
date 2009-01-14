# This Makefile is to be used with gmake.  On Linux systems, gmake is
# the default make.

#############################################################
# System dependent definitions (adjust to suit your system) #
#############################################################

AR = ar
ARFLAGS = r

system := $(shell uname)
hostname := $(shell uname -n)
# planck_compiler=g77
  planck_compiler=ifort
export FC FFLAGS AR ARFLAGS

# FC:         Sequential FC compile command
# FFLAGS:      Flags for sequential FC compilations
# FC_MPI:     Parallel FC compile command (may e.g. the F90 sequential compiler)
# FLINKER_MPI: Parallel link editor command (is useful if FC_MPI for some reason
#                is not an MPI compiler)
# FFLAGS_MPI:  Flags for parallel FC compilations (May have to contain include
#                directory that contains mpif.h)
# In most cases the Fortran compilers used for sequential and MPI programs are
# in fact the same, but many implementations of MPI provide a script (often called
# something like mpif77 of hf77) to call the Fortran compiler with the correct
# options, libraries and includes for MPI programs.  If your script does not have
# an include to mpif.h then you will need to put the path of your system's mpif.h
# INCLUDE_MPI variable.  Typically the system mpif.h is in a directory such
# as /usr/local/mpi/include.

ifeq ($(findstring IRIX,$(system)),IRIX)
  FC = f77
  FFLAGS = -O3 -n32 -mips4 -DIRIX
# FFLAGS = -g -C -n32 -mips4 -DIRIX
  FC_MPI = mpif77
  FFLAGS_MPI = -O3 -n32 -mips4 -DIRIX -DMPI -DNON_STANDARD
# FFLAGS_MPI = -g -C -n32 -mips4 -DIRIX -DMPI -DNON_STANDARD
endif

# This was originally for the Cornell SP2, later modified for the NERSC SP3, compile node Seaborg.
ifeq ($(system),AIX)
  FC = xlf
  FFLAGS = -O2 -bmaxdata:0x70000000 -qnosave -qsigtrap -WF,-DAIX -qfixed=132
  FC_MPI = mpxlf
  FFLAGS_MPI = -O2 -qfixed=132 -bmaxdata:0x70000000 -qnosave -qsigtrap -WF,-DAIX,-DMPI,-DNON_STANDARD
endif

ifeq ($(system),Linux)

# Linux default is now ifort
  FC = ifort
# FFLAGS = -g -CB -d4 -zero -save -extend_source -w -r8 -fpp2 -O0 -ax
# FFLAGS = -zero -save -extend_source -w -r8 -fpp2 -O2 -ax
# FFLAGS = -zero -save -extend_source -w -r8 -O3 -axN -ipo -pad
  FFLAGS = -zero -save -extend_source -w -r8 -O3 -pad
# LDFLAGS = -Vaxlib -static
  LDFLAGS = -static -nothread
  FC_MPI = mpif77 -fc=ifort
# FFLAGS_MPI = -g -CB -zero -save -extend_source -w -r8 -O0 -ax
# FFLAGS_MPI = -g -CB -zero -save -extend_source -w -r8 -O2 -ax
# FFLAGS_MPI = -zero -save -extend_source -w -r8 -O2 -ax
  FFLAGS_MPI = -zero -save -extend_source -w -r8 -O3 -ax -pad
# LDFLAGS_MPI = -Vaxlib -static
  LDFLAGS_MPI = -static -nothread -lmpi

# Old AMD Athlon cluster of Arias group
ifeq ($(hostname),dft.physics.cornell.edu)
ifeq ($(planck_compiler),g77)
# FC = g77-3.2
# FC = g77-3.0
  FC = g77
# FFLAGS = -g -C -malign-double -ffixed-line-length-none -fno-automatic -Wall -fbounds-check
# FFLAGS = -O2 -malign-double -ffixed-line-length-none -fno-automatic -Wall -fbounds-check
  FFLAGS = -O3 -malign-double -ffixed-line-length-none -fno-automatic
# FFLAGS = -pg -O3 -malign-double -ffixed-line-length-none -fno-automatic
endif

ifeq ($(planck_compiler),ifort)
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
endif
endif

ifeq ($(hostname),gauss)
  #-zero Initializes to zero all local scalar variables of intrinsic type INTEGER, REAL, COMPLEX, or LOGICAL that are saved but not yet initialized
#-save ????
#-extend-source Determines the column that ends the statement field of each source line in a fixed-format file (default 132). 
#-w disable all warnings
#-r8 set default size of real to 8 bytes
#-O3 aggressive optimization
#-pad enable changing variable and arrau memory layout
#-mtune=core2 optimize specifically for core2 cpu
#-march=core2 generate code exclusively for core2cpu
#-xT generate code to run exclusively on Intel(R) Core(TM)2 processor family with SSSE3
  FFLAGS = -zero -save -extend_source -w -r8 -O3 -pad -mtune=core2 -march=core2 -xT
#-static prevents linking with shared libraries
#-nothread
#-no multithreaded libraries
  LDFLAGS = -static -nothread
  FC_MPI = /usr/bin/mpif90 -fc=ifort
  FFLAGS_MPI = -zero -save -extend_source -w -r8 -O3  -pad -mtune=core2 -march=core2 -xT
  LDFLAGS_MPI = -static -nothread -lmpi
endif

ifeq ($(hostname),nanolab.cnf.cornell.edu)
  FC = ifort
# FFLAGS = -g -CB -d4 -zero -save -extend_source -w -r8 -fpp2 -O0 -ax
# FFLAGS = -zero -save -extend_source -w -r8 -fpp2 -O2 -ax
# FFLAGS = -zero -save -extend_source -w -r8 -O3 -axN -ipo -pad
  FFLAGS = -zero -save -extend_source -w -r8 -O3 -pad
# LDFLAGS = -Vaxlib -static
  LDFLAGS = -static -nothread
  FC_MPI = /usr/local/bin/mpif77 -fc=ifort
# FFLAGS_MPI = -g -CB -zero -save -extend_source -w -r8 -O0 -ax
# FFLAGS_MPI = -g -CB -zero -save -extend_source -w -r8 -O2 -ax
# FFLAGS_MPI = -zero -save -extend_source -w -r8 -O2 -ax
  FFLAGS_MPI = -zero -save -extend_source -w -r8 -O3 -ax -pad
# LDFLAGS_MPI = -Vaxlib -static
  LDFLAGS_MPI = -static -nothread -lmpi
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
ifeq ($(hostname),amd-login1)
# AMD Cluster at OSC using Portland Group Compiler
  FC = pgf77
  FFLAGS = -fast -r8 -Mextend -Mdalign -Msave -Mbounds
  FC_MPI = mpif90
  FFLAGS_MPI = -O -DMPI -DAIX -YEXT_NAMES=LCS -I/usr/local/lam/include
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

# Argonne Blue-Gene
ifeq ($(hostname),login1)
# for debug front-end
  FC = xlf
  FFLAGS = -C -g -qlinedebug -qsigtrap -qfixed=132

# for nodes
  FC = blrts_xlf
  FFLAGS = -O5 -qarch=440 -qtune=440  -qfixed=132
# FFLAGS = -O2 -qnosave -qsigtrap -WF,-DAIX -qfixed=132
# FFLAGS = -g -O -qarch=440 -qtune=440 -qmaxmem=64000 -qnosave -qsigtrap -WF,-DAIX -qfixed=132
# FFLAGS = -g -O5 -qarch=440 -qtune=440 -qnosave -qsigtrap -WF,-DAIX -qfixed=132
# FFLAGS = -g -O5 -qarch=440d -qtune=440 -qnosave -qsigtrap -WF,-DAIX -qfixed=132
  LDFLAGS =
# FC_MPI = mpxlf
  FC_MPI = mpif90.ibm
# FC_MPI = blrts_xlf
# FFLAGS_MPI = -O2 -qnosave -qsigtrap -WF,-DAIX -qfixed=132
  FFLAGS_MPI = -g -O -qarch=440 -qtune=440 -qmaxmem=64000 -qnosave -qsigtrap -WF,-DAIX -qfixed=132
# FFLAGS_MPI = -g -O5 -qarch=440 -qtune=440 -qnosave -qsigtrap -WF,-DAIX -qfixed=132
# FFLAGS_MPI = -g -O5 -qarch=440d -qtune=440 -qnosave -qsigtrap -WF,-DAIX -qfixed=132
  LDFLAGS_MPI =

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
  LDFLAGS = -Vaxlib -L/usr/local/mkl-6.0-beta/lib -lmkl_vml_itp
  FFLAGS_NOOPT = -c -extend_source -fpp -O0 -w -pad -axi
  FC_MPI = mpif77
  FLINKER_MPI = mpif77
  FFLAGS_MPI = -DMPI -DNON_STANDARD -DAIX -fpp -w
  LDFLAGS_MPI = -Vaxlib -L/usr/local/mkl-6.0-beta/lib -lmkl_vml_itp
endif
endif

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

# CRAY X1 at OSC armstrong.sf.osc.edu
# real*8 numbers ranging from 1d-323 to 1d+308
# real*4 numbers ranging from 1e-44 to 1e+38
ifeq ($(hostname),armstrong)
  CYRUSLIB= -L$(LIBHOME)/lib -lcyrus
  BLAS    = -L$(LIBHOME)/lib2/blas -lblas
  LINPACK = -L$(LIBHOME)/lib2/linpack -llinpack
  QUENCH  = -L$(LIBHOME)/SimulatedAnnealing/quench_anneal/lib -lquench -lquench_seq
  PSPLINE = -L$(LIBHOME)/lib2/pspline/pspline -lpspline
  FC = ftn
  FFLAGS = -O2 -Rabc -Ossp -ev -N 132
  FC_MPI = ftn
  FFLAGS_MPI = -O2 -DMPI -Ossp -ev -N 132
endif

ifeq ($(system),OSF1)
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
#  FFLAGS_MPI = -fcase-upper -fno-underscoring -enable-stdcall-fixup -O2 -malign-double -ffixed-line-length-none -fno-automatic -Wall -fbounds-check -I/cygdrive/c/PROGRA~1/MPIPro164/include
#  LDFLAGS_MPI = -L/cygdrive/c/PROGRA~1/MPIPro164/lib -lMPIPro_stdc -lMPIPro -L/cygdrive/c/progra~1/mifd68~1/Lib
#endif
# -4Yportlib is needed while compiling and linking for flush(6) to work.  It still does not
#  work because we use -Qlowercase and the library has capital FLUSH.
ifeq ($(system),CYGWIN_NT-5.0)
  FC = win32fe ifl
  FFLAGS = -FI -4L132 -Qlowercase -w -4Yportlib -O2
# FFLAGS = -Qlowercase -w -O2 -4L132 /Qoption,link,/verbose
# BLAS = $(MKL) -lderf -Vaxlib
# LDFLAGS = -L/cygdrive/c/progra~1/mifd68~1/Lib /Qoption,link,/verbose mkl_c.lib mkl_p3.lib
  FC_MPI = win32fe ifl
# MPIPro164 for v2.  For v1,vplus use MPIPro and change library names
# Another way to deal with spaces in names is "Program~s" instead of "Program\ Files"
  FFLAGS_MPI = -FI -4L132 -Qlowercase -w -O2 -I/cygdrive/c/Program\ Files/MPIPro164/include
  LDFLAGS_MPI = -Vaxlib -L/cygdrive/c/Program\ Files/MPIPro164/lib MPIPro.lib MPIPro_cdec.lib
#-L/cygdrive/c/progra~1/mifd68~1/Lib mkl_c.lib mkl_p3.lib
endif

# Jaguar at NCCS Oak Ridge
ifeq ($(findstring jaguar,$(hostname)), jaguar)

# Portland compiler
#  FC    = pgf95
#  FFLAGS = -fast -r8 -Mextend -Mdalign -Msave
#  FC_MPI = mpif77.pgf
#  FFLAGS_MPI = -fast -r8 -DCLOBBER -DMPI

# Pathscale compiler
# module swap PrgEnv-pgi PrgEnv-pathscale && module load acml
  FC = ftn -target=linux
  FFLAGS = -O3 -r8 -extend-source -zerouv  -fno-second-underscore
  LDFLAGS = -O3 -r8 -extend-source -zerouv  -fno-second-underscore
  FC_MPI = ftn -target=linux
  FFLAGS_MPI = -O3 -r8 -DCLOBBER -DMPI -zerouv -fno-second-underscore

endif

# clastos in Paris
ifeq ($(hostname),clastos)
  FC = /usr/local/bin/mpif90_64
  FC_MPI = /usr/local/bin/mpif90_64
endif

# AMD Opteron cluster with Pathscale compiler at Lyngby, Denmark
ifeq ($(findstring slid.fysik.dtu.dk,$(hostname)), slid.fysik.dtu.dk)

 FC = pathf90
#FFLAGS = -cpp -DPATHSCALE -O2 -r8 -extend-source -zerouv
 FFLAGS = -O3 -r8 -extend-source -zerouv  -fno-second-underscore
 FC_MPI = pathf90
 FFLAGS_MPI = -O3 -r8 -DCLOBBER -DMPI -zerouv -fno-second-underscore

endif

# Julich BlueGene P
ifeq ($(hostname),jugene1)

  FC = mpixlf77
  FFLAGS = -O3 -qstrict -qarch=450 -qtune=450 -qsuffix=cpp=f -WF,-DAIX -qfixed=132 
  LDFLAGS = 
  FC_MPI = mpixlf77      
  FFLAGS_MPI = -O3 -qstrict -qarch=450 -qtune=450 -qsuffix=cpp=f -WF,-DAIX -qfixed=132 
  LDFLAGS_MPI =

endif

###################################
# Targets, Dependencies and Rules #
###################################

# PHONY is a Make keyword that tells it to execute the rule regardless of
# the dependencies.  We use it here for 2 reasons:
# 1) prevents it from getting confused if there happens to be a file named
#    "clean" or "clean_all" in the directory.
# 2) Also, for the libraries, linpack etc. it does the make even though there
#    are directories named linpack etc.
.PHONY: linpack lapack blas pspline clean clean_all

all: linpack lapack blas pspline
# all: linpack blas pspline

linpack:
	cd linpack ; $(MAKE) liblinpack.a

lapack:
	cd lapack ; $(MAKE) liblapack.a

blas:
	cd blas ; $(MAKE) libblas.a

pspline:
	cd pspline/pspline ; $(MAKE) libpspline.a

clean:
	-rm -f linpack/*.o lapack/*.o blas/*.o pspline/pspline/*.o
clean_all:
	-rm -f linpack/*.o lapack/*.o blas/*.o pspline/pspline/*.o \
            linpack/liblinpack.a lapack/liblapack.a pspline/pspline/libpspline.a
