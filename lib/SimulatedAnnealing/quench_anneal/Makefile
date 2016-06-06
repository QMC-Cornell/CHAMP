# I do not know why but when this is called from champ/qmc/Makefile and we are compiling there
# in MPI mode, using "make mpi" then it insists on using mpif90 instead of resetting the compiler
# to gfortran.  Anyway, this does not seem to matter but we need the "FFLAGS = -DNOMPI" to prevent
# quench_anneal from being compiled in MPI mode when CHAMP is being compiled in MPI mode.
# Alternatively one can do a "make libs" from here before doing a "make mpi" in champ/qmc.

FC = gfortran
F77 = gfortran
F77_MPI = gfortran
F95 = gfortran
F95_MPI = gfortran
#FFLAGS=-finit-local-zero -O3 -fdefault-real-8 -ffixed-line-length-132 -x f77-cpp-input -Wall -fbounds-check
FFLAGS=-finit-local-zero -O3 -ffixed-line-length-132 -x f77-cpp-input -Wall -fbounds-check
# FFLAGS=-finit-local-zero -O3 -fdefault-real-8 -ffixed-line-length-132 -x f77-cpp-input -Wall
export FC F77 FFLAGS LDFLAGS F77_MPI FFLAGS_MPI FLINKER_MPI
 
# ifdef F77
# F77 = gfortran
# FFLAGS = -finit-local-zero -O3 -ffixed-line-length-132 -x f77-cpp-input -Wall -fbounds-check -DNOMPI
# # FC = gfortran
# # F77_MPI = gfortran
# # F95 = gfortran
# # F95_MPI = gfortran
# endif
# 
# ifndef F77

##########
# README #
##########

# This Makefile is to be used with gmake.  On Linux systems, gmake is
# the default make.
# Some versions of gmake give error messages at lines starting with:
# $(shell cp -p parameters/

# (i) If you are working with one of the following 3 operating systems,
# IRIX, AIX, or Unix, and if you wish to use the standard fortran 77 compiler
# that comes with these operating systems then you need not make this change.
# If you are not using one of these operating systems, you need to add lines
# below to create a system-dependent include/parameters.h file and define F77,
# FFLAGS, F77_MPI and FFLAGS_MPI for the new operating system, using the lines
# below as a model.  Also, even if you are using one of these operating systems
# but with a different compiler or different compiler options you will need to
# either modify the existing lines or add new ones.

# (For most compilers, preprocessor constants may be defined using the -D option.)

# (ii) MPI is a preprocessor flag used in main/func_mpi.F.  If MPI is
# defined, the parts of main/func_mpi.F that use MPI calls are compiled;
# if MPI is undefined, these parts are omitted.  The variable FFLAGS_MPI
# must always include a definition of MPI; the variable FFLAGS must
# *not* define MPI.  The MPI preprocessor flag is also used in tests/test.F.

# (iii) NON_STANDARD is a preprocessor flag used to include code that is
# non-standard and will not compile on all systems.  At present it is used
# to include subroutine allgatherv_2d2, allgatherv_2d1 in main/func_mpi.F.
# These routines use pointers, so if you are using a compiler that does not
# support pointers, you will be unable to use the MPI version.  This could
# be easily remedied though.

# (iv) CLOBBER is a preprocessor flag used to indicate if it is OK to open
# files using "status=old" instead of "status=read" or "status=readonly" on
# systems where we do not know how to open the files in read only mode.

#############################################################
# System dependent definitions (adjust to suit your system) #
#############################################################

system := $(shell uname)
hostname := $(shell uname -n)
# planck_compiler=g77
planck_compiler=ifort
$(shell mkdir -p lib)
# system = AIX
export FC F77 FFLAGS LDFLAGS F77_MPI FFLAGS_MPI FLINKER_MPI

# F77:         Sequential F77 compile command
# FFLAGS:      Flags for sequential F77 compilations
# F77_MPI:     Parallel F77 compile command (may e.g. the F90 sequential compiler)
# FLINKER_MPI: Parallel link editor command (is useful if F77_MPI for some reason
#                is not an MPI compiler)
# FFLAGS_MPI:  Flags for parallel F77 compilations (May have to contain include
#                directory that contains mpif.h)
# In most cases the Fortran compilers used for sequential and MPI programs are
# in fact the same, but many implementations of MPI provide a script (often called
# something like mpif77 of hf77) to call the Fortran compiler with the correct
# options, libraries and includes for MPI programs.  If your script does not have
# an include to mpif.h then you will need to put the path of your system's mpif.h
# INCLUDE_MPI variable.  Typically the system mpif.h is in a directory such
# as /usr/local/mpi/include.


ifeq ($(system),AIX)
  $(shell cp -p parameters/parameters.AIX include/parameters.h)
  F77 = xlf
  FFLAGS = -O2 -qnosave -qsigtrap -WF,-DAIX
  F77_MPI = mpxlf
  FLINKER_MPI = mpxlf
  FFLAGS_MPI = -O2 -qnosave -qsigtrap -WF,-DAIX,-DMPI,-DNON_STANDARD
endif

ifeq ($(system),Linux)
  $(shell cp -p parameters/parameters.Linux include/parameters.h)
# Linux default is now ifort
  F77 = ifort
# FFLAGS = -g -CB -d4 -zero -save -extend_source -w -r8 -fpp2 -O0 -ax
# FFLAGS = -zero -save -extend_source -w -r8 -fpp2 -O2 -ax
# FFLAGS = -zero -save -extend_source -w -r8 -O3 -axN -ipo -pad
  FFLAGS = -zero -save -extend_source -w -r8 -O3 -pad
# LDFLAGS = -Vaxlib -static
  LDFLAGS = -static -nothread
  F77_MPI = mpif77 -fc=ifort
# FFLAGS_MPI = -g -CB -zero -save -extend_source -w -r8 -O0 -ax
# FFLAGS_MPI = -g -CB -zero -save -extend_source -w -r8 -O2 -ax
# FFLAGS_MPI = -zero -save -extend_source -w -r8 -O2 -ax
  FFLAGS_MPI = -zero -save -extend_source -w -r8 -O3 -ax -pad
# LDFLAGS_MPI = -Vaxlib -static
  LDFLAGS_MPI = -static -nothread -lmpi

# Old AMD Athlon cluster of Arias group
ifeq ($(hostname),dft.physics.cornell.edu)
ifeq ($(planck_compiler),g77)
# to link g77 mpi libraries, make sure /usr/bin is ahead of /usr/local/bin
# F77 = g77-3.2
  F77 = g77-3.0
# F77 = g77
# FFLAGS = -g -C -malign-double -ffixed-line-length-none -fno-automatic -Wall -fbounds-check
# FFLAGS = -O2 -malign-double -ffixed-line-length-none -fno-automatic -Wall -fbounds-check
  FFLAGS = -O3 -malign-double -ffixed-line-length-none -fno-automatic
# FFLAGS = -pg -O3 -malign-double -ffixed-line-length-none -fno-automatic
# LDFLAGS = -static
  F77_MPI = mpif77 -fc=g77
# FFLAGS_MPI = -O2 -malign-double -ffixed-line-length-none -fno-automatic -Wall -fbounds-check -DCLOBBER -DMPI
  FFLAGS_MPI = -O3 -malign-double -ffixed-line-length-none -fno-automatic -DCLOBBER -DMPI
# LDFLAGS_MPI = -static
  FLINKER_MPI = mpif77
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
  F77 = ifort
# FFLAGS = -g -C -d4 -extend_source -w -r8 -fpp2 -O0 -axi
# FFLAGS = -extend_source -w -r8 -fpp2 -O2 -axi
# FFLAGS = -extend_source -w -r8 -O3 -axiN -ipo -pad
  FFLAGS = -extend_source -w -r8 -O3 -axi -pad
# LDFLAGS = -Vaxlib -static
  LDFLAGS = -static
  F77_MPI = mpif77 -fc=ifort
# FFLAGS_MPI = -g -C -extend_source -w -r8 -O2 -axi
# FFLAGS_MPI = -extend_source -w -r8 -O2 -axi
  FFLAGS_MPI = -extend_source -w -r8 -O3 -axi -pad -DCLOBBER -DMPI
# LDFLAGS_MPI = -Vaxlib -static
  LDFLAGS_MPI = -static -lmpi
  FLINKER_MPI = mpif77
endif
endif

#Gauss in new Arias cluster
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
  F77_MPI = /usr/bin/mpif90.openmpi -fc=ifort
  FFLAGS_MPI = -zero -save -extend_source -w -r8 -O3  -pad -mtune=core2 -march=core2 -xT
  LDFLAGS_MPI = -static -nothread -lmpi
  FLINKER_MPI = /usr/bin/mpif90.openmpi
endif

#new nodes at cornell
ifeq ($(findstring node,$(hostname)), node)
  F77 = gfortran
  FFLAGS = -finit-local-zero -O3 -fdefault-real-8 -ffixed-line-length-132 -Wall
  F77_MPI = gfortran
  FFLAGS_MPI = -finit-local-zero -O3 -fdefault-real-8 -ffixed-line-length-132 -Wall
endif

#another new node at cornell
ifeq ($(hostname),qmc8)
  F77 = gfortran
  FFLAGS = -finit-local-zero -O3 -fdefault-real-8 -ffixed-line-length-132 -Wall
  F77_MPI = gfortran
  FFLAGS_MPI = -finit-local-zero -O3 -fdefault-real-8 -ffixed-line-length-132 -Wall
endif

ifeq ($(hostname),nanolab.cnf.cornell.edu)
  F77 = ifort
# FFLAGS = -g -CB -d4 -zero -save -extend_source -w -r8 -fpp2 -O0 -ax
# FFLAGS = -zero -save -extend_source -w -r8 -fpp2 -O2 -ax
# FFLAGS = -zero -save -extend_source -w -r8 -O3 -axN -ipo -pad
  FFLAGS = -zero -save -extend_source -w -r8 -O3 -pad
# LDFLAGS = -Vaxlib -static
  LDFLAGS = -static -nothread
  F77_MPI = mpif77
# FFLAGS_MPI = -g -CB -zero -save -extend_source -w -r8 -O0 -ax
# FFLAGS_MPI = -g -CB -zero -save -extend_source -w -r8 -O2 -ax
# FFLAGS_MPI = -zero -save -extend_source -w -r8 -O2 -ax
#  FFLAGS_MPI = -zero -save -extend_source -w -r8 -O3 -ax -pad
  FFLAGS_MPI = -O3
# LDFLAGS_MPI = -Vaxlib -static
  LDFLAGS_MPI = -static -nothread -lmpi
endif

ifeq ($(hostname),tyche.acc.uri.edu)
  F77 = f90
  FFLAGS = -O -YEXT_NAMES=UCS
  F77_MPI = f90
  FLINKER_MPI = mpif77
  FFLAGS_MPI = -O -DMPI -DAIX -DNON_STANDARD -YEXT_NAMES=UCS -I$(HOME)/mpich/include
endif
ifeq ($(hostname),charming)
  F77 = pgf90
  FFLAGS = -O2
  FFLAGS = -fast -r8
  F77_MPI = mpif90.pgf
  FLINKER_MPI = mpif90.pgf
  FFLAGS_MPI = -fast -r8 -DCLOBBER -DMPI
endif
ifeq ($(hostname),mesdag)
  F77 = pgf90
  FFLAGS = -O2
  FFLAGS = -fast -r8
  F77_MPI = mpif90.pgf
  FLINKER_MPI = mpif90.pgf
  FFLAGS_MPI = -fast -r8 -DCLOBBER -DMPI
endif
# ifeq ($(hostname),oscbw01)
# ifeq ($(hostname),amd-login1)
# AMD Cluster at OSC using Portland Group Compiler
#   F77 = pgf77
#   FFLAGS = -fast -r8 -Mextend -Mdalign -Msave -Mbounds
#   F77_MPI = mpif77
##  FFLAGS_MPI = -fast -r8 -Mextend -Mdalign -Msave -Mbounds
#   FFLAGS_MPI = -fast -r8 -Mextend -Mdalign -Msave
# endif
ifeq ($(hostname),amd-login1)
  F77 = ifc
  FFLAGS = -FI -extend_source -O2 -w -r8 -fpp2
  LDFLAGS = -Vaxlib
  F77_MPI = mpif77 -fc=ifc
  FFLAGS_MPI = -FI -extend_source -O2 -w -r8
  LDFLAGS_MPI = -Vaxlib
endif
ifeq ($(hostname),ia64)
# Intel Itanium cluster at OSC using MPI and Intel Compiler efc
# LIBHOME = $(QMCHOME)
  F77 = efc
# FFLAGS = -FI -extend_source -O2 -w -r8 -Dderf=INTELERF -Dderfc=INTELERF77 -fpp2
  FFLAGS = -FI -extend_source -O2 -w -r8
  F77_MPI = mpif77
  FFLAGS_MPI = -FI -extend_source -O2 -w -r8
# BLAS = $(MKL) -lderf -Vaxlib
  BLAS = $(MKL) -L$(HOME)/QMC/qmc/lib -laux -Vaxlib
endif
ifeq ($(hostname),mck-login1)
# Intel Itanium cluster at OSC using MPI and Intel Compiler efc
# LIBHOME = $(QMCHOME)
  F77 = efc
  FFLAGS = -FI -extend_source -O2 -w -r8
  F77_MPI = mpif77
  FFLAGS_MPI = -FI -extend_source -O2 -w -r8
  BLAS = $(MKL) -L$(HOME)/QMC/qmc/lib -laux -Vaxlib
endif

# Argonne Blue-Gene
ifeq ($(hostname),login1)

# for debug front-end
  F77 = xlf
  FFLAGS = -C -g -qlinedebug -qsigtrap -qfixed=132

# for nodes
  F77 = blrts_xlf
  FFLAGS = -O5 -qarch=440 -qtune=440 -WF,-DAIX -qfixed=132
#  FFLAGS = -g -O -qarch=440 -qtune=440 -qmaxmem=64000 -qnosave -qsigtrap -WF,-DAIX -qfixed=132
# FFLAGS = -O2 -qnosave -qsigtrap -WF,-DAIX -qfixed=132
# FFLAGS = -g -O5 -qarch=440 -qtune=440 -qnosave -qsigtrap -WF,-DAIX -qfixed=132
# FFLAGS = -g -O5 -qarch=440d -qtune=440 -qnosave -qsigtrap -WF,-DAIX -qfixed=132
  LDFLAGS =
# F77_MPI = mpxlf
  F77_MPI = mpif90.ibm
# F77_MPI = blrts_xlf
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
  F77 = mpif77
  FFLAGS = -c -extend_source -fpp -O2 -w -pad -axi
# -CA -CB -CU -g
  FFLAGS_NOOPT = -c -extend_source -fpp -O0 -w -pad -axi
  F77_MPI = mpif77
  FLINKER_MPI = mpif77
  FFLAGS_MPI = -DMPI -DNON_STANDARD -DAIX -fpp -w
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
  F77 = ftn
  FFLAGS = -O2 -Rabc -Ossp -ev -N 132
  F77_MPI = ftn
  FFLAGS_MPI = -O2 -DMPI -Ossp -ev -N 132
endif


# Not sure if these are correct for SunOS
# If hostname is sasn100 then cross-compile at Sandia on SunOS for Linux
# for pgf compiler at Sandia, use cif90 rather than cif77 because cif77
# does not understand automatic arrays
# I have not checked if -DNON_STANDARD works with either compiler.
ifeq ($(system),SunOS)
ifeq ($(hostname),sasn100)
  $(shell cp -p parameters/parameters.Linux include/parameters.h)
  F77 = cif90
  FFLAGS = -O2
# FFLAGS = -g -C
  F77_MPI = cif90
  FLINKER_MPI = cif90
  FFLAGS_MPI = -O2 -DCLOBBER -DMPI -DNON_STANDARD
else
  $(shell cp -p parameters/parameters.IRIX_Solaris include/parameters.h)
  F77 = f77
  FFLAGS = -O2 -stackvar
  F77_MPI = mpif77
  FLINKER_MPI = mpif77
  FFLAGS_MPI = -O2 -stackvar -DCLOBBER -DMPI -DNON_STANDARD
endif
endif

ifeq ($(system),OSF1)
# DEC/Compaq compiler
# I have not checked if -DNON_STANDARD works
ifeq ($(hostname),eigernorth.sandia.gov)
  $(shell cp -p parameters/parameters.Linux include/parameters.h)
# F77 = $(CPLANTHOME)/bin/f90
  F77 = /usr/local/cplant/west/current/bin/f77
  FFLAGS = -g -extend_source -align dcommons -fpe1
  F77_MPI = /usr/local/cplant/west/current/bin/f77
  FFLAGS_MPI = -g -extend_source -align dcommons -fpe1 -DMPI -DNON_STANDARD
endif
ifeq ($(hostname),eigersouth.sandia.gov)
  $(shell cp -p parameters/parameters.Linux include/parameters.h)
# F77 = $(CPLANTHOME)/bin/f90
  F77 = /usr/local/cplant/west/current/bin/f77
  FFLAGS = -g -extend_source -align dcommons -fpe1
  F77_MPI = /usr/local/cplant/west/current/bin/f77
  FFLAGS_MPI = -g -extend_source -align dcommons -fpe1 -DMPI -DNON_STANDARD
endif
endif

#ifeq ($(system),CYGWIN_NT-5.0)
#  $(shell cp -p parameters/parameters.Linux include/parameters.h)
#  F77 = f77.exe
#  FFLAGS = -unix -fpp:"/n"
#  F77_MPI = f77.exe
#  FLINKER_MPI = link.exe
#  FFLAGS_MPI = -unix -fpp:"/n" -DCLOBBER -DMPI -Ic:/Progra~1/MPIPro/include
#  AR = ar
#endif
#ifeq ($(system),CYGWIN_NT-5.0)
# $(shell cp -p parameters/parameters.Linux include/parameters.h)
#  F77 = g77
#  FFLAGS =  -fcase-upper -fno-underscoring -enable-stdcall-fixup -O3 -malign-double -ffixed-line-length-none -fno-automatic -L/cygdrive/c/progra~1/mifd68~1/Lib
#  F77_MPI = g77
#  FFLAGS_MPI = -fcase-upper -fno-underscoring -enable-stdcall-fixup -O2 -malign-double -ffixed-line-length-none -fno-automatic -Wall -fbounds-check -I/cygdrive/c/PROGRA~1/MPIPro164/include
#  LDFLAGS_MPI = -L/cygdrive/c/PROGRA~1/MPIPro164/lib -lMPIPro_stdc -lMPIPro -L/cygdrive/c/progra~1/mifd68~1/Lib
#endif
ifeq ($(system),CYGWIN_NT-5.0)
  $(shell cp -p parameters/parameters.Linux include/parameters.h)
  F77 = win32fe ifl
  FFLAGS = -FI -4L132 -Qlowercase -w -O2 -Qfpp
# FFLAGS = -FI -4L132 -Qlowercase -w -O2 -Qfpp /Qoption,link,/verbose
# BLAS = $(MKL) -lderf -Vaxlib
# LDFLAGS = -L/cygdrive/c/progra~1/mifd68~1/Lib /Qoption,link,/verbose mkl_c.lib mkl_p3.lib
  F77_MPI = win32fe ifl
# MPIPro164 for v2.  For v1,vplus use MPIPro and change library names
# Another way to deal with spaces in names is "Program~s" instead of "Program\ Files"
  FFLAGS_MPI = -FI -4L132 -Qlowercase -w -O2 -Qfpp -I/cygdrive/c/Program\ Files/MPIPro164/include
  LDFLAGS_MPI = -L/cygdrive/c/Program\ Files/MPIPro164/lib MPIPro.lib MPIPro_cdec.lib
# -L/cygdrive/c/progra~1/mifd68~1/Lib mkl_c.lib mkl_p3.lib
endif

# Jaguar at NCCS Oak Ridge
ifeq ($(findstring jaguar,$(hostname)), jaguar)

# Portland compiler
#  F77    = pgf95
#  FFLAGS = -fast -r8 -Mextend -Mdalign -Msave
#  F77_MPI = mpif77.pgf
#  FFLAGS_MPI = -fast -r8 -DCLOBBER -DMPI

# Pathscale compiler
# module swap PrgEnv-pgi PrgEnv-pathscale && module load acml
  F77 = ftn -target=linux
  FFLAGS = -O3 -r8 -extend-source -zerouv  -fno-second-underscore
  LDFLAGS = -O3 -r8 -extend-source -zerouv  -fno-second-underscore
  F77_MPI = ftn -target=linux
  FFLAGS_MPI = -O3 -r8 -DCLOBBER -DMPI -zerouv -fno-second-underscore

endif

ifeq ($(hostname),clastos)
  F77 = /usr/local/bin/mpif90_64
  F77_MPI = /usr/local/bin/mpif90_64
  FLINKER_MPI = /usr/local/bin/mpif90_64
endif

#ifeq ($(hostname),hydrogene.lct.jussieu.fr)
#FC = ifort
#F77 = ifort
#F95 = ifort
#FFLAGS=-zero -r8 -pad -align all -O3 -ip -diag-disable vec -fpp -shared-intel
#endif

# AMD Opteron cluster with Pathscale compiler at Lyngby, Denmark
ifeq ($(findstring slid.fysik.dtu.dk,$(hostname)), slid.fysik.dtu.dk)

 F77 = pathf90
#FFLAGS = -cpp -DPATHSCALE -O2 -r8 -extend-source -zerouv
 FFLAGS = -O3 -r8 -extend-source -zerouv  -fno-second-underscore
 F77_MPI = pathf90
 FFLAGS_MPI = -O3 -r8 -DCLOBBER -DMPI -zerouv -fno-second-underscore

endif

# Julich BlueGene P
ifeq ($(hostname),jugene1)

  F77 = mpixlf77
  FFLAGS = -O3 -qstrict -qarch=450 -qtune=450 -qsuffix=cpp=F -WF,-DAIX -qfixed=132
  LDFLAGS =
  F77_MPI = mpixlf77
  FFLAGS_MPI = -O3 -qstrict -qarch=450 -qtune=450 -qsuffix=cpp=F -WF,-DAIX -qfixed=132
  LDFLAGS_MPI =

endif

#BM, personnal computer 1
ifeq ($(hostname),antares)
  F77=gfortran
  FFLAGS=-finit-local-zero -O3 -ffixed-line-length-132 -x f77-cpp-input -Wall -fbounds-check
endif

#BM, personnal computer 2
ifeq ($(hostname),antaresICS)
  F77=gfortran
  FFLAGS=-finit-local-zero -O3 -ffixed-line-length-132 -x f77-cpp-input -Wall -fbounds-check
endif

# The Intel compiler ifc causes trouble in dlamch unless the latter is
# compiled with -O0
# ifneq ($(F77),ifc)
#   FFLAGS_NOOPT = $(FFLAGS)
# endif
ifndef FFLAGS_NOOPT
  FFLAGS_NOOPT = $(FFLAGS)
endif

#endif # End of ifndef at top of Makefile

#######################
# Using this Makefile #
#######################

# In addition to the self-explanatory linpack, lapack, and blas archives,
# the quench-anneal library consists of three archive (.a) files.  These are:
# ./lib/libquench.a, which contains object files used by both the sequential
# and MPI versions of the library;
# ./lib/libquench_seq.a, which contains object files used only by the sequential
# version of the library; and
# ./lib/libquench_mpi.a, which contains object files used only by the MPI version
# of the library.
# The default make target, called "all", makes both the sequential and the
# mpi versions.  "gmake seq" makes ./lib/libquench.a and
# ./lib/libquench_seq.a.  It does not use MPI and does not require that the
# MPI library be installed.  To make the parallel version of the
# library, use "gmake mpi", which makes ./lib/libquench.a and
# ./lib/libquench_mpi.a.  (Of course this will not work unless MPI is
# already installed on your system.)  Note that running "make mpi" uses
# both sequential and MPI f77 compilers, so all four of the variables
# F77, FFLAGS, F77_MPI, FFLAGS_MPI, and FLINKER_MPI must be defined;
# also see INCLUDE_MPI above.

# On most systems the Makefile should not leave unwanted object files hanging
# around.  Note however that if you change a compiler option, then you need to
# first get rid of any libraries that already exist by typing "gmake clean_all"
# before remaking them, otherwise make will think they are up to date.

###################################################################
# The rest of this Makefile is system independent (We hope).  You #
# should not need to read or modify it.                           #
###################################################################


########################
# Variable definitions #
########################

# Shorthand names for library targets.
LIB_COMMON	= lib/libquench.a
LIB_SEQ 	= lib/libquench_seq.a
LIB_MPI 	= lib/libquench_mpi.a
LIB_BLAS	= lib/libblas.a
LIB_LAPACK	= lib/liblapack.a
LIB_LINPACK	= lib/liblinpack.a

# Includes
INC          = include/parameters.h include/dimen.h
INC_MPI	     = include/parameters.h include/dimen.h include/IO_stop.h

# Objects that are identical in MPI and sequential modes
OBJS         = $(LIB_COMMON)(main/anneal.o)\
               $(LIB_COMMON)(main/anneal_2.o)\
               $(LIB_COMMON)(main/derivs.o)
OBJS_MATH    = $(LIB_COMMON)(math/check_arith.o)\
               $(LIB_COMMON)(math/chlsky.o)\
               $(LIB_COMMON)(math/cholesky_piv.o)\
               $(LIB_COMMON)(math/gaushn.o)\
               $(LIB_COMMON)(math/lugaus_new.o)\
               $(LIB_COMMON)(math/lxb.o)\
               $(LIB_COMMON)(math/ransr.o)\
               $(LIB_COMMON)(math/svd_gaus_test.o)\
               $(LIB_COMMON)(math/uxb.o)
OBJS_BLAS	=\
		$(LIB_BLAS)(lib_standard/blas/daxpy.o)\
		$(LIB_BLAS)(lib_standard/blas/dcopy.o)\
		$(LIB_BLAS)(lib_standard/blas/ddot.o)\
		$(LIB_BLAS)(lib_standard/blas/dgemm.o)\
		$(LIB_BLAS)(lib_standard/blas/dgemv.o)\
		$(LIB_BLAS)(lib_standard/blas/dger.o)\
		$(LIB_BLAS)(lib_standard/blas/dnrm2.o)\
		$(LIB_BLAS)(lib_standard/blas/drot.o)\
		$(LIB_BLAS)(lib_standard/blas/drotg.o)\
		$(LIB_BLAS)(lib_standard/blas/dscal.o)\
		$(LIB_BLAS)(lib_standard/blas/dswap.o)\
		$(LIB_BLAS)(lib_standard/blas/dtrmm.o)\
		$(LIB_BLAS)(lib_standard/blas/dtrmv.o)\
		$(LIB_BLAS)(lib_standard/blas/lsame.o)\
		$(LIB_BLAS)(lib_standard/blas/xerbla.o)

OBJS_LINPACK	=$(LIB_LINPACK)(lib_standard/linpack/dsvdc.o)

OBJS_LAPACK	=\
		$(LIB_LAPACK)(lib_standard/lapack/dbdsqr.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dgebd2.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dgebrd.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dgelq2.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dgelqf.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dgeqr2.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dgeqrf.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dgesvd.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dlabrd.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dlacpy.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dlamch.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dlange.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dlapy2.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dlarf.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dlarfb.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dlarfg.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dlarft.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dlartg.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dlas2.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dlascl.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dlaset.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dlasq1.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dlasq2.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dlasq3.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dlasq4.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dlasq5.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dlasq6.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dlasr.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dlasrt.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dlassq.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dlasv2.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dorg2r.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dorgbr.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dorgl2.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dorglq.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dorgqr.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dorm2r.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dormbr.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dorml2.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dormlq.o)\
		$(LIB_LAPACK)(lib_standard/lapack/dormqr.o)\
		$(LIB_LAPACK)(lib_standard/lapack/ieeeck.o)\
		$(LIB_LAPACK)(lib_standard/lapack/ilaenv.o)
# Objects for sequential mode only.
OBJS_SEQ     = $(LIB_SEQ)(main/func_mpi.o)

# Objects for MPI mode only.
OBJS_MPI     = $(LIB_MPI)(main/func_mpi.o)\
               $(LIB_MPI)(main/IO_stop.o)\
	       $(LIB_MPI)(main/start_stop_mpi.o)


###################################
# Targets, Dependencies and Rules #
###################################

default: $(OBJS) $(OBJS_MATH) $(OBJS_SEQ)

all: seq mpi

libs: $(LIB_SEQ) $(LIB_COMMON)

seq: $(LIB_SEQ) $(LIB_COMMON) $(LIB_LINPACK) $(LIB_LAPACK) $(LIB_BLAS)
	cd tests ; $(MAKE) test_seq

mpi: $(LIB_MPI) $(LIB_COMMON) $(LIB_LINPACK) $(LIB_LAPACK) $(LIB_BLAS)
	cd tests ; $(MAKE) test_mpi

#$(LIB_COMMON): $(OBJS) $(OBJS_MATH) $(OBJS_LINPACK) $(OBJS_LAPACK) $(OBJS_BLAS)
$(LIB_COMMON): $(OBJS) $(OBJS_MATH)

$(LIB_SEQ): $(OBJS_SEQ)

$(LIB_MPI): $(OBJS_MPI)

$(LIB_COMMON)(main/anneal.o):             $(INC)
$(LIB_COMMON)(main/anneal_2.o):           $(INC)
$(LIB_COMMON)(main/derivs.o):             $(INC)
$(LIB_COMMON)(math/svd_gaus_test.o):      $(INC)
$(LIB_COMMON)(lib_standard/blas/dnrm2.o): $(INC)

$(LIB_SEQ)(main/func_mpi.o):              $(INC)

$(LIB_MPI)(main/func_mpi.o):              $(INC_MPI)
$(LIB_MPI)(main/IO_stop.o):               $(INC_MPI)
$(LIB_MPI)(main/start_stop_mpi.o):        $(INC_MPI)

# MPI files need special compilation, so we have to write special rules
# for them.  Ordinary (sequential) code is compiled using the default
# rules defined further down.
$(LIB_MPI)(main/func_mpi.o) : main/func_mpi.F
	echo using MPI compiler $(F77_MPI)
	$(F77_MPI) $(FFLAGS_MPI) -Iinclude -c $< -o $*.o
	$(AR) $(ARFLAGS) $@ $*.o
#	-rm $*.o
$(LIB_MPI)(main/IO_stop.o): main/IO_stop.F
	$(F77_MPI) $(FFLAGS_MPI) -Iinclude -c $< -o $*.o
	$(AR) $(ARFLAGS) $@ $*.o
#	-rm $*.o
$(LIB_MPI)(main/start_stop_mpi.o): main/start_stop_mpi.F
	$(F77_MPI) $(FFLAGS_MPI) -Iinclude -c $< -o $*.o
	$(AR) $(ARFLAGS) $@ $*.o
#	-rm $*.o

# Define suffix rule for .F files, since not all versions of make have
# one built in.  Also redefine .f rule, just to be sure.  Note that .x.a
# rules are special.  Instead of specifying how to make an archive file
# archive.a from a source file called archive.x, they say how to make an
# archive *member* archive.a(member.o) from a source file called
# member.x.  In writing .x.a rules, the make variables $< (the first
# prerequisite) and $* (the stem of $< with the matching suffix stripped
# off) refer to the member file, not the archive, while $@ is the name
# of the archive itself.
.SUFFIXES: .F .f .for
.F.a:
	echo using sequential compiler $(F77)
	$(F77) $(FFLAGS) -Iinclude -c $< -o $*.o
	$(AR) $(ARFLAGS) $@ $*.o
#	-rm $*.o
.f.a:
	echo using sequential compiler $(F77)
	$(F77) $(FFLAGS) -Iinclude -c $< -o $*.o
	$(AR) $(ARFLAGS) $@ $*.o
#	-rm $*.o
.for.a:
	echo using sequential compiler $(F77)
	$(F77) $(FFLAGS_NOOPT) -Iinclude -c $< -o $*.o
	$(AR) $(ARFLAGS) $@ $*.o
#	-rm $*.o

# PHONY is a Make keyword that prevents it from getting confused if there
# happens to be a file named "clean" in the directory.
# We used to remove include/parameters.h as part of clean_all because it can be machine dependent
# but now I just leave it as parameters/parameters.Linux.
.PHONY: clean clean_all
clean:
	-rm -f *.o main/*.o lib_standard/blas/*.o lib_standard/linpack/*.o math/*.o tests/*.o
clean_all:
	-rm -f *.o main/*.o lib_standard/blas/*.o lib_standard/linpack/*.o math/*.o tests/*.o \
            lib/libquench.a lib/libquench_seq.a lib/libquench_mpi.a lib/libblas.a lib/liblinpack.a lib/liblapack.a \
            tests/test_seq tests/test_mpi
