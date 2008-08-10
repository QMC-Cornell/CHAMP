# This Makefile is to be used with gmake.  On Linux systems, gmake is
# the default make.

system := $(shell uname)
hostname := $(shell uname -n)
# planck_compiler=g77
planck_compiler=ifort
export FC FFLAGS

# set default
FC=g77
FFLAGS = -O3 -malign-double -ffixed-line-length-none -fno-automatic

ifeq ($(findstring IRIX,$(system)),IRIX)
  FC = f77
  FFLAGS = -O2 -n32 -mips4 -DIRIX
# FFLAGS = -g -C -n32 -mips4 -DIRIX
  FC_MPI = mpif77
  FFLAGS_MPI = -O2 -n32 -mips4 -DIRIX -DMPI -DNON_STANDARD
# FFLAGS_MPI = -g -C -n32 -mips4 -DIRIX -DMPI -DNON_STANDARD
endif

ifeq ($(system),AIX)
  FC = xlf
  FFLAGS = -O2 -qnosave -qsigtrap -WF,-DAIX
  FC_MPI = mpxlf
  FFLAGS_MPI = -O2 -qnosave -qsigtrap -WF,-DAIX,-DMPI,-DNON_STANDARD
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
  FC = g77-3.0
# FC = g77
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

#gauss in Arias cluster
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
  FC_MPI = mpif90 -fc=ifort
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
#ifeq ($(hostname),oscbw01)
#ifeq ($(hostname),amd-login1)
## AMD Cluster at OSC using Portland Group Compiler
#  FC = pgf90
## FFLAGS = -fast -r8 -Mextend -Mdalign -Msave -Mbounds
#  FFLAGS = -fast -r8 -Mextend -Mdalign -Msave
#  FC_MPI = mpif90
## FFLAGS_MPI = -fast -r8 -Mextend -Mdalign -Msave -Mbounds
#  FFLAGS_MPI = -fast -r8 -Mextend -Mdalign -Msave
#endif
ifeq ($(hostname),amd-login1)
# AMD Cluster at OSC using Intel Compiler
 FC = ifc
# FFLAGS = -g -C -d4 -FI -extend_source -O0 -w -r8 -fpp2
# FFLAGS = -FI -extend_source -O2 -w -r8 -fpp2
  LDFLAGS = -Vaxlib
  FC_MPI = mpif77 -fc=ifc
  FFLAGS_MPI = -g -C -FI -extend_source -O2 -w -r8
# FFLAGS_MPI = -FI -extend_source -O2 -w -r8
  LDFLAGS_MPI = -Vaxlib
endif
ifeq ($(hostname),ia64)
# Intel Itanium cluster at OSC using MPI and Intel Compiler efc
# -FI -extend_source or -FI -4L132 do fixed format with 132 columns
# LIBHOME = $(QMCHOME)
  FC = efc
  FFLAGS = -FI -extend_source -O2 -w -r8 -Dderf=INTELERF -Dderfc=INTELERFC -fpp2
  FC_MPI = mpif77
  FFLAGS_MPI = -FI -extend_source -O2 -w -r8
# BLAS = $(MKL) -lderf -Vaxlib
  BLAS = $(MKL) -L$(HOME)/QMC/qmc/lib -laux -Vaxlib
endif
ifeq ($(hostname),mck-login1)
# Intel Itanium cluster at OSC using MPI and Intel Compiler efc
# LIBHOME = $(QMCHOME)
  FC = efc
  FFLAGS = -FI -extend_source -O2 -w -r8
  FC_MPI = mpif77
  FFLAGS_MPI = -FI -extend_source -O2 -w -r8
  BLAS = $(MKL) -L$(HOME)/QMC/qmc/lib -laux -Vaxlib
endif

# Argonne Blue-Gene
ifeq ($(hostname),login1)

# for debug front-end
  FC = xlf
  FFLAGS = -C -g -qlinedebug -qsigtrap -qfixed=132

# for nodes
  FC = blrts_xlf
  FFLAGS = -O5 -qarch=440 -qtune=440 -qfixed=132
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
ifeq ($(hostname),do)
  URI_CLUSTER = true
endif
ifeq ($(hostname),marx)
  URI_CLUSTER = true
endif
# -ax{i|M|K|W} i=Pentium Pro, M-Pentium MMX, K=Pentium3, W=Pentium4,Xeon
ifeq ($(URI_CLUSTER),true)
  FC = ifc
  FFLAGS = -c -extend_source -fpp -O2 -w -pad -axi
# -CA -CB -CU -g
  FFLAGS_NOOPT = -c -extend_source -fpp -O0 -w -pad -axi
  FC_MPI = mpif77
  FLINKER_MPI = mpif77
  FFLAGS_MPI = -DMPI -DNON_STANDARD -DAIX -fpp -w
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
ifeq ($(system),CYGWIN_NT-5.0)
  FC = win32fe ifl
  FFLAGS = -FI -4L132 -Qlowercase -w -O2
# FFLAGS = -Qlowercase -w -O2 -4L132 /Qoption,link,/verbose
# BLAS = $(MKL) -lderf -Vaxlib
  LDFLAGS = -L/cygdrive/c/progra~1/mifd68~1/Lib /Qoption,link,/verbose mkl_c.lib mkl_p3.lib
  FC_MPI = win32fe ifl
# MPIPro164 for v2.  For v1,vplus use MPIPro and change library names
# Another way to deal with spaces in names is "Program~s" instead of "Program\ Files"
  FFLAGS_MPI = -FI -4L132 -Qlowercase -w -O2 -I/cygdrive/c/Program\ Files/MPIPro164/include
  LDFLAGS_MPI = -L/cygdrive/c/Program\ Files/MPIPro164/lib MPIPro.lib MPIPro_cdec.lib
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

# OBJS = my_second.o ranf.o shellnn.o intpol.o intpol_monot.o grater.o \
# hartre.o deriv.o deriv2.o divergence.o lapla3.o laplac.o fderiv.o simson.o bode.o \
# hartre_grid2.o deriv_grid2.o deriv2_grid2.o divergence_grid2.o laplac_grid2.o fderiv_grid2.o \
# bode_grid2.o

OBJS = my_second.o intpol.o simson.o

libcyrus.a: $(OBJS)
	ar r $@ $(OBJS)
	ranlib $@

clean:
	-rm *.o

clean_all: clean
	-rm *.a
