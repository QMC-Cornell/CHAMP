##########################################################################
# MACHINE-DEPENDENT COMPILATION OPTIONS FOR CHAMP
##########################################################################

system := $(shell uname)
hostname := $(shell uname -n)
uname := $(shell uname)
machine := $(shell uname -m)

####################### DEFAULTS ########################################

# Libraries --------------------------------------------------------------

LIBHOME=../lib
# LIBHOME=/home/cyrus/champ/lib

QUENCH=-L$(LIBHOME)/SimulatedAnnealing/quench_anneal/lib -lquench -lquench_seq
CYRUSLIB=-L$(LIBHOME)/lib -lcyrus

LINPACK=-L$(LIBHOME)/lib2/linpack -llinpack
#LINPACK=-L/usr/lib -llinpack

LAPACK=-L$(LIBHOME)/lib2/lapack -llapack
#LAPACK=-L/usr/lib/ -llapack

BLAS=-L$(LIBHOME)/lib2/blas -lblas
#BLAS=-L/usr/lib/ -lblas

PSPLINE=-L$(LIBHOME)/lib2/pspline/pspline -lpspline

EINSPLINE=-L$(LIBHOME)/lib2/einspline/lib -leinspline

LIBS=$(CYRUSLIB) $(LINPACK) $(LAPACK) $(BLAS) $(QUENCH) $(LINPACK) $(BLAS) $(PSPLINE) $(EINSPLINE) $(QUENCH)

LIBS_MAKE = $(LIBHOME)/lib/libcyrus.a $(LIBHOME)/lib2/blas/libblas.a $(LIBHOME)/lib2/lapack/liblapack.a $(LIBHOME)/lib2/linpack/liblinpack.a $(LIBHOME)/lib2/einspline/lib/libeinspline.a $(LIBHOME)/lib2/pspline/pspline/libpspline.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench_seq.a

# C compilation ------------------------------------------------------------

CC=icc
CXX=icpc
CC_MPI=icc

CC_FAST_FLAGS=-O3 -DUNDERSCORE
CC_PROF_FLAGS=-O3 -DUNDERSCORE -pg
CC_DEBUG_FLAGS=-DUNDERSCORE

CC_FLAGS=$(CC_FAST_FLAGS)
CC_MPI_FLAGS=$(CC_FLAGS)

# Fortran compilation ------------------------------------------------------

# These are defaults using the ifort compiler and are used e.g. on the old AMD cluster dft, and on nanolab
F95=ifort
F95_MPI=mpif77

# F95_FAST_FLAGS= -zero -r8 -pad -align all -fpp -ip -O3
F95_FAST_FLAGS= -zero -r8 -pad -align all -O3 -ax,SSE3,SSE4.2 -ip -diag-disable vec -fpp -shared-intel -DNOEINSPLINE
F95_PROF_FLAGS= -zero -r8 -pad -align all -fpp -pg -O3
F95_DEBUG_FLAGS=-zero -r8 -pad -align all -fpp -ip -O0 -check all -check nooutput_conversion -debug extended -debug-parameters all -g -traceback -error_limit 1 -DDEBUG

# by default, compilation for f77 is the same than for f95
F77=$(F95)
F77_MPI=$(F95_MPI)
F77_FAST_FLAGS=$(F95_FAST_FLAGS) -extend_source
F77_DEBUG_FLAGS=$(F95_DEBUG_FLAGS) -extend_source
F77_PROF_FLAGS=$(F95_PROF_FLAGS) -extend_source

# by default, fast options
F95_FLAGS=$(F95_FAST_FLAGS)
F77_FLAGS=$(F77_FAST_FLAGS) 

# by default, same options for MPI
F95_MPI_FLAGS=$(F95_FLAGS) -DMPI
F77_MPI_FLAGS=$(F77_FLAGS) -DMPI
F95_DEBUG_MPI_FLAGS=$(F95_DEBUG_FLAGS) -DMPI
F77_DEBUG_MPI_FLAGS=$(F77_DEBUG_FLAGS) -DMPI

# by default, include and module directories options
INCCMD=-I
MODCMD=-J

# Linking  ------------------------------------------------------------------

# by default, linker is f95 compiler
LD=$(F95)
LD_MPI=$(F95_MPI)

LD_FAST_FLAGS=$(F95_FLAGS)
LD_PROF_FLAGS=-O3 -pg
# LD_DEBUG_FLAGS=-align all
LD_DEBUG_FLAGS=$(LD_FAST_FLAGS)

# by default, fast options
LD_FLAGS=$(LD_FAST_FLAGS)

LD_END=
#LD_END_MPI=-lmpi -lpthread
LD_END_MPI=

# Executables -------------------------------------------------------------

EXE=champ.exe
MPI_EXE=champ_mpi.exe
DEBUG_EXE=champ.exe
PROF_EXE=champ_prof.exe

# Export these for compiling the libraries
export F77 F77_FLAGS LD_FLAGS CC CXX EINSPLINE USE_GCC 

####################### MACHINES ########################################
ifeq ($(machine),x86_64)
ifeq ($(uname),Linux)
     CC=icc
     CXX=icpc
     CC_MPI=mpicc
     F95=ifort
     F95_MPI=mpif90
#    F95_FAST_FLAGS=-zero -r8 -pad -align all -O3 -xST -ip -diag-disable vec -fpp -DNOEINSPLINE -shared-intel
     F95_FAST_FLAGS=-zero -r8 -pad -align all -O3 -ax,SSE3,SSE4.2 -ip -diag-disable vec -fpp -shared-intel -DNOEINSPLINE
     F95_DEBUG_FLAGS=-zero -r8 -pad -align all -O0 -fpp -ip -check all -check nooutput_conversion -debug extended -debug-parameters all -g -traceback -error_limit 1 -DDEBUG -DNOEINSPLINE -shared-intel
     EINSPLINE=
     LIBS_MAKE= $(LIBHOME)/lib/libcyrus.a $(LIBHOME)/lib2/blas/libblas.a $(LIBHOME)/lib2/lapack/liblapack.a $(LIBHOME)/lib2/linpack/liblinpack.a $(LIBHOME)/lib2/pspline/pspline/libpspline.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench_seq.a
endif
endif


# New dft cluster at Cornell
# We used to use -w to suppress warnings because otherwise the compilation of parser_tools fails after priting 10 error msgs. such as
#  Error: parser_tools_mod.f90, line 581: The dimensions of this array have been defined more than once.   [VALUE_LIST]
#    character(len=*), allocatable, intent(out) :: value_list (:)
# -------------------------------------------------^
# For some reason these error msgs. have disappeared, so remove the -w
ifeq ($(hostname),gauss)
  CC_MPI=mpicc
  F95_MPI=mpif90
  CC_FAST_FLAGS=-O3 -DUNDERSCORE -mtune=core2 -march=core2 -xT 
# F95_FAST_FLAGS=-zero -w -r8 -pad -align all -fpp -ip -O3 -mtune=core2 -march=core2 -xT -DNOEINSPLINE
  F95_FAST_FLAGS=-zero -r8 -pad -align all -fpp -ip -O3 -mtune=core2 -march=core2 -xT -DNOEINSPLINE
  F77_FAST_FLAGS=$(F95_FAST_FLAGS) -extend_source 
  F95_DEBUG_FLAGS=$(F95_FAST_FLAGS) -check all -check nooutput_conversion -debug extended -debug-parameters all -g -traceback -error_limit 1 -DDEBUG
  F77_DEBUG_FLAGS=$(F77_FAST_FLAGS) -check all -check nooutput_conversion -debug extended -debug-parameters all -g -traceback -error_limit 1 -DDEBUG
  LD_FAST_FLAGS=$(F95_FAST_FLAGS)
  EINSPLINE=
  LIBS_MAKE = $(LIBHOME)/lib/libcyrus.a $(LIBHOME)/lib2/blas/libblas.a $(LIBHOME)/lib2/lapack/liblapack.a $(LIBHOME)/lib2/linpack/liblinpack.a $(LIBHOME)/lib2/pspline/pspline/libpspline.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench_seq.a
endif

#new nodes at cornell
ifeq ($(findstring node,$(hostname)), node)
  CC=gcc
  CXX=g++
  CC_MPI=mpicc
  F95=gfortran
  F95_MPI=mpif90
  F95_FAST_FLAGS=-finit-local-zero -O3 -fdefault-real-8 -ffree-line-length-none -x f95-cpp-input -Wall  -DNOEINSPLINE
  F77_FAST_FLAGS=-finit-local-zero -O3 -fdefault-real-8 -ffixed-line-length-132 -x f77-cpp-input -Wall  -DNOEINSPLINE
  LD=gfortran
  LD_FLAGS=-O3
  LD_MPI=mpif90
  EINSPLINE=
endif

ifeq ($(findstring login,$(hostname)), login)
  CC=gcc
  CXX=g++
  CC_MPI=mpicc
  F95=gfortran
  F95_MPI=mpif90
# F95_FAST_FLAGS=-finit-local-zero -O3 -fdefault-real-8 -ffree-line-length-none -x f95-cpp-input -Wall  -DNOEINSPLINE
# F77_FAST_FLAGS=-finit-local-zero -O3 -fdefault-real-8 -ffixed-line-length-132 -x f77-cpp-input -Wall  -DNOEINSPLINE
  F95_FAST_FLAGS=-finit-local-zero -O3 -ffree-line-length-none -x f95-cpp-input -Wall  -DNOEINSPLINE
  F77_FAST_FLAGS=-finit-local-zero -O3 -ffixed-line-length-132 -x f77-cpp-input -Wall  -DNOEINSPLINE
  LD=gfortran
  LD_FLAGS=-O3
  LD_MPI=mpif90
  EINSPLINE=
endif

# nanolab at Cornell
ifeq ($(findstring cnf.cornell,$(hostname)),cnf.cornell)
  LD_FAST_FLAGS=$(F95_FLAGS)
  LD_PROF_FLAGS=-O3 -pg
  LD_DEBUG_FLAGS=-align all
  LD_IPO=-O3 -ipo -align all -fast
endif

#qmc1 at cornell
ifeq ($(hostname),qmc1)
  F95=ifort
  CC=gcc
  CXX=g++
  CC_MPI=mpicc

  F95_FAST_FLAGS= -zero -r8 -pad -align all -O3 -fpp -ip
  F95_PROF_FLAGS= -zero -r8 -pad -align all -O3 -fpp -pg
  F95_DEBUG_FLAGS=-zero -r8 -pad -align all -O0 -fpp -ip -check all -check nooutput_conversion -debug extended -debug-parameters all -g -traceback -error_limit 1 -DDEBUG

# by default, compilation for f77 is the same than for f95
  F77=$(F95)
  F77_FAST_FLAGS=$(F95_FAST_FLAGS) -extend_source
  F77_DEBUG_FLAGS=$(F95_DEBUG_FLAGS) -extend_source
  F77_PROF_FLAGS=$(F95_PROF_FLAGS) -extend_source

# by default, fast options
  F95_FLAGS=$(F95_FAST_FLAGS)
  F77_FLAGS=$(F77_FAST_FLAGS)

  USE_GCC=yes
endif

#qmc2 at cornell
ifeq ($(hostname),qmc2)
  CC=gcc
  CXX=g++
  CC_MPI=mpicc
  F95=gfortran
  F95_MPI=mpif90
  F95_FAST_FLAGS=-finit-local-zero -O3 -fdefault-real-8 -ffree-line-length-none -x f95-cpp-input -Wall -DNOEINSPLINE
  F77_FAST_FLAGS=-finit-local-zero -O3 -fdefault-real-8 -ffixed-line-length-132 -x f77-cpp-input -Wall -DNOEINSPLINE
  F95_DEBUG_FLAGS=$(F95_FAST_FLAGS) -g -fbounds-check
  F77_DEBUG_FLAGS=$(F77_FAST_FLAGS) -g -fbounds-check
  F95_PROF_FLAGS=$(F95_FAST_FLAGS) -pg
  F77_PROF_FLAGS=$(F77_FAST_FLAGS) -pg
  LD_FLAGS=-O3
  LD_DEBUG_FLAGS=-Wall -g -fbounds-check
  LD_PROF_FLAGS=-pg
  EINSPLINE=
  LIBS_MAKE= $(LIBHOME)/lib/libcyrus.a $(LIBHOME)/lib2/blas/libblas.a $(LIBHOME)/lib2/lapack/liblapack.a $(LIBHOME)/lib2/linpack/liblinpack.a $(LIBHOME)/lib2/pspline/pspline/libpspline.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench_seq.a
endif

#All qmc at cornell
ifeq ($(findstring qmc,$(hostname)), qmc)
  CC=gcc
  CXX=g++
  CC_MPI=mpicc
  F95=gfortran
  F95_MPI=mpif90
# F95_FAST_FLAGS=-finit-local-zero -O3 -fdefault-real-8 -ffree-line-length-none -x f95-cpp-input -Wall -DNOEINSPLINE
  F95_FAST_FLAGS=-finit-local-zero -O3 -ffree-line-length-none -x f95-cpp-input -Wall -DNOEINSPLINE
# F77_FAST_FLAGS=-finit-local-zero -O3 -fdefault-real-8 -ffixed-line-length-132 -x f77-cpp-input -Wall -DNOEINSPLINE
  F77_FAST_FLAGS=-finit-local-zero -O3 -ffixed-line-length-132 -x f77-cpp-input -Wall -DNOEINSPLINE
  F95_DEBUG_FLAGS=$(F95_FAST_FLAGS) -g -fbounds-check
  F77_DEBUG_FLAGS=$(F77_FAST_FLAGS) -g -fbounds-check
  F95_PROF_FLAGS=$(F95_FAST_FLAGS) -pg
  F77_PROF_FLAGS=$(F77_FAST_FLAGS) -pg
  LD_FLAGS=-O3
  LD_DEBUG_FLAGS=-Wall -g -fbounds-check
  LD_PROF_FLAGS=-pg
  EINSPLINE=
  LIBS_MAKE= $(LIBHOME)/lib/libcyrus.a $(LIBHOME)/lib2/blas/libblas.a $(LIBHOME)/lib2/lapack/liblapack.a $(LIBHOME)/lib2/linpack/liblinpack.a $(LIBHOME)/lib2/pspline/pspline/libpspline.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench_seq.a
endif

#qmc8 at cornell
ifeq ($(hostname),qmc8)
  CC=gcc
  CXX=g++
  CC_MPI=mpicc
  F95=gfortran
  F95_MPI=mpif90
# F95_FAST_FLAGS=-finit-local-zero -O3 -fdefault-real-8 -ffree-line-length-none -x f95-cpp-input -Wall -DNOEINSPLINE
  F95_FAST_FLAGS=-finit-local-zero -O3 -ffree-line-length-none -x f95-cpp-input -Wall -DNOEINSPLINE -fbounds-check
# F77_FAST_FLAGS=-finit-local-zero -O3 -fdefault-real-8 -ffixed-line-length-132 -x f77-cpp-input -Wall -DNOEINSPLINE
  F77_FAST_FLAGS=-finit-local-zero -O3 -ffixed-line-length-132 -x f77-cpp-input -Wall -DNOEINSPLINE -fbounds-check
  F95_DEBUG_FLAGS=$(F95_FAST_FLAGS) -g -fbounds-check
  F77_DEBUG_FLAGS=$(F77_FAST_FLAGS) -g -fbounds-check
  F95_PROF_FLAGS=$(F95_FAST_FLAGS) -pg
  F77_PROF_FLAGS=$(F77_FAST_FLAGS) -pg
  LD_FLAGS=-O3
  LD_DEBUG_FLAGS=-Wall -g -fbounds-check
  LD_PROF_FLAGS=-pg
  EINSPLINE=
  LIBS_MAKE= $(LIBHOME)/lib/libcyrus.a $(LIBHOME)/lib2/blas/libblas.a $(LIBHOME)/lib2/lapack/liblapack.a $(LIBHOME)/lib2/linpack/liblinpack.a $(LIBHOME)/lib2/pspline/pspline/libpspline.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench_seq.a
endif

#jasmine at cornell
ifeq ($(hostname),jasmine)
  CC=gcc
  CXX=g++
  CC_MPI=mpicc
  F95=gfortran
  F95_MPI=mpif90
  F95_FAST_FLAGS=-finit-local-zero -O3 -fdefault-real-8 -ffree-line-length-none -x f95-cpp-input -Wall -fbounds-check -DNOEINSPLINE
  F77_FAST_FLAGS=-finit-local-zero -O3 -fdefault-real-8 -ffixed-line-length-132 -x f77-cpp-input -Wall -fbounds-check -DNOEINSPLINE
  F95_SLOW_FLAGS=-finit-local-zero -O0 -fdefault-real-8 -ffree-line-length-none -x f95-cpp-input -Wall -fbounds-check -DNOEINSPLINE
  F77_SLOW_FLAGS=-finit-local-zero -O0 -fdefault-real-8 -ffixed-line-length-132 -x f77-cpp-input -Wall -fbounds-check -DNOEINSPLINE
  F95_DEBUG_FLAGS=$(F95_SLOW_FLAGS) -g -fbounds-check
  F77_DEBUG_FLAGS=$(F77_SLOW_FLAGS) -g -fbounds-check
  F95_PROF_FLAGS=$(F95_FAST_FLAGS) -pg
  F77_PROF_FLAGS=$(F77_FAST_FLAGS) -pg
  LD=gfortran
  LD_FLAGS=-O3
  LD_MPI=mpif90
  LD_DEBUG_FLAGS=-Wall -g -fbounds-check
  LD_PROF_FLAGS=-pg
  EINSPLINE=
  LIBS_MAKE= $(LIBHOME)/lib/libcyrus.a $(LIBHOME)/lib2/blas/libblas.a $(LIBHOME)/lib2/lapack/liblapack.a $(LIBHOME)/lib2/linpack/liblinpack.a $(LIBHOME)/lib2/pspline/pspline/libpspline.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench_seq.a
endif

# Theory Center's linux cluster
# ifeq ($(hostname),linuxlogin1.tc.cornell.edu)
ifeq ($(findstring linuxlogin,$(hostname)), linuxlogin)
 F95_MPI=mpif90
 LD_MPI=mpif90
endif

# clastos in Paris
# Remember to put . /usr/local/intel/Compiler/11.1/064/bin/intel64/ifortvars_intel64.sh in .bashrc 
ifeq ($(hostname),clastos.lct.jussieu.fr)

# CC=/usr/local/bin/mpicc_64
# CXX=/usr/local/bin/mpicxx_64
# F95=/usr/local/bin/mpif90_64
# F95_FAST_FLAGS=-zero -r8 -pad -align all -xW -ip -diag-disable vec -fpp -DNOEINSPLINE
# F95_DEBUG_FLAGS=-zero -r8 -pad -align all -O0 -fpp -ip -check all -check nooutput_conversion -debug extended -debug-parameters all -g -traceback -error_limit 1 -DDEBUG -DNOEINSPLINE
# EINSPLINE=
# F95_MPI=/usr/local/bin/mpif90_64
# CC_MPI=/usr/local/bin/mpicc_64
# LIBS_MAKE= $(LIBHOME)/lib/libcyrus.a $(LIBHOME)/lib2/blas/libblas.a $(LIBHOME)/lib2/lapack/liblapack.a $(LIBHOME)/lib2/linpack/liblinpack.a $(LIBHOME)/lib2/pspline/pspline/libpspline.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench_seq.a

#  export OMPI_F77=ifort
#  export OMPI_FC=ifort
  CC=icc
  CXX=icc
  CC_MPI=/usr/local/openmpi/bin/mpicc
  F95=ifort
  F95_MPI=/usr/local/openmpi/bin/mpif90
#  F95_FAST_FLAGS=-zero -r8 -pad -align all -O3 -xT -ip -diag-disable vec -fpp -DNOEINSPLINE -shared-intel
  F95_FAST_FLAGS=-zero -r8 -pad -align all -O3 -ip -diag-disable vec -fpp -DNOEINSPLINE -shared-intel
  F95_DEBUG_FLAGS=-zero -r8 -pad -align all -O0 -fpp -ip -check all -check nooutput_conversion -debug extended -debug-parameters all -g -traceback -error_limit 1 -DDEBUG -DNOEINSPLINE -shared-intel
  EINSPLINE=
  LIBS_MAKE= $(LIBHOME)/lib/libcyrus.a $(LIBHOME)/lib2/blas/libblas.a $(LIBHOME)/lib2/lapack/liblapack.a $(LIBHOME)/lib2/linpack/liblinpack.a $(LIBHOME)/lib2/pspline/pspline/libpspline.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench_seq.a

endif

# deuterium in Paris
ifeq ($(hostname),deuterium.lct.jussieu.fr)

# CC=/usr/local/bin/mpicc_64
# CXX=/usr/local/bin/mpicxx_64
# F95=/usr/local/bin/mpif90_64
# F95_FAST_FLAGS=-zero -r8 -pad -align all -O3 -xT -ip -diag-disable vec -fpp -DNOEINSPLINE
# F95_DEBUG_FLAGS=-zero -r8 -pad -align all -O0 -fpp -ip -check all -check nooutput_conversion -debug extended -debug-parameters all -g -traceback -error_limit 1 -DDEBUG -DNOEINSPLINE
# EINSPLINE=
# F95_MPI=/usr/local/bin/mpif90_64
# CC_MPI=/usr/local/bin/mpicc_64
# LIBS_MAKE= $(LIBHOME)/lib/libcyrus.a $(LIBHOME)/lib2/blas/libblas.a $(LIBHOME)/lib2/lapack/liblapack.a $(LIBHOME)/lib2/linpack/liblinpack.a $(LIBHOME)/lib2/pspline/pspline/libpspline.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench_seq.a

#  export OMPI_F77=ifort
#  export OMPI_FC=ifort
  CC=icc
  CXX=icc
  CC_MPI=/usr/local/openmpi_ifort11/bin/mpicc
  F95=ifort
  F95_MPI=/usr/local/openmpi_ifort11/bin/mpif90
  F95_FAST_FLAGS=-zero -r8 -pad -align all -O3 -xT -ip -diag-disable vec -fpp -DNOEINSPLINE -shared-intel
  F95_DEBUG_FLAGS=-zero -r8 -pad -align all -O0 -fpp -ip -check all -check nooutput_conversion -debug extended -debug-parameters all -g -traceback -error_limit 1 -DDEBUG -DNOEINSPLINE -shared-intel
  EINSPLINE=
  LIBS_MAKE= $(LIBHOME)/lib/libcyrus.a $(LIBHOME)/lib2/blas/libblas.a $(LIBHOME)/lib2/lapack/liblapack.a $(LIBHOME)/lib2/linpack/liblinpack.a $(LIBHOME)/lib2/pspline/pspline/libpspline.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench_seq.a

endif


# hydrogene in Paris
ifeq ($(hostname),hydrogene.lct.jussieu.fr)
# F95_FAST_FLAGS=-zero -r8 -pad -align all -O3 -ip -fpp -DNOEINSPLINE
# CC=gcc
# EINSPLINE=
# LIBS_MAKE= $(LIBHOME)/lib/libcyrus.a $(LIBHOME)/lib2/blas/libblas.a $(LIBHOME)/lib2/lapack/liblapack.a $(LIBHOME)/lib2/linpack/liblinpack.a $(LIBHOME)/lib2/pspline/pspline/libpspline.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench_seq.a
  CC=icc
  CXX=icc
  F95=ifort
  F95_FAST_FLAGS=-zero -r8 -pad -align all -O3 -ip -diag-disable vec -fpp -DNOEINSPLINE -shared-intel
  F95_DEBUG_FLAGS=-zero -r8 -pad -align all -O0 -fpp -ip -check all -check nooutput_conversion -debug extended -debug-parameters all -g -traceback -error_limit 1 -DDEBUG -DNOEINSPLINE -shared-intel
  EINSPLINE=
  LIBS_MAKE= $(LIBHOME)/lib/libcyrus.a $(LIBHOME)/lib2/blas/libblas.a $(LIBHOME)/lib2/lapack/liblapack.a $(LIBHOME)/lib2/linpack/liblinpack.a $(LIBHOME)/lib2/pspline/pspline/libpspline.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench_seq.a
endif

# NASA Ames Columbia Itanium cluster
ifeq ($(hostname),cfe2)
 CC=icc
 CXX=icc
 F95_MPI=ifort
 F77_MPI=ifort
 LD_MPI=ifort
 F95_FAST_FLAGS=-O3 -ip  -align all -fpp -DNOEINSPLINE
 LD_FAST_FLAGS=-O3 -ip -align all -fpp -DNOEINSPLINE
 F95_DEBUG_FLAGS=-zero -r8 -pad -align all -O0 -fpp -ip -check all -check nooutput_conversion -debug extended -debug-parameters all -g -traceback -error_limit 1 -DDEBUG -DNOEINSPLINE
 LD_END_MPI=-nothread -lmpi -lpthread
 EINSPLINE=
 LIBS_MAKE= $(LIBHOME)/lib/libcyrus.a $(LIBHOME)/lib2/blas/libblas.a $(LIBHOME)/lib2/lapack/liblapack.a $(LIBHOME)/lib2/linpack/liblinpack.a $(LIBHOME)/lib2/pspline/pspline/libpspline.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench_seq.a
endif

# Pleiades at NASA
ifeq ($(hostname),p4fe1)
  CC=icc
  CXX=icc
  F95_MPI=mpiifort
  F77_MPI=mpiifort
  LD_MPI=mpiifort
#  F95_MPI=ifort
#  F77_MPI=ifort
#  LD_MPI=ifort
#  F95_MPI=mpif90
#  F77_MPI=mpif90
#  LD_MPI=mpif90
  EINSPLINE=
  LIBS_MAKE= $(LIBHOME)/lib/libcyrus.a
  F95_FAST_FLAGS=-O3 -ip -align all -fpp -DNOEINSPLINE -mcmodel medium
-shared-intel
  LD_FAST_FLAGS=-O3 -ip -align all -fpp -DNOEINSPLINE -mcmodel medium
-shared-intel
#  F95_FAST_FLAGS=-O3 -ip  -fpe0 -traceback -align all -fpp -DNOEINSPLINE
#  LD_FAST_FLAGS=-O3 -ip -fpe0 -traceback -align all -fpp -DNOEINSPLINE
  LD_END_MPI=-nothread -lmpi -lpthread
endif

# NNIN cluster at Stanford
ifeq ($(hostname),exiton)
 CC=icc
 CXX=icc
 F95=ifort
 CC_MPI=mpicc
 F95_MPI=mpif90
 F95_FAST_FLAGS=-O3 -ip -align all -fpp -DNOEINSPLINE
 F95_DEBUG_FLAGS=$(F95_FAST_FLAGS) -check all -check nooutput_conversion -debug extended -debug-parameters all -g -traceback -error_limit 1 -DDEBUG
 EINSPLINE=
endif

ifeq ($(hostname),cpstar)
 CC=icc
 CXX=icc
 F95=ifort
 CC_MPI=mpicc
 F95_MPI=mpif90
 F95_FAST_FLAGS=-O3 -ip -align all -fpp -DNOEINSPLINE
 F95_DEBUG_FLAGS=$(F95_FAST_FLAGS) -check all -check nooutput_conversion -debug extended -debug-parameters all -g -traceback -error_limit 1 -DDEBUG
 EINSPLINE=
endif

# OSC Itanium cluster, requires Intel compiler >9.0
ifeq ($(hostname),ipf-login1)
  F95_MPI=mpif90
  LAPACK=-L/usr/local/mkl-7.0.1.6/mkl701/lib/64 -lmkl_lapack -lmkl_ipf -lguide -lpthread
  BLAS=-L/usr/local/mkl-7.0.1.6/mkl701/lib/64 -lmkl_ipf -lguide -lpthread
endif

# OSC Pentium IV cluster, requires Intel compiler >9.0
ifeq ($(hostname),piv-login1)
  F95_MPI=mpif90
  LAPACK=-L/usr/local/mkl-7.0.17/mkl70/lib/32 -lmkl_lapack -lmkl_ia32 -lguide -lpthread -Vaxlib
  BLAS=-L/usr/local/mkl-7.0.17/mkl70/lib/32 -lmkl_ia32 -lguide -lpthread
endif

# OSC Opteron cluster
ifeq ($(hostname),opt-login01.osc.edu)
# Serial (Intel compiler), requires ACML for Intel module
# module load acml-intel

# Parallel (Intel compiler), requires ACML & MPI for Intel
# module load acml-intel
# module swap mpi mvapich-0.9.8-intel

  F95_MPI = mpif77 -fc=ifort
  CC_MPI = mpicc -cc=icc
  LD_MPI = mpif77 -fc=ifort
  LD_FAST_FLAGS=-O3 -align all

# Serial (Portland compiler), requires ACML module
# module load acml
# F95 = pgf95
# CC = pgcc
# LD = pgf95
# F77_FAST_FLAGS = -Mpreprocess -Mextend -O2 -tp amd64e
# F95_FAST_FLAGS = -Mpreprocess -O2 -tp amd64 -DPGI
# CC_FAST_FLAGS = -DUNDERSCORE -O2 -tp amd64e
# LD_FAST_FLAGS = -O2 -tp amd64e

# Parallel (Portland compiler), requires ACML module
# module load acml
# F95_MPI = mpif90
# CC_MPI = mpicc
# LD_MPI = mpif90
# F77_MPI_FLAGS = -Mpreprocess -Mextend -DMPI -O2 -tp amd64e
# F95_MPI_FLAGS = -Mpreprocess -DMPI -O2 -tp amd64e -DPGI
# CC_MPI_FLAGS = -DMPI -DUNDERSCORE -O2 -tp amd64e

  LAPACK=${ACML}
  BLAS=${ACML}
endif

# Richard's Cornell machine
ifeq ($(hostname),kyra.ccmr.cornell.edu)
 LAPACK=-L/opt/intel/mkl_serial_8.1/lib/32/ -lmkl_lapack
 BLAS=-L/opt/intel/mkl_serial_8.1/lib/32/ -lmkl_ia32 -lguide -lpthread
endif

# Machines at OSU
ifeq ($(hostname),devon)
 LAPACK=-L${MKL_HOME}/lib/32 -lmkl_lapack
 BLAS=-L${MKL_HOME}/lib/32 -lmkl_ia32 -lguide -lpthread
 CC=icc
endif
ifeq ($(hostname),kunoy)
 LAPACK=-L${MKL_HOME}/lib/32 -lmkl_lapack
 BLAS=-L${MKL_HOME}/lib/32 -lmkl_ia32 -lguide -lpthread
endif
ifeq ($(hostname),banks)
 LAPACK=-L${MKL_HOME}/lib/32 -lmkl_lapack
 BLAS=-L${MKL_HOME}/lib/32 -lmkl_ia32 -lguide -lpthread
endif
ifeq ($(hostname),pico)
 LAPACK=-L${MKL_HOME}/lib/32 -lmkl_lapack
 BLAS=-L${MKL_HOME}/lib/32 -lmkl_ia32 -lguide -lpthread
endif
ifeq ($(hostname),vega)
 LAPACK=-L${MKL_HOME}/lib/32 -lmkl_lapack
 BLAS=-L${MKL_HOME}/lib/32 -lmkl_ia32 -lguide -lpthread
endif
ifeq ($(hostname),naxos)
 LAPACK=-L${MKL_HOME}/lib/32 -lmkl_lapack
 BLAS=-L${MKL_HOME}/lib/32 -lmkl_ia32 -lguide -lpthread
endif
ifeq ($(hostname),samos)
 LAPACK=-L${MKL_HOME}/mkl701/lib/64 -lmkl_lapack
 BLAS=-L${MKL_HOME}/mkl701/lib/64 -lmkl_ia64 -lguide -lpthread
endif

# Argonne Blue-Gene BG/L, login1
ifeq ($(hostname),login1a)

# serial
  F77 = xlf
  F77_FAST_FLAGS = -O5 -qarch=440 -qtune=440 -qsuffix=cpp=f -WF,-DUNDERSCORE -qfixed=132 -Ifit -Ivmc -Ivmc/MPI -Idmc -Idmc/dmc_elec -Idmc/dmc_elec/MPI -qsuppress=1512-039 -qsuppress=1512-047
  F95 = xlf95
  F95_FAST_FLAGS = -O5 -qarch=440 -qtune=440 -qsuffix=cpp=f90 -qsuppress=1512-039 -qsuppress=1512-047
  CC = xlc
  CC_FAST_FLAGS = -O5
  LD = xlf95
  LD_FAST_FLAGS =

# serial debug
  F77_DEBUG_FLAGS = -C -g -qsigtrap -qlinedebug -qsuffix=cpp=f -WF,-DDEBUG,-DUNDERSCORE -qfixed=132 -Ifit -Ivmc -Ivmc/MPI -Idmc -Idmc/dmc_elec -Idmc/dmc_elec/MPI -qsuppress=1512-039 -qsuppress=1512-047
  F95_DEBUG_FLAGS = -C -g -qsigtrap -qlinedebug -qsuffix=cpp=f90 -WF,-DDEBUG -qsuppress=1512-039 -qsuppress=1512-047
  CC_DEBUG_FLAGS =
  LD_DEBUG_FLAGS =

# MPI debug
  F95_MPI = mpif90.ibm
  F77_MPI_FLAGS = -C -g -qsigtrap -qlinedebug -qflttrap=enable:inexact:invalid:overflow:underflow:zerodivide -qarch=440 -qtune=440 -qsuffix=cpp=f -WF,-DMPI -qfixed=132 -Ifit -Ivmc -Ivmc/MPI -Idmc -Idmc/dmc_elec -Idmc/dmc_elec/MPI -qsuppress=1512-039 -qsuppress=1512-047
  F95_MPI_FLAGS = -C -g -qsigtrap -qlinedebug -qflttrap=enable:inexact:invalid:overflow:underflow:zerodivide -qarch=440 -qtune=440 -qsuffix=cpp=f90  -WF,-DMPI -qsuppress=1512-039 -qsuppress=1512-047
  CC_MPI = mpicc.ibm
  CC_MPI_FLAGS =  -C -g -qsigtrap -qlinedebug -qflttrap=enable:inexact:invalid:overflow:underflow:zerodivide
  LD_MPI = mpif90.ibm

# MPI
  F95_MPI = mpif90.ibm
  F77_MPI_FLAGS = -O5 -qarch=440 -qtune=440 -qsuffix=cpp=f -WF,-DMPI,-DUNDERSCORE -qfixed=132 -Ifit -Ivmc -Ivmc/MPI -Idmc -Idmc/dmc_elec -Idmc/dmc_elec/MPI -qsuppress=1512-039 -qsuppress=1512-047
  F95_MPI_FLAGS = -O5 -qarch=440 -qtune=440 -qsuffix=cpp=f90  -WF,-DMPI -qsuppress=1512-039 -qsuppress=1512-047
  CC_MPI = mpicc.ibm
  CC_MPI_FLAGS = -O5
  LD_MPI = mpif90.ibm
endif

# Jaguar at NCCS Oak Ridge
ifeq ($(findstring jaguar,$(hostname)), jaguar)

# Portland compiler (requires ACML module?)
## module load acml

# Serial compilation
# the following configuration compiles but the code segm. faults
# F95 = pgf95
# CC = pgcc
# LD = pgf95
# F77_FAST_FLAGS = -Mpreprocess -Mextend -O2 -tp amd64e
# F95_FAST_FLAGS = -Mpreprocess -O2 -tp amd64 -DPGI
# CC_FAST_FLAGS = -DUNDERSCORE -O2 -tp amd64e
# LD_FAST_FLAGS = -O2 -tp amd64e

# MPI compilation
# F95_MPI = mpif90
# CC_MPI = mpicc
# LD_MPI = mpif90
# F77_MPI_FLAGS = -Mpreprocess -Mextend -DMPI -O2 -tp amd64e
# F95_MPI_FLAGS = -Mpreprocess -DMPI -O2 -tp amd64e -DPGI
# CC_MPI_FLAGS = -DMPI -DUNDERSCORE -O2 -tp amd64e

#  LAPACK=${ACML}
#  BLAS=${ACML}

# Pathscale compiler
# must unload pgi and load pathscale and acml first:
# module swap PrgEnv-pgi PrgEnv-pathscale && module load acml
# The "-Wl,-z,muldefs" tells it to not stop if it encounters multiple definitions of e.g. __exp

# Serial compilation
 F95 = ftn -target=linux
 CC = cc -target=linux
 F77_FAST_FLAGS = -cpp -DPATHSCALE -O2 -r8 -extend-source -zerouv
 F95_FAST_FLAGS = -cpp -DPATHSCALE -O2 -r8 -zerouv
 CC_FAST_FLAGS =  -DPATHSCALE -DUNDERSCORE -O2
 LD = ftn -target=linux -Wl,-z,muldefs
 LD_FAST_FLAGS = -O2

# MPI compilation
 F95_MPI = ftn -target=linux
 CC_MPI = cc -target=linux
 F77_MPI_FLAGS = -cpp -DPATHSCALE -DMPI -O2 -r8 -extend-source -zerouv
 F95_MPI_FLAGS = -cpp -DPATHSCALE -DMPI -O2 -r8 -zerouv
 CC_MPI_FLAGS = -DPATHSCALE -DMPI -DUNDERSCORE -O2
 LD_MPI = ftn -target=linux -O2 -Wl,-z,muldefs

# acml has both the LAPACK and the BLAS
 LAPACK=-lacml
 BLAS=

endif

# AMD Opteron cluster with Pathscale compiler on slid.fysik at Lyngby, Denmark
ifeq ($(findstring slid.fysik.dtu.dk,$(hostname)), slid.fysik.dtu.dk)

# Serial compilation
#F95 = ftn -target=linux
 F95 = pathf90
#CC = cc -target=linux
 CC = pathf90
 F77_FAST_FLAGS = -cpp -DPATHSCALE -O2 -r8 -extend-source -zerouv
 F95_FAST_FLAGS = -cpp -DPATHSCALE -O2 -r8 -zerouv
 CC_FAST_FLAGS =  -DPATHSCALE -DUNDERSCORE -O2
 LD = pathf90 -Wl,-z,muldefs
 LD_FAST_FLAGS = -O2

# MPI compilation
 F95_MPI = pathf90
 CC_MPI = pathf90
 F77_MPI_FLAGS = -cpp -DPATHSCALE -DMPI -O2 -r8 -extend-source -zerouv
 F95_MPI_FLAGS = -cpp -DPATHSCALE -DMPI -O2 -r8 -zerouv
 CC_MPI_FLAGS = -DPATHSCALE -DMPI -DUNDERSCORE -O2
 LD_MPI = pathf90 -O2 -Wl,-z,muldefs

# acml has both the LAPACK and the BLAS
 LAPACK=-L/opt/acml-3.6.0/pathscale64/lib -lacml
 BLAS=

endif

# AMD Opteron cluster with Intel compiler on school1.fysik at Lyngby, Denmark
ifeq ($(findstring school,$(hostname)), school)

#F77_FAST_FLAGS=-zero -extend_source -w -r8 -O3 -pad -ip -fpp -mtune=core2 -march=core2 -xT
#F77_FLAGS=$(F77_FAST_FLAGS)
#F95_FAST_FLAGS=-O3 -ip -align all -fpp -mtune=core2 -march=core2 -xT
#F95_FLAGS=$(F95_FAST_FLAGS)
#LD_FAST_FLAGS=-O3 -align all -mtune=core2 -march=core2 -xT
#LD_FLAGS=$(LD_FAST_FLAGS)
#CC_MPI=mpicc
#F95_MPI=mpif90
#F77_MPI=mpif90
#LD_MPI=mpif90
#F77_MPI_FLAGS=-zero -extend_source -w -r8 -O3 -pad -ip -fpp -DMPI -mtune=core2 -march=core2 -xT
#F95_MPI_FLAGS=-O3 -ip -align all -fpp -DMPI -mtune=core2 -march=core2 -xT

F95=ifort
F95_MPI=mpif77

F77_FAST_FLAGS=-zero -extend_source -w -r8 -O3 -pad -ip -fpp
F77_IPO=-zero -extend_source -w -r8 -O3 -pad -ipo -fpp
F77_DEBUG_FLAGS=-zero -extend_source -r8 -check all -check nooutput_conversion -debug extended -debug-parameters all -g -traceback -error_limit 1 -fpp -DDEBUG
F77_PROF_FLAGS=-zero -extend_source -w -r8 -O3 -pad -pg -fpp
F77_MPI_FLAGS=-zero -extend_source -w -r8 -O3 -pad -ip -fpp -DMPI
F77_DEBUG_MPI_FLAGS=-zero -extend_source -r8 -check all -check nooutput_conversion -debug extended -debug-parameters all -g -traceback -error_limit 1 -fpp -DMPI -DDEBUG
F77_FLAGS=$(F77_FAST_FLAGS)

F95_FAST_FLAGS=-O3 -ip -align all -fpp
F95_IPO=-O3 -ipo -align all -fast -fpp
F95_DEBUG_FLAGS=-check all -check nooutput_conversion -debug extended -debug-parameters all -g -traceback -error_limit 1 -fpp -DDEBUG
F95_PROF_FLAGS=-O3 -pg -fpp
F95_MPI_FLAGS=-O3 -ip -align all -fpp -DMPI
F95_DEBUG_MPI_FLAGS=-check all -check nooutput_conversion -debug extended -debug-parameters all -g -traceback -error_limit 1 -fpp -DMPI -DDEBUG
F95_FLAGS=$(F95_FAST_FLAGS)

# Linking  ------------------------------------------------------------------

LD=ifort
LD_MPI=mpif77

LD_FAST_FLAGS=-O3 -align all
LD_IPO=-O3 -ipo -align all -fast
LD_DEBUG_FLAGS=
LD_PROF_FLAGS=-O3 -pg
LD_FLAGS=$(LD_FAST_FLAGS)

LD_END=
#LD_END_MPI=-lmpi -lpthread
LD_END_MPI=

# acml has both the LAPACK and the BLAS
 LAPACK=-L/opt/acml-3.6.0/ifort64/lib -lacml
 BLAS=

endif

# IBM BlueGene/L at RPI CCNI
ifeq ($(hostname), bgl16fen)
# serial
 F77 = /opt/ibmcmp/xlf/bg/10.1/bin/blrts_xlf
 F95 = /opt/ibmcmp/xlf/bg/10.1/bin/blrts_xlf95
 CC = /opt/ibmcmp/vac/bg/8.0/bin/blrts_xlc
 CXX = /opt/ibmcmp/vacpp/bg/8.0/bin/blrts_xlC
 LD = /opt/ibmcmp/xlf/bg/10.1/bin/blrts_xlf95

 F77_FAST_FLAGS = -O3 -qstrict -qarch=440 -qtune=440 -qsuffix=cpp=f -WF,-DUNDERSCORE,-DNOQUENCH,-DNOEINSPLINE -qfixed=132
 F95_FAST_FLAGS = -O3 -qstrict -qarch=440 -qtune=440 -qsuffix=cpp=f90 -WF,-DNOQUENCH,-DNOEINSPLINE
 CC_FAST_FLAGS = -O3 -qstrict  -qarch=440 -qtune=440
 LD_FLAGS = -qarch=440 -qtune=440
 CC_FAST_FLAGS = -O3  -qstrict  -qarch=440 -qtune=440

#MPI
 F77_MPI = mpixlf90
 F95_MPI = mpixlf90
 LD_MPI = mpixlf90
 CC_MPI = mpixlc
#F77_MPI_FLAGS =  -O3 -qcheck  -C -qstrict -qarch=440 -qtune=440 -qsuffix=cpp=f -WF,-DUNDERSCORE,-DNOQUENCH,-DNOEINSPLINE,-DMPI -qfixed=132
 F77_MPI_FLAGS =  -O3 -qstrict -qarch=440 -qtune=440 -qsuffix=cpp=f -WF,-DUNDERSCORE,-DNOQUENCH,-DNOEINSPLINE,-DMPI -qfixed=132
 F95_MPI_FLAGS =  -O3 -qstrict -qarch=440 -qtune=440 -qsuffix=cpp=f90  -WF,-DNOQUENCH,-DNOEINSPLINE,-DMPI

#Libraries
 EINSPLINE=
 QUENCH=
 LIBS_MAKE=$(LIBHOME)/lib/libcyrus.a $(LIBHOME)/lib2/linpack/liblinpack.a $(LIBHOME)/lib2/pspline/pspline/libpspline.a $(LIBHOME)/lib2/lapack/liblapack.a

# Math libraries
 BLAS = -L/bgl/local/essl/lib/ -lesslbg

endif

# Julich BlueGene P
ifeq ($(hostname),jugene5)
# MPI
#  F95_MPI = mpixlf2003
#  F77_MPI_FLAGS = -O3 -qstrict -qarch=450 -qtune=450 -qsuffix=cpp=f -WF,-DMPI,-DUNDERSCORE -qfixed=132 -qsuppress=1512-039 -qsuppress=1512-047
#  F95_MPI_FLAGS = -O3 -qstrict -qarch=450 -qtune=450 -qsuffix=cpp=f90  -WF,-DMPI -qsuppress=1512-039 -qsuppress=1512-047
#  CC_MPI = mpixlc
#  CC_MPI_FLAGS = -O3 -qstrict

 F95=mpixlf2003_r
 F95_FAST_FLAGS=-O3 -qstrict -qarch=450 -qtune=450 -qsuffix=cpp=f90  -WF,-DMPI,-DNOEINSPLINE,-DNOQUENCH -qsuppress=1512-039 -qsuppress=1512-047
 F77_FAST_FLAGS=-O3 -qstrict -qarch=450 -qtune=450 -qsuffix=cpp=f -WF,-DMPI,-DUNDERSCORE,-DNOEINSPLINE,-DNOQUENCH -qfixed=132 -qsuppress=1512-039 -qsuppress=1512-047

# MPI
  F95_MPI = mpixlf2003_r
  F95_MPI_FLAGS = -O3 -qstrict -qarch=450 -qtune=450 -qsuffix=cpp=f90  -WF,-DMPI,-DNOEINSPLINE,-DNOQUENCH -qsuppress=1512-039 -qsuppress=1512-047
  F77_MPI_FLAGS = -O3 -qstrict -qarch=450 -qtune=450 -qsuffix=cpp=f -WF,-DMPI,-DUNDERSCORE,-DNOEINSPLINE,-DNOQUENCH -qfixed=132 -qsuppress=1512-039 -qsuppress=1512-047
  F95_DEBUG_MPI_FLAGS = -g -C -qarch=450 -qtune=450 -qsuffix=cpp=f90  -WF,-DMPI,-DNOEINSPLINE,-DNOQUENCH -qsuppress=1512-039 -qsuppress=1512-047
  F77_DEBUG_MPI_FLAGS = -g -C -qarch=450 -qtune=450 -qsuffix=cpp=f -WF,-DMPI,-DUNDERSCORE,-DNOEINSPLINE,-DNOQUENCH -qfixed=132 -qsuppress=1512-039 -qsuppress=1512-047
  CC_MPI = mpixlc
  CC_MPI_FLAGS = -O3 -qstrict -qarch=450 -qtune=450
# no lapack in Julich
  BLAS=-L/opt/ibmmath/essl/4.4/lib -lesslbg
  QUENCH=
  EINSPLINE=
  LIBS_MAKE=$(LIBHOME)/lib/libcyrus.a $(LIBHOME)/lib2/lapack/liblapack.a $(LIBHOME)/lib2/linpack/liblinpack.a $(LIBHOME)/lib2/pspline/pspline/libpspline.a
endif

# BlueGene P at IDRIS (Paris)
ifeq ($(hostname),babel1)

# -qarch=450d makes the code crash with BE_MPI (ERROR):   "killed with signal 7"
# -qhot or -O5 does not seem to improve performance at all

 F95=mpixlf2003_r
 F95_FAST_FLAGS=-O3 -qstrict -qarch=450 -qtune=450 -qsuffix=cpp=f90  -WF,-DMPI,-DNOEINSPLINE,-DNOQUENCH -qsuppress=1512-039 -qsuppress=1512-047
 F77_FAST_FLAGS=-O3 -qstrict -qarch=450 -qtune=450 -qsuffix=cpp=f -WF,-DMPI,-DUNDERSCORE,-DNOEINSPLINE,-DNOQUENCH -qfixed=132 -qsuppress=1512-039 -qsuppress=1512-047

# MPI
  F95_MPI = mpixlf2003_r
  F95_MPI_FLAGS = -O3 -qstrict -qarch=450 -qtune=450 -qsuffix=cpp=f90  -WF,-DMPI,-DNOEINSPLINE,-DNOQUENCH -qsuppress=1512-039 -qsuppress=1512-047
  F77_MPI_FLAGS = -O3 -qstrict -qarch=450 -qtune=450 -qsuffix=cpp=f -WF,-DMPI,-DUNDERSCORE,-DNOEINSPLINE,-DNOQUENCH -qfixed=132 -qsuppress=1512-039 -qsuppress=1512-047
  CC_MPI = mpixlc
  CC_MPI_FLAGS = -O3 -qstrict -qarch=450 -qtune=450
  LAPACK=-L/bglocal/prod/tools_ibm/lib -llapack -L/opt/ibmmath/essl/4.4/lib -lesslbg
  BLAS=-L/opt/ibmmath/essl/4.4/lib -lesslbg
  QUENCH=
  EINSPLINE=
  LIBS_MAKE=$(LIBHOME)/lib/libcyrus.a $(LIBHOME)/lib2/linpack/liblinpack.a $(LIBHOME)/lib2/pspline/pspline/libpspline.a
endif

# BlueGene P at MPG Garching, Munich
ifeq ($(hostname),genius1)

# -qarch=450d makes the code crash with BE_MPI (ERROR):   "killed with signal 7"
# -qhot or -O5 does not seem to improve performance at all

 F95=mpixlf2003_r
 F95_FAST_FLAGS=-O3 -qstrict -qarch=450 -qtune=450 -qsuffix=cpp=f90  -WF,-DMPI,-DNOEINSPLINE,-DNOQUENCH -qsuppress=1512-039 -qsuppress=1512-047
 F77_FAST_FLAGS=-O3 -qstrict -qarch=450 -qtune=450 -qsuffix=cpp=f -WF,-DMPI,-DUNDERSCORE,-DNOEINSPLINE,-DNOQUENCH -qfixed=132 -qsuppress=1512-039 -qsuppress=1512-047

# MPI
  F95_MPI = mpixlf2003_r
  F95_MPI_FLAGS = -O3 -qstrict -qarch=450 -qtune=450 -qsuffix=cpp=f90  -WF,-DMPI,-DNOEINSPLINE,-DNOQUENCH -qsuppress=1512-039 -qsuppress=1512-047
  F77_MPI_FLAGS = -O3 -qstrict -qarch=450 -qtune=450 -qsuffix=cpp=f -WF,-DMPI,-DUNDERSCORE,-DNOEINSPLINE,-DNOQUENCH -qfixed=132 -qsuppress=1512-039 -qsuppress=1512-047
  CC_MPI = mpixlc
  CC_MPI_FLAGS = -O3 -qstrict -qarch=450 -qtune=450
# no lapack in Garching
  BLAS=-L/opt/ibmmath/essl/4.4/lib -lesslbg
  QUENCH=
  EINSPLINE=
  LIBS_MAKE=$(LIBHOME)/lib/libcyrus.a $(LIBHOME)/lib2/lapack/liblapack.a $(LIBHOME)/lib2/linpack/liblinpack.a $(LIBHOME)/lib2/pspline/pspline/libpspline.a
endif

# Mac OS for Intel Core2 machines with Intel Fortran compiler
ifeq ($(uname),Darwin)
    CC=icc
    CXX=icpc
    CC_MPI=mpicc
    F95=ifort
    F95_MPI=mpif90
    F95_FAST_FLAGS=-zero -r8 -m64 -pad -align all -O3 -axSTP -ip -diag-disable vec -fpp -DNOEINSPLINE -shared-intel
    F95_DEBUG_FLAGS=-zero -r8 -m64 -pad -align all -O0 -fpp -ip -check all -check nooutput_conversion -debug extended -debug-parameters all -g -traceback -error_limit 1 -DDEBUG -DNOEINSPLINE -shared-intel
    EINSPLINE=
    LIBS_MAKE= $(LIBHOME)/lib/libcyrus.a $(LIBHOME)/lib2/blas/libblas.a $(LIBHOME)/lib2/lapack/liblapack.a $(LIBHOME)/lib2/linpack/liblinpack.a $(LIBHOME)/lib2/pspline/pspline/libpspline.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench_seq.a
endif

# grads-46 at Duke
ifeq ($(hostname),grads-46)
  F95=ifort
  CC=gcc
  CXX=g++
  CC_MPI=mpicc

  F95_FAST_FLAGS= -zero -r8 -pad -align all -O3 -fpp -ip
  F95_PROF_FLAGS= -zero -r8 -pad -align all -O3 -fpp -pg
  F95_DEBUG_FLAGS=-zero -r8 -pad -align all -O0 -fpp -ip -check all -check nooutput_conversion -debug extended -debug-parameters all -g -traceback -error_limit 1 -DDEBUG

# by default, compilation for f77 is the same than for f95
  F77=$(F95)
  F77_FAST_FLAGS=$(F95_FAST_FLAGS) -extend_source
  F77_DEBUG_FLAGS=$(F95_DEBUG_FLAGS) -extend_source
  F77_PROF_FLAGS=$(F95_PROF_FLAGS) -extend_source

# by default, fast options
  F95_FLAGS=$(F95_FAST_FLAGS)
  F77_FLAGS=$(F77_FAST_FLAGS)

  USE_GCC=yes
endif

# baranger12 at Duke
ifeq ($(hostname),baranger12)
  CC=gcc
  CXX=g++
  CC_MPI=mpicc
  F95=gfortran
  F95_MPI=mpif90
  F95_FAST_FLAGS=-finit-local-zero -O3 -fdefault-real-8 -ffree-line-length-none -x f95-cpp-input -Wall -fbounds-check -DNOEINSPLINE
  F77_FAST_FLAGS=-finit-local-zero -O3 -fdefault-real-8 -ffixed-line-length-132 -x f77-cpp-input -Wall -fbounds-check -DNOEINSPLINE
  F95_SLOW_FLAGS=-finit-local-zero -O0 -fdefault-real-8 -ffree-line-length-none -x f95-cpp-input -Wall -fbounds-check -DNOEINSPLINE
  F77_SLOW_FLAGS=-finit-local-zero -O0 -fdefault-real-8 -ffixed-line-length-132 -x f77-cpp-input -Wall -fbounds-check -DNOEINSPLINE
  F95_DEBUG_FLAGS=$(F95_SLOW_FLAGS) -g -fbounds-check
  F77_DEBUG_FLAGS=$(F77_SLOW_FLAGS) -g -fbounds-check
  F95_PROF_FLAGS=$(F95_FAST_FLAGS) -pg
  F77_PROF_FLAGS=$(F77_FAST_FLAGS) -pg
  LD=gfortran
  LD_FLAGS=-O3
  LD_MPI=mpif90
  LD_DEBUG_FLAGS=-Wall -g -fbounds-check
  LD_PROF_FLAGS=-pg
  EINSPLINE=
  LIBS_MAKE= $(LIBHOME)/lib/libcyrus.a $(LIBHOME)/lib2/blas/libblas.a $(LIBHOME)/lib2/lapack/liblapack.a $(LIBHOME)/lib2/linpack/liblinpack.a $(LIBHOME)/lib2/pspline/pspline/libpspline.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench_seq.a
endif

# lqcd nodes at Duke
ifeq ($(findstring lqcd,$(hostname)),lqcd)
  CC=gcc
  CXX=g++
  CC_MPI=mpicc
  F95=gfortran
  F95_MPI=mpif90
  F95_FAST_FLAGS=-finit-local-zero -O3 -fdefault-real-8 -ffree-line-length-none -x f95-cpp-input -Wall -fbounds-check -DNOEINSPLINE
  F77_FAST_FLAGS=-finit-local-zero -O3 -fdefault-real-8 -ffixed-line-length-132 -x f77-cpp-input -Wall -fbounds-check -DNOEINSPLINE
  F95_SLOW_FLAGS=-finit-local-zero -O0 -fdefault-real-8 -ffree-line-length-none -x f95-cpp-input -Wall -fbounds-check -DNOEINSPLINE
  F77_SLOW_FLAGS=-finit-local-zero -O0 -fdefault-real-8 -ffixed-line-length-132 -x f77-cpp-input -Wall -fbounds-check -DNOEINSPLINE
  F95_DEBUG_FLAGS=$(F95_SLOW_FLAGS) -g -fbounds-check
  F77_DEBUG_FLAGS=$(F77_SLOW_FLAGS) -g -fbounds-check
  F95_PROF_FLAGS=$(F95_FAST_FLAGS) -pg
  F77_PROF_FLAGS=$(F77_FAST_FLAGS) -pg
  LD=gfortran
  LD_FLAGS=-O3
  LD_MPI=mpif90
  LD_DEBUG_FLAGS=-Wall -g -fbounds-check
  LD_PROF_FLAGS=-pg
  EINSPLINE=
  LIBS_MAKE= $(LIBHOME)/lib/libcyrus.a $(LIBHOME)/lib2/blas/libblas.a $(LIBHOME)/lib2/lapack/liblapack.a $(LIBHOME)/lib2/linpack/liblinpack.a $(LIBHOME)/lib2/pspline/pspline/libpspline.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench_seq.a
endif

# condor submit host at Duke
ifeq ($(hostname),condor.internal.phy.duke.edu)
  CC=condor_compile gcc
  CXX=condor_compile g++
  CC_MPI=condor_compile mpicc
  F95=condor_compile gfortran
  F95_MPI=condor_compile mpif90
  F95_FAST_FLAGS=-finit-local-zero -O3 -fdefault-real-8 -ffree-line-length-none -x f95-cpp-input -Wall -fbounds-check -DNOEINSPLINE
  F77_FAST_FLAGS=-finit-local-zero -O3 -fdefault-real-8 -ffixed-line-length-132 -x f77-cpp-input -Wall -fbounds-check -DNOEINSPLINE
  F95_SLOW_FLAGS=-finit-local-zero -O0 -fdefault-real-8 -ffree-line-length-none -x f95-cpp-input -Wall -fbounds-check -DNOEINSPLINE
  F77_SLOW_FLAGS=-finit-local-zero -O0 -fdefault-real-8 -ffixed-line-length-132 -x f77-cpp-input -Wall -fbounds-check -DNOEINSPLINE
  F95_DEBUG_FLAGS=$(F95_SLOW_FLAGS) -g -fbounds-check
  F77_DEBUG_FLAGS=$(F77_SLOW_FLAGS) -g -fbounds-check
  F95_PROF_FLAGS=$(F95_FAST_FLAGS) -pg
  F77_PROF_FLAGS=$(F77_FAST_FLAGS) -pg
  LD=condor_compile gfortran
  LD_FLAGS=-O3
  LD_MPI=condor_compile mpif90
  LD_DEBUG_FLAGS=-Wall -g -fbounds-check
  LD_PROF_FLAGS=-pg
  EINSPLINE=
  LIBS_MAKE= $(LIBHOME)/lib/libcyrus.a $(LIBHOME)/lib2/blas/libblas.a $(LIBHOME)/lib2/lapack/liblapack.a $(LIBHOME)/lib2/linpack/liblinpack.a $(LIBHOME)/lib2/pspline/pspline/libpspline.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench.a $(LIBHOME)/SimulatedAnnealing/quench_anneal/lib/libquench_seq.a
endif

#JT, personnal computer
ifeq ($(hostname),europa)
  CC=gcc
  CXX=g++
  CC_MPI=mpicc
  F95=gfortran
  F95_MPI=mpif90
  F95_FAST_FLAGS=-finit-local-zero -O3 -ffree-line-length-none -x f95-cpp-input -DNOEINSPLINE
  F77_FAST_FLAGS=-finit-local-zero -O3 -ffixed-line-length-132 -x f77-cpp-input -DNOEINSPLINE
  LD=gfortran
  LD_FLAGS=-O3
  LD_MPI=mpif90
  EINSPLINE=
endif

#BM, personnal computer 1
ifeq ($(hostname),antares)
  CC=gcc
  CXX=g++
  CC_MPI=mpicc
  F95=gfortran
  F95_MPI=mpif90
  F95_FAST_FLAGS=-finit-local-zero -O3 -ffree-line-length-none -x f95-cpp-input -DNOEINSPLINE
  F77_FAST_FLAGS=-finit-local-zero -O3 -ffixed-line-length-132 -x f77-cpp-input -DNOEINSPLINE
  LD=gfortran
  LD_FLAGS=-O3
  LD_MPI=mpif90
  EINSPLINE=
endif

#BM, personnal computer 2
ifeq ($(hostname),antaresICS)
  CC=gcc
  CXX=g++
  CC_MPI=mpicc
  F95=gfortran
  F95_MPI=mpif90
  F95_FAST_FLAGS=-finit-local-zero -O3 -ffree-line-length-none -x f95-cpp-input -DNOEINSPLINE
  F77_FAST_FLAGS=-finit-local-zero -O3 -ffixed-line-length-132 -x f77-cpp-input -DNOEINSPLINE
  LD=gfortran
  LD_FLAGS=-O3
  LD_MPI=mpif90
  EINSPLINE=
endif

