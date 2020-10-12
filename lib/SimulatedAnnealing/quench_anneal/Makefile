# I do not know why but when this is called from champ/qmc/Makefile and we are compiling there
# in MPI mode, using "make mpi" then it insists on using mpif90 instead of resetting the compiler
# to gfortran.  Anyway, this does not seem to matter but we need the "FFLAGS = -DNOMPI" to prevent
# quench_anneal from being compiled in MPI mode when CHAMP is being compiled in MPI mode.
# Alternatively one can do a "make libs" from here before doing a "make mpi" in champ/qmc.

include ../../../makefile.inc

ifndef FFLAGS
FC = gfortran
#FFLAGS=-finit-local-zero -O0 -ffixed-line-length-132 -x f77-cpp-input -Wall -fbounds-check -g -fbacktrace
FFLAGS=-finit-local-zero -O3 -ffixed-line-length-132 -x f77-cpp-input -Wall -fbounds-check
endif

export FC FFLAGS LDFLAGS FC_MPI FFLAGS_MPI LDFLAGS_MPI
 
##########
# README #
##########

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

$(shell mkdir -p lib)
# system = AIX
export F77 FFLAGS LDFLAGS F77_MPI FFLAGS_MPI FLINKER_MPI

# export FC F77 FFLAGS LDFLAGS F77_MPI FFLAGS_MPI FLINKER_MPI

# In most cases the Fortran compilers used for sequential and MPI programs are
# in fact the same, but many implementations of MPI provide a script (often called
# something like mpif77 of hf77) to call the Fortran compiler with the correct
# options, libraries and includes for MPI programs.  If your script does not have
# an include to mpif.h then you will need to put the path of your system's mpif.h
# INCLUDE_MPI variable.  Typically the system mpif.h is in a directory such
# as /usr/local/mpi/include.


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
