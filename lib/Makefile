include ../makefile.inc

FC     = $(CHAMP_F77)
FFLAGS = $(CHAMP_F77_FLAGS)
LDFLAGS= $(CHAMP_LD_FLAGS)
CC     = $(CHAMP_C)
CXX    = $(CHAMP_CXX)

export FC FFLAGS LDFLAGS CC CXX

# PHONY is a Make keyword that tells it to execute the rule regardless of
# the dependencies.  We use it here for 2 reasons:
# 1) prevents it from getting confused if there happens to be a file named
#    "clean" or "clean_all" in the directory.
# 2) Also, for the libraries, linpack etc. it does the make even though there
#    are directories named linpack etc.
.PHONY: clean clean_all


# In quench_anneal lib/libquench_seq.a also needs to be made, but it tags along with lib/libquench.a
all: lib/libcyrus.a lib2/blas/libblas.a lib2/lapack/liblapack.a lib2/linpack/liblinpack.a lib2/pspline/pspline/libpspline.a lib2/einspline/lib/libeinspline.a SimulatedAnnealing/quench_anneal/lib/libquench.a

lib/libcyrus.a:
	cd lib ; $(MAKE)

lib2_dir:
	cd lib2 ; $(MAKE)

lib2/blas/libblas.a:
	cd lib2/blas ; $(MAKE) libblas.a

lib2/lapack/liblapack.a:
	cd lib2/lapack ; make liblapack.a

lib2/linpack/liblinpack.a:
	cd lib2/linpack ; make liblinpack.a

lib2/pspline/pspline/libpspline.a:
	cd lib2/pspline/pspline ; make libpspline.a

ifdef EINSPLINE
ifndef USE_GCC
lib2/einspline/lib/libeinspline.a:
#	cd einspline ; configure_my ; $(MAKE) ; $(MAKE) install
#	cd einspline ; ./configure "CC=$(CC) -std=c99" "CXX=$(CXX) -std=c99" "FC=$(CHAMP_F77)" "F77=$(CHAMP_F77)" --prefix=$(shell pwd)/einspline/ "--enable-sse" ; $(MAKE) ; $(MAKE) install
	cd lib2/einspline ; ./configure "CC=$(CC) -std=c99" "CXX=$(CXX) -std=c99" "FC=$(CHAMP_F77)" "F77=$(CHAMP_F77)" --prefix=$(shell pwd)/einspline/ ; $(MAKE) ; $(MAKE) install ; rm lib/*.so*
endif
ifdef USE_GCC
lib2/einspline/lib/libeinspline.a:
	cd lib2/einspline ; ./configure "CC=$(CC)" "FC=$(CHAMP_F77)" "F77=$(CHAMP_F77)" --prefix=$current_dir --exec-prefix=$current_dir --prefix=$(shell pwd)/einspline/ ; $(MAKE) ; $(MAKE) install ; rm lib/*.so*
endif
endif

SimulatedAnnealing/quench_anneal/lib/libquench.a:
	cd SimulatedAnnealing/quench_anneal ; $(MAKE) libs


# clean just .o files
clean: clean_lib_dir clean_lib2_dir clean_quench

clean_lib_dir:
	cd lib ; $(MAKE) clean

clean_lib2_dir:
	cd lib2/blas ; $(MAKE) clean
	cd lib2/lapack ; $(MAKE) clean
	cd lib2/linpack ; $(MAKE) clean
	cd lib2/pspline/pspline ; $(MAKE) clean
	cd lib2/einspline ; $(MAKE) clean

clean_quench:
	cd SimulatedAnnealing/quench_anneal ; $(MAKE) clean


# clean .o and .a files
clean_all: clean_all_lib_dir clean_all_lib2_dir clean_all_quench

clean_all_lib_dir:
	cd lib ; $(MAKE) clean_all

clean_all_lib2_dir:
	cd lib2/blas ; $(MAKE) clean_all
	cd lib2/lapack ; $(MAKE) clean_all
	cd lib2/linpack ; $(MAKE) clean_all
	cd lib2/pspline/pspline ; $(MAKE) clean_all
	cd lib2/einspline ; $(MAKE) distclean ; rm lib/*.a

ccccc:
	cd lib2/einspline ; $(MAKE) distclean ; rm lib/*.a

clean_all_quench:
	cd SimulatedAnnealing/quench_anneal ; $(MAKE) clean_all
