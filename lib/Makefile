include ../makefile.inc

FC     = $(CHAMP_F77)
FFLAGS = $(CHAMP_F77_FLAGS)
LDFLAGS= $(CHAMP_LD_FLAGS)
CC     = $(CHAMP_C)
CXX    = $(CHAMP_CXX)

export FC FFLAGS LDFLAGS CC CXX

# PHONY is a Make keyword that prevents it from getting confused if there
# happens to be a file named "clean" in the directory.
.PHONY: clean clean_all


#all: lib_dir lib2_dir quench_dir
all: lib_dir lib2_dir quench

lib_dir:
	cd lib ; $(MAKE)

lib2_dir:
	cd lib2 ; $(MAKE)

quench:
	cd SimulatedAnnealing/quench_anneal ; $(MAKE) libs

clean: clean_lib_dir clean_lib2_dir clean_quench

clean_lib_dir:
	cd lib ; $(MAKE) clean

clean_lib2_dir:
	cd lib2 ; $(MAKE) clean

clean_quench:
	cd SimulatedAnnealing/quench_anneal ; $(MAKE) clean
