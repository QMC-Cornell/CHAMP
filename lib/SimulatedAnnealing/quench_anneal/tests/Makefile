# This makefile is normally executed by the one in the parent directory.

LIBS = ../lib/libquench.a ../lib/libquench_seq.a ../lib/liblinpack.a ../lib/liblapack.a ../lib/libblas.a
LIBS_MPI = ../lib/libquench.a ../lib/libquench_mpi.a ../lib/liblinpack.a ../lib/liblapack.a ../lib/libblas.a

test_seq: test.F functions.f $(LIBS)
	$(FC) -c $(FFLAGS) test.F functions.f -I../include
	$(FC) $(LDFLAGS) -o $@ test.o functions.o -I../include -L../lib -lquench -lquench_seq -llinpack -llapack -lblas
	-rm test.o functions.o
	@echo Made test_seq

test_mpi: test.F functions.f $(LIBS_MPI)
	$(FC_MPI) -c test.F functions.f -I../include $(FFLAGS_MPI)
	$(FLINKER_MPI) $(LDFLAGS) -o $@ test.o functions.o -I../include -L../lib -lquench -lquench_mpi -llinpack -llapack -lblas
	-rm test.o functions.o
	@echo Made test_mpi

clean:
	rm -f *.o

clean_all:
	rm -f *.o *.a
