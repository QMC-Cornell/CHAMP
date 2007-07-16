# Warning: sort_strings.F is compiler dependent because malloc can be
# either call by value or call by reference.
#FC     = pgf90
#FFLAGS = -fast -r8 -Mdalign -Mbounds -Mextend
#FFLAGS = -fast -r8

 FC     = g77
 FFLAGS = -O2 -malign-double -ffixed-line-length-none -Wall
#FFLAGS = -O3 -malign-double -ffixed-line-length-none

#FC     = xlf
#FFLAGS = -gC -qextchk -qsave -qspnans -qdpc -bloadmap:map
#FFLAGS =  -O2 -qextchk -qsave -qspnans -qdpc -bloadmap:map

OBJ    = io_test.o io.o sort_strings.o

io_test: $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) -o $@
