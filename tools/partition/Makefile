#FC     = pgf77
#FFLAGS = -fast -r8 -Mextend -Msave -Mdalign -Mbounds
#FFLAGS = -fast -r8 -Mextend -Msave

 FC     = g77
#FFLAGS = -O2 -malign-double -ffixed-line-length-none -fno-automatic -Wall
 FFLAGS = -O3 -malign-double -ffixed-line-length-none -fno-automatic

#FC     = xlf
#FFLAGS = -gC -qextchk -qsave -qspnans -qdpc -bloadmap:map
#FFLAGS =  -O2 -qextchk -qsave -qspnans -qdpc -bloadmap:map

partition: partition.o
	$(FC) $(FFLAGS) partition.o -o $@
