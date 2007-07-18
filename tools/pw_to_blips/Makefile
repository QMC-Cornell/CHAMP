#------------------------------------------------------------#
# Makefile for the blip converter.                           #
#                                                            #
# MDT 9.2003                                                 #
# Gutted for CHAMP by William Parker 1.2007		     #
#------------------------------------------------------------#

include ../../qmc/makefile.inc

default: pw2b

clean:
	rm -rf pw2b *.o *.mod

OBJ = singleton.o  blip.o 

pw2b : $(OBJ)
	$(CHAMP_LD) $(LDFLAGS) -o pw2b \
        $(OBJ) $(LDLIBS)

$(OBJ) : $(PATTERN_RULE)
	$(CHAMP_F95) -c $(FFLAGS) $(@F:.o=.f90)

