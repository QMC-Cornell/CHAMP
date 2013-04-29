#------------------------------------------------------------#
# Makefile for the blip converter.                           #
#                                                            #
# modified by wparker 2005.12.15
# MDT 9.2003                                                 #
#------------------------------------------------------------#

# System dependent include file
#include ./$(QMC_ARCH).inc
# User overrides
#include ../../zmakes/users/$(QMC_ARCH)/$(QMC_ID).inc

MPF90 = ifort
MPLD = ifort
BINDIR   = ./
LIBDIR   = ./
#LIBDIR   = ./zlib/$(QMC_ARCH)

default: dirs $(BINDIR)/blip

dirs:
	@if [ ! -d $(BINDIR) ]; then mkdir -p $(BINDIR) ; fi
	@if [ ! -d $(LIBDIR) ]; then mkdir -p $(LIBDIR) ; fi

clean:
	rm -rf zlib/$(QMC_ARCH) $(BINDIR)/blip *~

OBJ =   $(LIBDIR)/singleton.o  $(LIBDIR)/blip.o

$(BINDIR)/blip : $(OBJ)
	$(MPLD) $(LDFLAGS) -o $(BINDIR)/blip \
	$(OBJ) $(LDLIBS)

$(OBJ) : $(PATTERN_RULE)
	cd $(LIBDIR) ; $(MPF90) -c $(FFLAGS) $(@F:.o=.f90)

