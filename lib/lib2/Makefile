AR = ar
ARFLAGS = r

###################################
# Targets, Dependencies and Rules #
###################################

# PHONY is a Make keyword that tells it to execute the rule regardless of
# the dependencies.  We use it here for 2 reasons:
# 1) prevents it from getting confused if there happens to be a file named
#    "clean" or "clean_all" in the directory.
# 2) Also, for the libraries, linpack etc. it does the make even though there
#    are directories named linpack etc.
.PHONY: linpack lapack blas pspline einspline clean clean_all

all: linpack lapack blas pspline einspline

linpack:
	cd linpack ; $(MAKE) liblinpack.a

lapack:
	cd lapack ; $(MAKE) liblapack.a

blas:
	cd blas ; $(MAKE) libblas.a

pspline:
	cd pspline/pspline ; $(MAKE) libpspline.a

einspline:
	cd einspline ; configure_my ; $(MAKE) ; $(MAKE) install

clean:
	-rm -f linpack/*.o lapack/*.o blas/*.o pspline/pspline/*.o einspline/src/*.o
clean_all:
	-rm -f linpack/*.o lapack/*.o blas/*.o pspline/pspline/*.o einspline/src/*.o \
            linpack/liblinpack.a lapack/liblapack.a pspline/pspline/libpspline.a einspline/lib/lib*
