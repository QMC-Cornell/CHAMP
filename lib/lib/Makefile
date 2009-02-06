# OBJS = my_second.o ranf.o shellnn.o intpol.o intpol_monot.o grater.o \
# hartre.o deriv.o deriv2.o divergence.o lapla3.o laplac.o fderiv.o simson.o bode.o \
# hartre_grid2.o deriv_grid2.o deriv2_grid2.o divergence_grid2.o laplac_grid2.o fderiv_grid2.o \
# bode_grid2.o

OBJS = my_second.o intpol.o simson.o

libcyrus.a: $(OBJS)
	ar r $@ $(OBJS)
	ranlib $@

clean:
	-rm -f *.o
	-rm -f *.a
