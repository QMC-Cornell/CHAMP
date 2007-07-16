FC     = g77
FFLAGS = -g -C -O2 -malign-double -ffixed-line-length-none -fno-automatic -Wall -fbounds-check

OBJS   = main.o comp2real.o filem.o order.o

comp2real: $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o $@
