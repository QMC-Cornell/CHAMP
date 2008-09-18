# FC     = g77
# FFLAGS = -g -C -O2 -malign-double -ffixed-line-length-none -fno-automatic -Wall -fbounds-check

  FC = ifort
  FFLAGS = -zero -extend_source -r8 -check all -check nooutput_conversion -debug extended -debug-parameters all -g -traceback -error_limit 1 -fpp -DDEBUG
# FFLAGS = -zero -save -extend_source -w -r8 -O3 -pad

OBJS   = main.o comp2real.o filem.o order.o

comp2real: $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o $@
