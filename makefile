# makefile for tish
PROGS = tish
FC = gfortran
FFLAGS	= -O

SRC = calmat.f90 trialf.f90 others.f90 dclisb.f90 dclisb3.f90 tish.f90 parameters.f90
OBJS	= $(SRC:.f90=.o)
.SUFFIXES:

all:$(PROGS)

tish: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

parameters.mod: parameters.o parameters.f90
	$(FC) -c $(FFLAGS) parameters.f90

parameters.o: parameters.f90
	$(FC) -c $(FFLAGS) parameters.f90

%.f90:
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJS) $(PROGS)
