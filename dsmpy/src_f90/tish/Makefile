# makefile for tish
PROGS = tish
FC = gfortran
FFLAGS	= -O

SRC = parameters.f90 calmat.f90 trialf.f90 others.f90 \
      dclisb.f90 tish.f90
OBJS	= $(SRC:.f90=.o)
.SUFFIXES: .f90 .o

all: $(PROGS)

tish: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

.f90.o:
	$(FC) $(FFLAGS) -c $< 

%.o: %.mod

clean:
	rm -f $(OBJS) $(PROGS) parameters.mod
