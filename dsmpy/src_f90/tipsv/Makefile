# makefile .fr tipsv
PROGS = tipsv
FC= gfortran
FFLAGS = -O
SRC = parameters.f90 calmat.f90 trialf.f90 others.f90 rk3.f90 \
      glu2.f90  dcsymbdl.f90 tipsv.f90
OBJS	= $(SRC:.f90=.o)
.SUFFIXES: .f90 .o

all: $(PROGS) 

tipsv: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

.f90.o:
	$(FC) $(FFLAGS) -c $< 

%.o: %.mod

clean:
	rm -f $(OBJS) $(PROGS) parameters.mod
