# makefile .fr tipsv
PROGS = mpi-tipsv
FC	= mpif90
FC2= gfortran
FFLAGS	= -O
SRC = calmat.f90 trialf.f90 others.f90 rk3.f90 glu2.f90  dcsymbdl.f90 \
	mpi-tipsv.f90 formpi.f90
OBJS	= $(SRC:.f90=.o)
.SUFFIXES: .f90 .o

all:$(PROGS) 

tipsv: tipsv.f90
	$(FC2) $(FFLAGS) tipsv.f90 -o tipsv $(OBJS)

mpi-tipsv: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

mpi-tipsv.o: mpi-tipsv.f90
	$(FC) $(FFLAGS) -c mpi-tipsv.f90 -o $@

.f90.o:
	$(FC2) $(FFLAGS) -c $< 

clean:
	rm -f $(OBJS) $(PROGS) work
