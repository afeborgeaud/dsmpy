# makefile for tish
PROGS = mpi-tish
FC	= mpif90
FC2=ifort
FFLAGS	= -O

SRC = calmat.f trialf.f others.f dclisb.f dclisb3.f formpi.f mpi-tish.f
SRC2 = calmat.f trialf.f others.f dclisb.f dclisb3.f tish.f
OBJS	= $(SRC:.f=.o)
OBJS2 =$(SRC2:.f=.o)
.SUFFIXES: .f .o

all:$(PROGS) tish

mpi-tish: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) 

tish: $(OBJS2)
	$(FC2) $(FFLAGS) -o $@ $(OBJS2)

mpi-tish.o: mpi-tish.f
	$(FC) $(FFLAGS) -c mpi-tish.f -o $@


.f.o:
	$(FC2) $(FFLAGS) -c $< 

clean:
	rm -f $(OBJS) $(OBJS2) $(PROGS) tish work
