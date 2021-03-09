# name of compiler
FC = gfortran

# optimization flags
#gfortran:
OPTIM	= -O3
DEBUG  =  -O0 -std=f2003 -g -fimplicit-none -Wall -pedantic
DEBUGMIN  =  -O0 -g -Wall -pedantic

OPTIONS = $(OPTIM)
#OPTIONS = $(DEBUG)
#OPTIONSMIN = $(OPTIM)
#OPTIONSMIN = $(DEBUGMIN)

LLIBS = -lblas -llapack


OBJS    = diag_mod.o nrtype.o
MODS    = nrtype.mod diag_mod.mod

swave_sample: swave_sample.o $(MODS)
	$(FC) -o swave_sample swave_sample.o $(OBJS)  $(LLIBS)

swave_sample.o: swave_sample.f90 $(MODS) Makefile
	$(FC) $(OPTIONS) -c swave_sample.f90

diag_mod.o diag_mod.mod : diag_mod.f90 nrtype.mod Makefile
	$(FC) $(OPTIONS) -c diag_mod.f90

nrtype.o nrtype.mod : nrtype.f90 Makefile
	$(FC) $(OPTIONS) -c nrtype.f90

clean : 
	rm $(OBJS) $(MODS)

