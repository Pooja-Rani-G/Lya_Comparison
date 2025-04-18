# Makefile for HOP.  Follow the checklist:
#
# 1) Compiler choice:
CC = gcc
FC = gfortran -ffree-line-length-none -cpp

# 2) Compiler flags: optimization, debugging, profiling, etc.
# For debugging:
#CFLAGS	=   -g
# For optimization:
# SGI:
#CFLAGS	=   -O3
# Solaris:
#CFLAGS	=   -fast -64 -fsimple=2 -fsingle
CFLAGS	=   -O -mcmodel=medium
# 3) Libraries: (if you have some replacement for math.h)
LIBS	=   -lm

# 4) Preprocessor definitions: Useful to toggle the DIFFERENT_MASSES flag
DEFS = -DDIFFERENT_MASSES

# 5) Return to the command line and type "make"

###############

default:	hop regroup poshalo evol isol

clean:
	rm -f hop.o hop_input.o kd.o smooth.o regroup.o slice.o poshalo.o

isol: halo_isolation.o
	g++ halo_isolation.o -o isol

halo_isolation.o: halo_isolation.cc
	g++ -c halo_isolation.cc

evol: halo_evolution.o
	$(FC) halo_evolution.o -o evol

halo_evolution.o : halo_evolution.f90
	$(FC) -c halo_evolution.f90

poshalo: poshalo.o
	$(FC) poshalo.o -o poshalo

poshalo.o: poshalo.f90
	$(FC) -c poshalo.f90

hop: hop.o kd.o smooth.o hop_input.o
	$(CC) $(CFLAGS) -o hop hop.o kd.o smooth.o hop_input.o $(LIBS)

hop.o: hop.c kd.h smooth.h
	$(CC) $(CFLAGS) $(DEFS) -c hop.c

kd.o: kd.c kd.h
	$(CC) $(CFLAGS) $(DEFS) -c kd.c

smooth.o: smooth.c kd.h smooth.h
	$(CC) $(CFLAGS) $(DEFS) -c smooth.c

hop_input.o: hop_input.c kd.h
	$(CC) $(CFLAGS) $(DEFS) -c hop_input.c

regroup: regroup.o slice.o
	$(CC) $(CFLAGS) -o regroup regroup.o slice.o $(LIBS)

regroup.o: regroup.c slice.h
	$(CC) $(CFLAGS) -c regroup.c

slice.o: slice.c slice.h
	$(CC) $(CFLAGS) -c slice.c
