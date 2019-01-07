# ======================================================================
# This is a make file
# ======================================================================

# The compiler
FC = gfortran

# flags for debugging or for maximum performance, comment as necessary
FCFLAGS = -g -fbounds-check
FCFLAGS = -O2

# flags for all (e.g. look for system .mod files, required in gfortran)
FCFLAGS += -I/usr/include


# List of executables to be built within the package
PROGRAMS = dsm

# "make" builds all
all: $(PROGRAMS)

io.o: constants.o
parameters.o: constants.o io.o
utilities.o: constants.o io.o parameters.o
qustate.o :constants.o io.o parameters.o utilities.o
wenst.o :constants.o io.o parameters.o utilities.o qustate.o
dsm.o: constants.o io.o parameters.o utilities.o qustate.o wenst.o
dsm: dsm.o constants.o io.o parameters.o utilities.o qustate.o wenst.o

# executable depends
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	rm -rf *.o *.mod *.MOD
	mv dsm ../obj/dsm


all: clean

