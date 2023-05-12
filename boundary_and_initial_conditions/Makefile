FC = mpif90
LIBS=-L$(NETCDF)/lib/ -lnetcdf -lnetcdff -L../partmc/src/libmosaic.a -L$(GSL_HOME)/lib -lgsl -lgslcblas -lm
INC= -I$(NETCDF)/include/ -I../partmc/src -I$(GSL_HOME)/include
FFLAGS_DEBUG = -g -fimplicit-none -W -Wall -Wconversion -Wunderflow -Wno-unused -fbounds-check
FFLAGS = -O3 #$(FFLAGS_DEBUG)

MODULES =  \
  ../partmc/src/util.o \
  ../partmc/src/rand_gsl.o \
  ../partmc/src/sort.o \
  ../partmc/src/constants.o \
  ../partmc/src/netcdf.o  \
  ../partmc/src/aero_dist.o \
  ../partmc/src/aero_mode.o \
  ../partmc/src/aero_data.o \
  ../partmc/src/env_state.o \
  ../partmc/src/fractal.o   \
  ../partmc/src/spec_file.o \
  ../partmc/src/spec_line.o \
  ../partmc/src/aero_component.o \
  ../partmc/src/mpi.o \
  ../partmc/src/gas_data.o \
  ../partmc/src/gas_state.o \
  ../partmc/src/rand.o \
  ../partmc/src/bin_grid.o \
  ../partmc/src/aero_weight.o \
  ../partmc/src/aero_particle.o \

%.o: %.mod

all: make_boundary_conditions make_initial_conditions

make_boundary_conditions: $(MODULES)
	$(FC) make_bcs.F90 $(MODULES) $(INC) $(LIBS) $(FFLAGS) -o create_bcs

make_initial_conditions: $(MODULES)
	$(FC) make_ics.F90 $(MODULES) $(INC) $(LIBS) $(FFLAGS) -o create_ics
#	$(FC) make_ics_ideal.F90 $(MODULES) $(INC) $(LIBS) $(FFLAGS) -o create_ics_ideal
clean:
	rm -rf create_bcs create_ics
