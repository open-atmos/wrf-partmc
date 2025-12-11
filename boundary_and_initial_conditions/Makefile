FC = mpif90

# Find NetCDF installation. NETCDF should be set to compile WRF.
# Use nc-config and nf-config to get includes and libraries
NFCONFIG := $(NETCDF)/bin/nf-config
ifeq ($(wildcard $(NFCONFIG)),)
  $(error nf-config not found in specified path. Check NFCONFIG location.)
endif
NFFLAGS := $(shell "$(NFCONFIG)" --fflags)
NFINCLUDES := $(shell "$(NFCONFIG)" --includedir)
NFLIBS := $(shell "$(NFCONFIG)" --flibs)

NCCONFIG := $(NETCDF)/bin/nc-config
ifeq ($(wildcard $(NCCONFIG)),)
  $(error nc-config not found in specified path. Check NCCONFIG location.)
endif
NCLIBS := $(shell $(NCCONFIG) --libs)
NCINCLUDES := $(shell $(NCCONFIG) --includedir)

# Set GSL
GSL_HOME ?= /usr
GSL_INCLUDES = -I$(GSL_HOME)/include
# Default to empty
GSL_LIBDIR :=
ifneq ($(wildcard $(GSL_HOME)/lib64),)
  GSL_LIBDIR := $(GSL_HOME)/lib64
else ifneq ($(wildcard $(GSL_HOME)/lib),)
  GSL_LIBDIR := $(GSL_HOME)/lib
else
  $(error No GSL library directory found. Checked $(GSL_HOME)/lib64 and $(GSL_HOME)/lib)
endif
GSL_LIBS := -L$(GSL_LIBDIR) -lgsl -lgslcblas

MOSAIC_LIB := ../mosaic/libmosaic.a

LIBS := $(NFLIBS) $(NCLIBS) $(MOSAIC_LIB) $(GSL_LIBS)
INC := -I$(NCINCLUDES) -I../partmc/src -I$(GSL_HOME)/include -I$(NFINCLUDES)
FFLAGS_DEBUG := -g -fimplicit-none -W -Wall -Wconversion -Wunderflow -Wno-unused -fbounds-check
FFLAGS := -O3 $(NFFLAGS)

MODULES =  \
  ../partmc/src/util.o \
  ../partmc/src/sys.o \
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
	$(FC) -c ic_bc_helper.F90 -I../partmc/src
	$(FC) make_bcs.F90 ic_bc_helper.o $(MODULES) $(INC) $(LIBS) $(FFLAGS) -o create_bcs

make_initial_conditions: $(MODULES)
	$(FC) make_ics.F90 ic_bc_helper.o $(MODULES) $(INC) $(LIBS) $(FFLAGS) -o create_ics
clean:
	rm -rf create_bcs create_ics
