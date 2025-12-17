F90 = mpif90
FC = gfortran

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

JSON_FORTRAN_HOME ?= /usr
JSONF_LIBS := $(JSON_FORTRAN_HOME)/lib -ljsonfortran
JSONF_INCLUDES := $(JSON_FORTRAN_HOME)/include

MOSAIC_LIB := ../mosaic/libmosaic.a

RPATH := -Wl,-rpath,$(JSON_FORTRAN_HOME)/lib

LIBS := $(NFLIBS) $(NCLIBS) $(MOSAIC_LIB) $(GSL_LIBS) -L$(JSONF_LIBS) $(RPATH)
INC := -I$(NCINCLUDES) -I../partmc/src -I$(GSL_HOME)/include -I$(NFINCLUDES) -I$(JSONF_INCLUDES)
FFLAGS_DEBUG := -g -fimplicit-none -W -Wall -Wconversion -Wunderflow -Wno-unused -fbounds-check
FFLAGS := -O3 $(NFFLAGS)
#FFLAGS= -g -fimplicit-none -W -Wall -Wconversion -Wunderflow -Wimplicit-interface -Wno-unused -fbounds-check -fconvert=big-endian

OBJ=  \
  ../partmc/src/util.o \
  ../partmc/src/sys.o \
  ../partmc/src/constants.o \
  ../partmc/src/netcdf.o  \
  ../partmc/src/aero_component.o \
  ../partmc/src/aero_dist.o \
  ../partmc/src/aero_mode.o \
  ../partmc/src/aero_data.o \
  ../partmc/src/fractal.o   \
  ../partmc/src/env_state.o \
  ../partmc/src/spec_file.o \
  ../partmc/src/spec_line.o \
  ../partmc/src/mpi.o \
  ../partmc/src/sort.o \
  ../partmc/src/rand_gsl.o \
  ../partmc/src/chamber.o \
  ../partmc/src/scenario.o \
  ../partmc/src/aero_state.o \
  ../partmc/src/gas_data.o \
  ../partmc/src/gas_state.o \
  ../partmc/src/aero_weight.o \
  ../partmc/src/rand.o \
  ../partmc/src/bin_grid.o \
  ../partmc/src/aero_particle.o \
  ../partmc/src/aero_particle_array.o \
  ../partmc/src/aero_weight_array.o \
  ../partmc/src/aero_info_array.o \
  ../partmc/src/aero_info.o \
  ../partmc/src/aero_binned.o \
  ../partmc/src/aero_sorted.o \
  ../partmc/src/integer_varray.o \
  ../partmc/src/integer_rmap2.o  \
  ../partmc/src/stats.o

WRF_SRC_ROOT_DIR=../WRFV3/
LIB_BUNDLED     = \
                      $(WRF_SRC_ROOT_DIR)/external/fftpack/fftpack5/libfftpack.a \
                      $(WRF_SRC_ROOT_DIR)/external/io_grib1/libio_grib1.a \
                      $(WRF_SRC_ROOT_DIR)/external/io_grib_share/libio_grib_share.a \
                      $(WRF_SRC_ROOT_DIR)/external/io_int/libwrfio_int.a \
                      $(ESMF_IO_LIB) \
                      $(WRF_SRC_ROOT_DIR)/external/RSL_LITE/librsl_lite.a \
                      $(WRF_SRC_ROOT_DIR)/frame/module_internal_header_util.o \
                      $(WRF_SRC_ROOT_DIR)/frame/pack_utils.o  \
                      $(WRF_SRC_ROOT_DIR)/../mosaic/libmosaic.a

%.o: %.mod

make_emission:
	$(F90) make_emissions.F90 $(INC) $(LIBS) $(OBJ) $(FFLAGS) -o create_emissions
	$(F90) compute_emission_statistics.F90 $(INC) $(LIBS) $(OBJ) $(FFLAGS) -o compute_stats
convert_emiss: convert_emiss.o
	$(RANLIB) ../WRFV3/main/libwrflib.a
	$(FC) -o convert_emiss.exe convert_emiss.o ../WRFV3/main/libwrflib.a $(LIB_BUNDLED)

clean:
	rm -rf create_emissions
