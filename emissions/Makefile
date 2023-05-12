F90 = mpif90
FC = gfortran

LIBS=-L$(NETCDF)/lib/ -lnetcdf -lnetcdff /data/keeling/a/jcurtis2/d/camp/camp_deps/jsonfortran-gnu-6.1.0/lib/libjsonfortran.a -L$(GSL_HOME)/lib -lgsl -lgslcblas -lm
INC= -I$(NETCDF)/include/ -I../partmc/src  -I/data/keeling/a/jcurtis2/d/camp/camp_deps/jsonfortran-gnu-6.1.0/include -I$(GSL_HOME)/include
FFLAGS= -g -fimplicit-none -W -Wall -Wconversion -Wunderflow -Wimplicit-interface -Wno-unused -fbounds-check -fconvert=big-endian

OBJ=  \
  ../partmc/src/util.o \
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
	$(F90) make_emissions_ideal.F90 $(INC) $(LIBS) $(OBJ) $(FFLAGS) -o create_emissions_ideal
	$(F90) compute_emission_statistics.F90 $(INC) $(LIBS) $(OBJ) $(FFLAGS) -o compute_stats
convert_emiss: convert_emiss.o
	$(RANLIB) ../WRFV3/main/libwrflib.a
	$(FC) -o convert_emiss.exe convert_emiss.o ../WRFV3/main/libwrflib.a $(LIB_BUNDLED)

clean:
	rm -rf create_emissions
