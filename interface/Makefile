#

LN      =       ln -sf
MAKE    =       make -i -r
RM      =       rm -f

PMC_MODULES =                           \
  ../partmc/src/aero_binned.o                 \
  ../partmc/src/aero_component.o              \
  ../partmc/src/aero_data.o                   \
  ../partmc/src/aero_dist.o                   \
  ../partmc/src/aero_info.o                   \
  ../partmc/src/aero_info_array.o             \
  ../partmc/src/aero_mode.o                   \
  ../partmc/src/aero_particle.o               \
  ../partmc/src/aero_particle_array.o         \
  ../partmc/src/aero_state.o                  \
  ../partmc/src/aero_sorted.o                 \
  ../partmc/src/aero_weight.o                 \
  ../partmc/src/aero_weight_array.o           \
  ../partmc/src/bin_grid.o                    \
  ../partmc/src/chamber.o                     \
  ../partmc/src/coag_kernel.o                 \
  ../partmc/src/coag_kernel_additive.o        \
  ../partmc/src/coag_kernel_brown.o           \
  ../partmc/src/coag_kernel_brown_cont.o      \
  ../partmc/src/coag_kernel_brown_free.o      \
  ../partmc/src/coag_kernel_constant.o        \
  ../partmc/src/coag_kernel_sedi.o            \
  ../partmc/src/coag_kernel_zero.o            \
  ../partmc/src/coagulation.o                 \
  ../partmc/src/coagulation_dist.o            \
  ../partmc/src/condense.o                    \
  ../partmc/src/constants.o                   \
  ../partmc/src/env_state.o                   \
  ../partmc/src/fractal.o                     \
  ../partmc/src/gas_data.o                    \
  ../partmc/src/gas_state.o                   \
  ../partmc/src/integer_rmap.o                \
  ../partmc/src/integer_rmap2.o               \
  ../partmc/src/integer_varray.o              \
  ../partmc/src/mosaic.o                      \
  ../partmc/src/mpi.o                         \
  ../partmc/src/netcdf.o                      \
  ../partmc/src/nucleate.o                    \
  ../partmc/src/output.o                      \
  ../partmc/src/rand.o                        \
  ../partmc/src/rand_gsl.o                    \
  ../partmc/src/scenario.o                    \
  ../partmc/src/sort.o                        \
  ../partmc/src/spec_file.o                   \
  ../partmc/src/spec_line.o                   \
  ../partmc/src/stats.o                       \
  ../partmc/src/sys.o                         \
  ../partmc/src/util.o

MODULES =                           \
  wrf_pmc_dep.o                      \
  wrf_pmc_dep_aero.o                 \
  wrf_pmc_trans_aero.o               \
  wrf_pmc_driver.o                   \
  wrf_pmc_trans.o                    \
  wrf_pmc_init.o

WRF_MODULES = \
  ../WRFV3/frame/module_wrf_error.o \
  ../WRFV3/frame/module_domain_type.o \
  ../WRFV3/frame/module_domain.o \
  ../WRFV3/share/module_model_constants.o \
  ../WRFV3/frame/module_configure.o \
  ../WRFV3/frame/module_state_description.o \
  ../WRFV3/frame/module_domain.o \
  ../WRFV3/frame/module_configure.o

OBJS = mem.o \
 
LIBTARGET    =  partmc_interface
TARGETDIR    =  ./
$(LIBTARGET) :  $(MODULES) $(PMC_MODULES) $(OBJS)
		$(AR) $(ARFLAGS) ../WRFV3/main/$(LIBWRFLIB) $(MODULES) $(PMC_MODULES) \
		$(OBJS)

include ../WRFV3/configure.wrf

all: partmc_interface

partmc_interface: $(MODULES) $(OBJS)

clean:
	@ echo 'use the clean script'

# DEPENDENCIES : only dependencies after this line (don't remove the word DEPENDENCIES)

wrf_pmc_trans_aero.o: $(WRF_MODULES)

wrf_pmc_driver.o: $(WRF_MODULES)

wrf_pmc_trans.o: wrf_pmc_trans_aero.o wrf_pmc_dep.o $(WRF_MODULES)

wrf_pmc_init.o: wrf_pmc_driver.o $(WRF_MODULES)

wrf_pmc_dep.o: wrf_pmc_dep_aero.o $(WRF_MODULES)

wrf_pmc_dep_aero.o:

mem.o:
