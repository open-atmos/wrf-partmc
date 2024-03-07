WRF-PartMC
----------

WRF-PartMC is a coupled model that combines the strengths of the Weather Research
and Forecast (WRF) model (state-of-the-art mesoscale numerical weather prediction) 
with the strengths of the particle-resolved aerosol model PartMC
(state-of-the-art aerosol dynamics and chemistry capable of resolving complex aerosol composition).
This enables the use of high detailed aerosol representation at the regional scale.

Acquiring and building WRF-PartMC
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#. To clone the WRF-PartMC repository:

      git clone https://github.com/open-atmos/wrf-partmc.git

#. Change into the clone repository:

      cd wrf-partmc

#. To checkout the PartMC submodule:

      git submodule update --init partmc

#. And if you have access to MOSAIC, you may check it out it by:

      git submodule update --init mosaic

#. In order to compile with MOSAIC support, set the environmental variable:

      export MOSAIC=1

   This flag will signal to the WRF-PartMC build process to attempt to build MOSAIC.

#. WRF-PartMC requires WRF-Chem to be built so set the environmental flag (``WRF_CHEM``)
   for WRF-Chem to be included by:

      export WRF_CHEM=1

#. To configure and compile, follow the process for a standard WRF installation which is:

      cd WRFV3

#. Set the netCDF path as done for WRF, typically can be done by:

      export NETCDF=$(nc-config --prefix)
  
#. You should be ready to configure the model by:

      ./configure
   
   and select the compiler options as you would do in WRF. WRF-PartMC does not
   currently support domain nesting. You can still build WRF with nesting enabled but
   you'll be limited to one domain when enabling PartMC.

#. To compile WRF-PartMC,

     ./compile em_real

Other choices may include em_rotational, em_les and em_scm_xy. See individual directories in ``test/`` for any specific instructions.

Description of the directories found in this repo
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``WPS/`` is the WRF preprocessor that creates the inputs for WRF (version 3.9.1)
* ``WRFV3/`` is the meteorology model and is a modified version of 3.9.1. NOTE:
  WPS and WRF should match in version number.
* ``partmc/`` is a submodule that points to a specific version that is for WRF-PartMC.
* ``mozbc/`` takes wrfinput and wrfbdy from ``./real.exe`` and populates these
  files with values from the MOZART output.
* ``emissions/`` contains the program to create WRF-PartMC emissions from SMOKE
  source-apportioned data
* ``boundary_and_initial_conditions/`` contains programs to create initial and
  boundary conditions from the ``wrfinput`` and ``wrfbdy`` files using MOZART global model output.
* ``interface/`` is all the code required to interface WRF with PartMC
* ``mosaic/`` is the chemistry code which is available by request.
