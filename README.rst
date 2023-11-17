WRF-PartMC
----------

WRF-PartMC is a coupled model that combines the strengths of the Weather Research and Forecast (WRF) model
(state-of-the-art mesoscale numerical weather prediction) with the strengths of the particle-resolved aerosol
model PartMC (state-of-the-art aerosol dynamics and chemistry capable of resolving complex aerosol composition).
This enables the use of high detailed aerosol representation at the regional scale.

Acquiring and building WRF-PartMC
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To clone the WRF-PartMC repository:

   git clone https://github.com/open-atmos/wrf-partmc.git

To get the PartMC submodule:

   git submodule update --init partmc

And if you have access to MOSAIC:

   git submodule update --init mosaic

In order to compile with MOSAIC, set the environmental variable:

   export MOSAIC=1
   export WRF_CHEM=1

such that MOSAIC is included in the WRF-PartMC build process.
To configure and compile, follow the process for WRF which is:

  cd WRFV3
  
Configure by:
   
   ./configure
   
and select the options as you would in WRF. To compile WRF:

   ./compile em_real

Description of the directories found in this repo
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* WPS/ is the WRF preprocessor that creates the inputs for WRF (version 3.9.1)
* WRFV3/ is the meteorology model and is a modified version of 3.9.1. NOTE: WPS and WRF should match in version number.
* partmc/ is a submodule that points to a specific version that is for WRF-PartMC.
* mozbc/ takes wrfinput and wrfbdy from ``./real.exe`` and populates these files with values from the MOZART output.
* emissions/ contains the program to create WRF-PartMC emissions from SMOKE source-apportioned data
* boundary_and_initial_conditions/ contains programs to create initial and boundary conditions from the wrfinput and wrfbdy files using MOZART global model output.
* interface/ is all the code required to interface WRF with PartMC
* mosaic/ is the chemistry code which is available by request.
