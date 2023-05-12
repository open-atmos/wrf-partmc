#!/bin/csh

echo Setting up quarter_ss case by linking data files into this directory

echo linking to RRTM_DATA in ../../run directory
ln -sf ../../run/RRTM_DATA .
ln -sf ../../run/LANDUSE.TBL .
ln -sf ../../run/VEGPARM.TBL .
ln -sf ../../run/SOILPARM.TBL .
ln -sf ../../run/GENPARM.TBL .
ln -sf ../../run/RRTM_DATA .
echo done
