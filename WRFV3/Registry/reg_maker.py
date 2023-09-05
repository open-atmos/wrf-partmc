# - species naming convention from aero_data and gas_data
# - WRF will loop over all variables in each group by incrementing through list
# of dname (I think) or whatever variable denotes position in the parent array
#    - WRF will probably make that variable lower case so we must adjust it 
#    - there will be an exception made for SOA. the mass code will be adjusted to
#      take SOA as a variable despite it not being a species and do the lumping

reg_format = '{:9} {:6} {:16} {:8} {:12} {:3} {:8} {:8} "{:16}" "{:16}" "{:16}"\n'

hist = "h2"
hist_3 = "h3"
ss = ["001","003","006","010"]

f = open('registry.partmc_process','w')

# Add any dimspec variables
string = "dimspec    bins    -     namelist=num_bins               c          n_bins \n"
f.write(string)
string = "dimspec    edges   -     namelist=num_edges             c          n_edges \n"
f.write(string)

# 1D array that are constant everywhere
var = "state"
type = "real"
dim = "{bins}"
var_name = "bin_centers"
name = "bin_centers"
description = ""
units = ""
string = reg_format.format(var, type, var_name,
       dim, 'pmc', 1, '-', hist_3, name, description, units)
f.write(string)

var = "state"
type = "real"
dim = "{edges}"
var_name = "bin_edges"
name = "bin_edges"
description = ""
units = ""
string = reg_format.format(var,type, var_name,
       dim, 'pmc', 1, '-', hist_3, name, description, units)
f.write(string)

# 3d variables that are really arrays
var = "state"
type = "real"
dim = "ikjf"
species = ["num","mass"]
units = ["#","mass"]

package = []

for i_spec in range(len(species)):
 package = []
 array_name = "%s_dist" %(species[i_spec])
 string = reg_format.format(var,type, "-",
       dim, array_name, '-', '-', '-','-', '-','-')
 f.write(string)
 for ii in range(1,101):
   var_name = "%s_a%03i" %(species[i_spec], ii)
   package.append(var_name)
   name = "%s" %var_name
   description="%s, aerosol bin %03i" %(species[i_spec], ii)
   string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist_3, name, description, units[i_spec])
   f.write(string)
 string = "package partmc_process_%s chem_opt==777 - %s:%s" % (array_name,array_name,package[0])
 for ii in range(1,len(package)):
    string = "%s,%s" %(string,package[ii])
 string = "%s \n" %string
 f.write(string)

# Bulk variables
dim = "ikj"
units = ""
array_name = "pmc"
var_name = "density_dry_air"
name = "density_dry_air"
description = "density of dry air"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "temperature"
name = "temperature"
description = "temperature"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "rel_humid"
name = "rel_humid"
description = "rel_humid"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

# Bulk variables
dim = "ikj"
units = ""
array_name = "pmc"
var_name = "tot_mass_conc"
name = "tot_mass_conc"
description = "total mass concentration"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "tot_num_conc"
name = "tot_num_conc"
description = "total number concentration"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "tot_wet_num_conc"
name = "tot_wet_num_conc"
description = "total wet number concentration"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "tot_hydrophobic_mass_conc"
name = "tot_hydrophobic_mass_conc"
description = "total hydrophobic mass concentration"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "tot_hydrophylic_mass_conc"
name = "tot_hydrophylic_mass_conc"
description = "total hydrophylic mass concentration"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

### PM totals
dim = "ikj"
units = ""
array_name = "pmc"
var_name = "PM1_mass_conc"
name = "PM1_mass_conc"
description = "total mass concentration less than 1 micron"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "PM25_mass_conc"
name = "PM25_mass_conc"
description = "total mass concentration less than 2.5 micron"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "PM10_mass_conc"
name = "PM10_mass_conc"
description = "total mass concentration less than 10 micron"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

### Optical properties
dim = "ikj"
units = ""
array_name = "pmc"
var_name = "ext_aer_550"
name = "ext_aer_550"
description = "extinction coefficient at 550 nm for particle-resolved"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "ext_aer_550_internal"
name = "ext_aer_550_internal"
description = "extinction coefficient at 550 nm for internal mixed"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "ext_aer_550_external"
name = "ext_aer_550_external"
description = "extinction coefficient at 550 nm for external mixed"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "scat_aer_550"
name = "scat_aer_550"
description = "scattering coefficient at 550 nm"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "scat_aer_550_internal"
name = "scat_aer_550_internal"
description = "scattering coefficient at 550 nm for internal mixed"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "scat_aer_550_external"
name = "scat_aer_550_external"
description = "scattering coefficient at 550 nm for external mixed"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

mode = ["a1","a2","a3"]
vars = ["num_conc","mass_conc"]
label = ["number concentration", "mass concentration"]
for i_var in range(len(vars)):
   for j_var in range(len(mode)):
      dim = "ikj"
      units = ""
      array_name = "pmc"
      var_name = "%s_%s" %(vars[i_var],mode[j_var])
      name = var_name
      description = "%s for mode %s" \
           %(label[i_var], mode[j_var])
      string = reg_format.format(var, type, var_name,
             dim, array_name, 1, '-', hist, name, description, units)
      f.write(string)    


vars = ["scat_aer_550","ext_aer_550"]
label = ["scattering", "extinction"]
assumption = ["pr","internal","external"]
mode = ["a1","a2","a3"]
for i_var in range(len(vars)):
   for j_var in range(len(mode)):
    for k_var in range(len(assumption)):
      dim = "ikj"
      units = ""
      array_name = "pmc"
      var_name = "%s_%s_%s" %(vars[i_var],assumption[k_var],mode[j_var])
      name = var_name
      description = "%s coefficient at 550 nm for %s mode %s" \
           %(label[i_var], assumption[k_var], mode[j_var])
      string = reg_format.format(var, type, var_name,
             dim, array_name, 1, '-', hist, name, description, units)
      f.write(string)

dim = "ikjf"
units = "m^-3"
for i_assumption in range(len(assumption)):
   array_name = "pmc_ccn_conc_modes_%s" %(assumption[i_assumption])
   string = reg_format.format(var,type, "-",
       dim, array_name, '-', '-', '-','-', '-','-')
   f.write(string)
   package = []
   for i_spec in range(len(ss)):
      for i_mode in range(len(mode)):
         var_name = "ccn_%s_%s_%s" %(assumption[i_assumption],ss[i_spec],mode[i_mode])
         package.append(var_name)
         name = "%s" %var_name
         description="%s, total mass" %(ss[i_spec])
         string = reg_format.format(var, type, var_name,
             dim, array_name, 1, '-', hist, name, description, units)
         f.write(string)

   string = "package partmc_process_cloud_modes_%s chem_opt==777 - pmc_ccn_conc_modes_%s:%s" % (
         assumption[i_assumption],assumption[i_assumption],package[0])
   for ii in range(1,len(package)):
      string = "%s,%s" %(string,package[ii])
   string = "%s \n" %string
   f.write(string)

#
dim = "ikj"
units = ""
array_name = "pmc"
var_name = "n_parts"
name = "n_parts"
description = "total computational particles"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "cell_vol"
name = "cell_vol"
description = "total volume of grid cell"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "n_components"
name = "n_components"
description = "total number of aerosol components"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "tot_num_conc_coagulated"
name = "tot_num_conc_coagulated"
description = "total number concentration of particles with many components"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "tot_bc_num_conc"
name = "tot_bc_num_conc"
description = "total number concentration of bc containing particles"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "tot_bc_num_conc_aged"
name = "tot_bc_num_conc_aged"
description = "total number concentration of aged bc containing particles"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "tot_coagulation_num_conc"
name = "tot_coagulation_num_conc"
description = "total number concentration of particles that coagulated within grid cell"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "d_alpha"
name = "d_alpha"
description = "alpha diversity"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "d_gamma"
name = "d_gamma"
description = "gamma diversity"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "chi"
name = "chi"
description = "mixing state parameter"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "d_alpha_ccn"
name = "d_alpha_ccn"
description = "d alpha ccn"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "d_gamma_ccn"
name = "d_gamma_ccn"
description = "d gamma ccn"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "chi_ccn"
name = "chi_ccn"
description = "mixing state parameter ccn"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "d_alpha_opt"
name = "d_alpha_opt"
description = "d alpha opt"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "d_gamma_opt"
name = "d_gamma_opt"
description = "d gamma opt"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "chi_opt"
name = "chi_opt"
description = "mixing state parameter opt"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

### submicron mixing state
dim = "ikj"
units = ""
array_name = "pmc"
var_name = "d_alpha_submicron"
name = "d_alpha_submicron"
description = "alpha diversity"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "d_gamma_submicron"
name = "d_gamma_submicron"
description = "gamma diversity"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "chi_submicron"
name = "chi_submicron"
description = "mixing state parameter"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "d_alpha_ccn_submicron"
name = "d_alpha_ccn_submicron"
description = "d alpha ccn"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "d_gamma_ccn_submicron"
name = "d_gamma_ccn_submicron"
description = "d gamma ccn"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "chi_ccn_submicron"
name = "chi_ccn_submicron"
description = "mixing state parameter ccn"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "d_alpha_opt_submicron"
name = "d_alpha_opt_submicron"
description = "d alpha opt"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "d_gamma_opt_submicron"
name = "d_gamma_opt_submicron"
description = "d gamma opt"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

dim = "ikj"
units = ""
array_name = "pmc"
var_name = "chi_opt_submicron"
name = "chi_opt_submicron"
description = "mixing state parameter opt"
string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
f.write(string)

# These should be 4D arrays eventually to incorporate a bunch of bins
chi = ['species','ccn','opt']
ranges = ['aitken','accumulation','coarse']
mode = ["a1","a2","a3"]
for i in range(3):
 for j in range(3):
   dim = "ikj"
   units = ""
   array_name = "pmc"
   var_name = "d_alpha_%s_%s" %(chi[j],mode[i])
   name = var_name 
   description = "d alpha %s in size range" %(chi[j])
   string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
   f.write(string)
   dim = "ikj"
   units = ""
   array_name = "pmc"
   var_name = "d_gamma_%s_%s" %(chi[j],mode[i])
   name = var_name 
   description = "d gamma %s in size range" %(chi[j])
   string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
   f.write(string)
   dim = "ikj"
   units = ""
   array_name = "pmc"
   var_name = "chi_%s_%s" %(chi[j],mode[i])
   name = var_name 
   description = "mixing state parameter %s in size range" %(chi[j])
   string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
   f.write(string)

# 3d variables that are really grouped scalars
# Aerosol species
species = ['SO4','NO3','Cl','NH4','MSA','ARO1','ARO2','ALK1','OLE1','API1',
     'API2','LIM1','LIM2','CO3','Na','Ca','OIN','OC','BC','H2O']
dim = "ikjf"
units = "total mass"
array_name = "aero_mass"
string = reg_format.format(var,type, "-",
       dim, array_name, '-', '-', '-','-', '-','-')
f.write(string)
package = []
for i_spec in range(len(species)):
   var_name = "pmc_%s" %(species[i_spec])
   package.append(var_name)
   name = "%s" %var_name
   description="%s, total mass" %(species[i_spec])
   string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
   f.write(string)

string = "package partmc_process_aero chem_opt==777 - aero_mass:%s" % package[0]
for ii in range(1,len(package)):
  string = "%s,%s" %(string,package[ii])
string = "%s\n" %string
f.write(string)

# 3d variables in a 4D array
dim = "ikjf"
units = "m^-3"
array_name = "pmc_ccn_conc"
string = reg_format.format(var,type, "-",
       dim, array_name, '-', '-', '-','-', '-','-')
f.write(string)
package = []
for i_spec in range(len(ss)):
   var_name = "ccn_%s" %(ss[i_spec])
   package.append(var_name)
   name = "%s" %var_name
   description="%s, total mass" %(ss[i_spec])
   string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
   f.write(string)

string = "package partmc_process_cloud chem_opt==777 - pmc_ccn_conc:%s" % package[0]
for ii in range(1,len(package)):
  string = "%s,%s" %(string,package[ii])
string = "%s\n" %string
f.write(string)

dim = "ikjf"
units = "m^-3"
array_name = "pmc_ccn_conc_internal"
string = reg_format.format(var,type, "-",
       dim, array_name, '-', '-', '-','-', '-','-')
f.write(string)
package = []
for i_spec in range(len(ss)):
   var_name = "ccn_internal_%s" %(ss[i_spec])
   package.append(var_name)
   name = "%s" %var_name
   description="%s, total mass" %(ss[i_spec])
   string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
   f.write(string)

string = "package partmc_process_cloud_internal chem_opt==777 - pmc_ccn_conc_internal:%s" % package[0]
for ii in range(1,len(package)):
  string = "%s,%s" %(string,package[ii])
string = "%s\n" %string
f.write(string)

dim = "ikjf"
units = "m^-3"
array_name = "pmc_ccn_conc_external"
string = reg_format.format(var,type, "-",
       dim, array_name, '-', '-', '-','-', '-','-')
f.write(string)
package = []
for i_spec in range(len(ss)):
   var_name = "ccn_external_%s" %(ss[i_spec])
   package.append(var_name)
   name = "%s" %var_name
   description="%s, total mass" %(ss[i_spec])
   string = reg_format.format(var, type, var_name,
       dim, array_name, 1, '-', hist, name, description, units)
   f.write(string)

string = "package partmc_process_cloud_external chem_opt==777 - pmc_ccn_conc_external:%s" % package[0]
for ii in range(1,len(package)):
  string = "%s,%s" %(string,package[ii])
string = "%s\n" %string
f.write(string)


# Gas species
#species = ['so2','o3','co','no','no2']
#dim = "ikjf"
#units = "ppb"
#array_name = "gas_mixrat"
#string = reg_format.format(var,type, "-",
#       dim, array_name, '-', '-', '-','-', '-','-')
#f.write(string)
#package = []
#for i_spec in range(len(species)):
#   var_name = "pmc_%s" %(species[i_spec])
#   package.append(var_name)
#   name = "%s" %var_name
#   description="%s, mixing ratio" %(species[i_spec])
#   string = reg_format.format(var, type, var_name,
#       dim, array_name, 1, '-', hist, name, description, units)
#   f.write(string)
#
#string = "package partmc_process_gas chem_opt==10 - gas_mixrat:%s" % package[0]
#for ii in range(1,len(package)):
#  string = "%s,%s" %(string,package[ii])
#string = "%s \n" %string
#f.write(string)

# Number concentrations for sources
max_sources = 150
package = []
var = "state"
type = "real"
dim = "ikjf"
units = ["#"]
array_name = "num_conc_source"
string = reg_format.format(var,type, "-",
       dim, array_name, '-', '-', '-','-', '-','-')
f.write(string)
package = []
for i_source in range(max_sources):
 var_name = "%s_%03i" %("num_conc_source", i_source)
 package.append(var_name)
 name = "%s" %var_name
 description="num conc source %03i" %(i_source)
 string = reg_format.format(var, type, var_name,
     dim, array_name, 1, '-', hist, name, description, units[0])
 f.write(string)
string = "package partmc_process_%s chem_opt==777 - %s:%s" % (array_name,array_name,package[0])
for ii in range(1,len(package)):
    string = "%s,%s" %(string,package[ii])
string = "%s\n" %string
f.write(string)

f.close()
