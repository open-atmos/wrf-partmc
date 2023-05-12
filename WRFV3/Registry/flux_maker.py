# - species naming convention from aero_data and gas_data
# - WRF will loop over all variables in each group by incrementing through list
# of dname (I think) or whatever variable denotes position in the parent array
#    - WRF will probably make that variable lower case so we must adjust it 
#    - there will be an exception made for SOA. the mass code will be adjusted to
#      take SOA as a variable despite it not being a species and do the lumping

reg_format = '{:9} {:6} {:16} {:8} {:12} {:3} {:8} {:8} "{:16}" "{:16}" "{:16}" \n'

hist = "h3"

f = open('registry.partmc_trans','w')

var = "state"
type = "real"
dim = "ikjf"
species = ["u_flux","v_flux", "w_flux"]
stagger = ["X","Y","Z"]
units = [" "," ", " "]

package = []

for i_spec in range(len(species)):
 package = []
 array_name = "%s" %(species[i_spec])
 string = reg_format.format(var,type, "-",
       dim, array_name, '-', stagger[i_spec], '-','-', '-','-')
 f.write(string)
 for ii in range(1,41):
   var_name = "%s_a%03i" %(species[i_spec], ii)
   package.append(var_name)
   name = "%s" %var_name
   description="%s, class %03i" %(species[i_spec], ii)
   string = reg_format.format(var, type, var_name,
       dim, array_name, 1, stagger[i_spec], hist, name, description, units[i_spec])
   f.write(string)
 string = "package partmc_trans_%s chem_opt==777 - %s:%s" % (array_name,array_name,package[0])
 for ii in range(1,len(package)):
    string = "%s,%s" %(string,package[ii])
 string = "%s \n" %string
 f.write(string)

f.close()
