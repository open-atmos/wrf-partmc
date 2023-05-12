import os
import numpy as np

L = 40*2000.0
dt_standard = 10
dx = 10
nxs = [21,41,81]
n_parts = [100.0,250.0,500.0,1000.0,10000.0,100000.0]
for nx in nxs:
  for n_part in n_parts:
    filename_in = "namelist.input"
    filename_out = "namelist.input.%03i.%08i" % (nx, n_part)
    print("filename_out: %s" % filename_out)
    f_in = open(filename_in, 'r')
    f_out = open(filename_out, 'w')
    dx = L / (nx-1)
    dt = dt_standard * (40.0 / (nx-1))
    for line in f_in:
        line = line.replace('%%DT%%', '%3i' % np.floor(dt))
        if (dt == np.floor(dt)):
            line = line.replace('%%FRACT_NUM%%', '%3i' % 0)
            line = line.replace('%%FRACT_DEN%%', '%3i' % 1)
        else:
            decimal = dt % 1
            val = (str(decimal).split('.')[1])
            denom = np.ceil(np.log10(float(val)))
            line = line.replace('%%FRACT_NUM%%', '%3i' % int(val))
            line = line.replace('%%FRACT_DEN%%', '%5i' % 10**denom)
        line = line.replace('%%AUXHIST2%%', "'aerosols_%03i_%08i_d<domain>_<date>'" % (nx,n_part))
        line = line.replace('%%NX%%', '%3i' % nx)
        line = line.replace('%%NY%%', '%3i' % nx)
        line = line.replace('%%DX%%', '%3i' % dx)
        line = line.replace('%%DY%%', '%3i' % dx)    
        line = line.replace('%%N_PART%%', '%12.1f' % n_part)
        f_out.write(line)

    f_in.close()
    f_out.close()
