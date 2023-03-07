#######################################################
#
#
#
#######################################################

import sys
import os
import numpy as np
import struct
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import def_HFM_homo


verbose=1.         # could be zero for less output
order=2.           # minimal stencil order
xorg=0.;zorg=0.    # square grid 
hcarre=0.01        # grid stepping
nx=100.;nz=100.    # grid dimensions 2D
slowness=1.        # cost
xsrc=0.5;zsrc=0.5  # seed
tsrc=0.            # time at seed
radius=-1.         # automatic factorization
export_time=1.     # output time field

field=def_HFM_homo.hfm_py(verbose,order,xorg,zorg,hcarre,nx,nz,slowness,xsrc,zsrc,tsrc,radius,export_time)

time=np.asarray(field)

nxc=int(nx)
nzc=int(nz)
time=time.reshape(nxc,nzc)
plt.imshow(time,origin='upper',aspect='equal',interpolation='bilinear',resample='true',cmap=cm.seismic,clip_on=True)
plt.show()
